"""Run a G-LoSA alignment end to end and collect the results in one directory.

G-LoSA is three programs, not one: `AssignChemicalFeatures` (Java) turns each PDB into the
chemical-feature file that the scorer needs, and `glosa` (C++) consumes both structures plus both
feature files.  This module hides that sequence behind `align()` so the CLI and the web service
drive the identical pipeline.

Two properties of the `glosa` binary shape the design:

* It writes its outputs -- `ali_struct.pdb`, `matrix.txt`, `ali_struct_with.pdb` -- to the current
  working directory under fixed names, so concurrent runs in a shared directory would overwrite
  each other.  Every call therefore gets its own scratch directory.
* It also drops two intermediate files there, `pairs.rst` and `product_graph.rst`.  For a
  binding-site pair those are trivial, but for a whole-protein pair `product_graph.rst` reaches
  tens of gigabytes -- the maximum-clique product graph grows with the square of the residue
  count.  They are deleted as soon as the run finishes.
"""

from __future__ import annotations

import os
import re
import shutil
import signal
import subprocess
import tempfile
import time
from dataclasses import dataclass, field
from pathlib import Path

HERE = Path(__file__).resolve().parent
PACKAGE_ROOT = HERE.parent

# In the image these all sit under /opt/glosa; in a source checkout they sit beside the .cpp and
# .java files they were built from.  Env vars win so the container does not depend on the layout.
GLOSA_BIN = Path(os.environ.get("GLOSA_BIN") or PACKAGE_ROOT / "glosa")
CLASSES_ACF = Path(os.environ.get("GLOSA_CLASSES_ACF") or PACKAGE_ROOT / "classes" / "acf")
CLASSES_ASS = Path(os.environ.get("GLOSA_CLASSES_ASS") or PACKAGE_ROOT / "classes" / "ass")
JAVA_BIN = os.environ.get("JAVA_BIN") or shutil.which("java") or "java"

# A whole-protein pair can run for hours and spill a very large product graph, which is a denial
# of service on a shared endpoint.  Callers may raise it for deliberate offline runs.
DEFAULT_TIMEOUT = int(os.environ.get("GLOSA_TIMEOUT", "900"))

# Scratch files glosa leaves in its working directory.  product_graph.rst is the big one.
# checkpoint.rst is the progress file written by the local patch to glosa.cpp; it is redundant
# once matrix.txt exists, so it is cleaned up with the rest.
CHECKPOINT_FILE = "checkpoint.rst"
SCRATCH_FILES = ("pairs.rst", "product_graph.rst", CHECKPOINT_FILE, "checkpoint.tmp")

# The convergence history, appended by the same patch.  Deliberately not scratch: it is the only
# record of how the score got where it did, and it survives the run.
CONVERGENCE_FILE = "convergence.tsv"

_GA_SCORE = re.compile(r"GA-score:\s*([0-9.eE+-]+)")


class GlosaError(RuntimeError):
    """A step of the pipeline failed, or the input was not usable.

    `partial` carries whatever could be salvaged from the last checkpoint when a run is stopped
    part-way, so a cancelled job still hands back the best alignment it had reached.
    """

    def __init__(self, message: str, partial: dict | None = None):
        super().__init__(message)
        self.partial: dict = partial or {}


@dataclass
class AlignResult:
    ga_score: float
    stdout: str
    workdir: Path
    matrix: dict[str, list[float]] | None = None
    outputs: dict[str, Path] = field(default_factory=dict)
    # Which chains of combined.pdb hold the reference, the aligned copy and anything carried
    # along.  Not fixed letters: the reference keeps whatever chains it arrived with.
    chain_roles: dict[str, list[str]] = field(default_factory=dict)


def _run(cmd: list[str], cwd: Path, timeout: int,
         pid_file: Path | None = None) -> subprocess.CompletedProcess:
    """Run one step to completion, optionally publishing its PID while it lives.

    Popen rather than `subprocess.run` only so the PID can be written out: the scorer prints
    nothing until it is finished, so the only way to report on a job in flight is to look at the
    process and the files it is leaving behind, and both need the PID.
    """
    argv = [str(c) for c in cmd]
    popen = subprocess.Popen(argv, cwd=str(cwd), stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE, text=True)
    if pid_file is not None:
        pid_file.write_text(str(popen.pid))
    try:
        stdout, stderr = popen.communicate(timeout=timeout)
    except subprocess.TimeoutExpired as exc:
        popen.kill()
        popen.communicate()
        raise GlosaError(f"{Path(cmd[0]).name} exceeded the {timeout}s time limit") from exc
    finally:
        if pid_file is not None:
            pid_file.unlink(missing_ok=True)

    proc = subprocess.CompletedProcess(argv, popen.returncode, stdout, stderr)
    if proc.returncode != 0:
        detail = (proc.stderr or proc.stdout or "").strip()[-2000:]
        name = Path(cmd[0]).name
        if proc.returncode < 0:
            # subprocess reports a fatal signal as a negative return code.  glosa has no input
            # validation of its own and dies this way on anything it cannot parse, so the raw
            # number on its own tells the user nothing.
            signal_name = signal.Signals(-proc.returncode).name
            raise GlosaError(
                f"{name} was killed by {signal_name}. This usually means one of the structures "
                f"is not something it can score. {detail}".strip()
            )
        raise GlosaError(f"{name} exited {proc.returncode}: {detail}")
    return proc


# Where glosa publishes the PID of the running scorer, and how long a scratch file may go
# untouched before it counts as finished rather than still being written.
PID_FILE = "glosa.pid"
WRITING_WINDOW = 15.0


def _read_proc_stat(pid: int) -> dict | None:
    """CPU seconds and resident size for a live PID, straight from /proc.

    Reading /proc avoids shelling out to ps on every poll.  Returns None if the process has
    already exited, which is the normal race between a poll and a job finishing.
    """
    try:
        fields = (Path("/proc") / str(pid) / "stat").read_text().rsplit(") ", 1)[1].split()
        statm = (Path("/proc") / str(pid) / "statm").read_text().split()
    except (OSError, IndexError):
        return None
    ticks = os.sysconf("SC_CLK_TCK")
    page = os.sysconf("SC_PAGE_SIZE")
    # Fields 11/12 after the comm field are utime/stime; statm[1] is the resident page count.
    return {
        "cpu_seconds": round((int(fields[11]) + int(fields[12])) / ticks, 1),
        "rss_mb": round(int(statm[1]) * page / 1024 / 1024, 1),
    }


def read_checkpoint(workdir: Path) -> dict | None:
    """Read the best-so-far alignment that the patched scorer publishes each iteration.

    See the `writeCheckpoint` block in glosa.cpp.  Upstream writes nothing until the whole run
    is over, so this file is the only window into a search in progress; the matrix in it is
    u_best/t_best, which is complete and usable at any point.
    """
    try:
        lines = (workdir / CHECKPOINT_FILE).read_text().splitlines()
    except OSError:
        return None

    checkpoint: dict = {}
    rotation: list[list[float]] = []
    translation: list[float] = []
    for line in lines:
        key, _, rest = line.partition("\t")
        fields = rest.split()
        try:
            if key == "stage":
                checkpoint["stage"] = rest.strip()
            elif key == "iteration" and len(fields) == 2:
                checkpoint["iteration"], checkpoint["iterations"] = int(fields[0]), int(fields[1])
            elif key in ("score", "best_score"):
                checkpoint[key] = float(fields[0])
            elif key == "matrix" and len(fields) == 4:
                values = [float(f) for f in fields]
                translation.append(values[0])
                rotation.append(values[1:])
        except (ValueError, IndexError):
            continue  # a torn read; the next poll gets a whole file

    if len(translation) == 3:
        checkpoint["matrix"] = {"translation": translation, "rotation": rotation}
    return checkpoint or None


def read_convergence(workdir: Path) -> list[dict]:
    """The score history appended by the patched scorer, one row per search iteration.

    A whole run produces a dozen or so rows, so this is read in full rather than tailed.
    """
    try:
        lines = (workdir / CONVERGENCE_FILE).read_text().splitlines()
    except OSError:
        return []

    rows = []
    for line in lines[1:]:  # skip the header
        fields = line.split("\t")
        if len(fields) != 6:
            continue  # a torn append; the row is skipped rather than guessed at
        try:
            rows.append({
                "seconds": float(fields[0]),
                "stage": fields[1],
                "iteration": int(fields[2]),
                "total": int(fields[3]),
                "score": float(fields[4]),
                "best_score": float(fields[5]),
            })
        except ValueError:
            continue
    return rows


def convergence_plot(rows: list[dict], width: int = 700, height: int = 180) -> dict | None:
    """Lay out the convergence history as coordinates for an inline SVG.

    Returning geometry rather than drawing a chart keeps this dependency-free: the page needs no
    plotting library, which on a closed network is one less thing that can fail to load, and the
    whole plot is a polyline and a dozen circles.

    The x axis is the iteration index, not elapsed seconds.  Time is recorded too, but the
    points are wildly uneven along it -- the clique search can spend an hour between two
    consecutive rows -- so plotting against time collapses everything into the left edge.
    """
    if len(rows) < 2:
        return None

    pad_l, pad_r, pad_t, pad_b = 46, 10, 12, 26
    plot_w, plot_h = width - pad_l - pad_r, height - pad_t - pad_b

    ymax = max(max(r["best_score"] for r in rows), 0.001) * 1.08
    last = len(rows) - 1

    def x(i: int) -> float:
        return pad_l + (plot_w * i / last)

    def y(value: float) -> float:
        return pad_t + plot_h * (1 - value / ymax)

    return {
        "width": width,
        "height": height,
        "baseline": pad_t + plot_h,
        "left": pad_l,
        "right": width - pad_r,
        "ymax": ymax,
        # The running best, which only ever climbs.
        "best": " ".join(f"{x(i):.1f},{y(r['best_score']):.1f}" for i, r in enumerate(rows)),
        # Each individual attempt, most of which are far below the best.
        "attempts": [
            {"x": round(x(i), 1), "y": round(y(r["score"]), 1), "score": r["score"],
             "seconds": r["seconds"], "stage": r["stage"],
             "iteration": r["iteration"], "total": r["total"],
             "improved": i == 0 or r["best_score"] > rows[i - 1]["best_score"]}
            for i, r in enumerate(rows)
        ],
        "elapsed": rows[-1]["seconds"],
        "final": rows[-1]["best_score"],
    }


def transform_pdb(pdb_text: str, matrix: dict) -> str:
    """Apply a G-LoSA translation/rotation to every coordinate record.

    This is the same arithmetic matrix.txt documents, X = t[i] + u[i]·r, and it is what lets the
    front end draw the current alignment without waiting for the scorer to write ali_struct.pdb
    at the very end of the run.
    """
    t, u = matrix["translation"], matrix["rotation"]
    out = []
    for line in pdb_text.splitlines():
        if line[:6].startswith(("ATOM", "HETATM")) and len(line) >= 54:
            try:
                r = [float(line[30:38]), float(line[38:46]), float(line[46:54])]
            except ValueError:
                out.append(line)
                continue
            moved = [t[i] + sum(u[i][j] * r[j] for j in range(3)) for i in range(3)]
            line = line[:30] + "".join(f"{v:8.3f}" for v in moved) + line[54:]
        out.append(line)
    return "\n".join(out) + "\n"


def write_matrix(dest: Path, matrix: dict, note: str = "") -> None:
    """Write a matrix.txt in glosa's layout, so the same parser reads either file.

    Columns are one place wider than upstream's `%15.10f`, which is exactly wide enough that a
    value like -85.6423730452 fills the field and runs into its neighbour; the extra column
    guarantees the whitespace that `parse_matrix` splits on.
    """
    t, u = matrix["translation"], matrix["rotation"]
    lines = ["*** Rotation matrix to superpose a structure to reference structure ***"]
    if note:
        lines.append(note)
    lines.append("i         t(i)           u(i,1)        u(i,2)         u(i,3)")
    for i in range(3):
        lines.append(f"{i + 1:<4}" + "".join(f"{v:>16.10f}" for v in (t[i], *u[i])))
    lines += [
        "",
        "X = t[1] + u[1][1]*x + u[1][2]*y + u[1][3]*z",
        "Y = t[2] + u[2][1]*x + u[2][2]*y + u[2][3]*z",
        "Z = t[3] + u[3][1]*x + u[3][2]*y + u[3][3]*z",
    ]
    dest.write_text("\n".join(lines) + "\n")


def salvage_checkpoint(outdir: Path) -> dict:
    """Turn the last checkpoint into real outputs, for a run that did not reach the end.

    A stopped job would otherwise leave nothing at all, even though the search had already
    found an alignment worth keeping.  On the whole-protein pair the best score is reached
    within the first of three iterations and the remaining hours only confirm it, so the
    salvaged result is usually the same answer the full run would have produced.
    """
    checkpoint = read_checkpoint(outdir)
    reference, mobile = outdir / "s1.pdb", outdir / "s2.pdb"
    if not checkpoint or "matrix" not in checkpoint:
        return {}
    if not (reference.is_file() and mobile.is_file()):
        return {}

    aligned = outdir / "ali_struct.pdb"
    aligned.write_text(transform_pdb(mobile.read_text(), checkpoint["matrix"]))
    write_matrix(outdir / "matrix.txt", checkpoint["matrix"],
                 note="*** PARTIAL: best alignment when the run was stopped ***")
    roles = write_combined(reference, aligned, outdir / "combined.pdb")

    return {
        "ga_score": checkpoint.get("best_score"),
        "chain_roles": roles,
        "matrix": checkpoint["matrix"],
        "outputs": {"aligned": aligned.name,
                    "combined": "combined.pdb",
                    "matrix": "matrix.txt"},
        "iteration": checkpoint.get("iteration"),
        "iterations": checkpoint.get("iterations"),
    }


def kill_scorer(workdir: Path) -> int | None:
    """Kill the scorer recorded in `workdir`'s PID file, returning the PID if one was killed.

    glosa is spawned as a child of the web process but does not die with it: killing the server
    mid-job leaves the scorer running, and because it is CPU-bound it will sit on a core for
    hours against a job nobody is waiting for any more.

    The PID is checked against /proc before signalling.  A recorded PID can be stale by minutes
    on a restart, long enough for the kernel to have handed the number to something else, and
    the comparison against the binary's name is what stops this killing an unrelated process.
    """
    pid_file = workdir / PID_FILE
    try:
        pid = int(pid_file.read_text().strip())
    except (OSError, ValueError):
        return None

    try:
        running = (Path("/proc") / str(pid) / "comm").read_text().strip()
    except OSError:
        running = ""

    killed = None
    if running and running == GLOSA_BIN.name[:len(running)]:
        try:
            os.kill(pid, signal.SIGTERM)
            killed = pid
        except OSError:
            pass
    pid_file.unlink(missing_ok=True)
    return killed


def probe_progress(workdir: Path) -> dict:
    """Describe what a running alignment is currently doing.

    glosa writes nothing to stdout until it is completely finished, so there is no progress to
    read.  What it does do is leave a predictable trail in its working directory, and the phase
    can be recovered from which files exist and whether they are still being appended to:

        chemical features -> pairs.rst -> product_graph.rst -> clique search -> results

    The last transition is the one that matters.  Building the product graph is quick and its
    file grows visibly; the maximum-clique search that follows touches no file at all and is
    where essentially all the time goes -- nearly three hours of a three-hour run, against a
    product_graph.rst that stopped changing in the first two minutes.  A job whose scratch has
    gone quiet is therefore not stuck, and saying so is the whole point of this function.

    There is no percentage here on purpose.  The search is branch-and-bound over a product
    graph, so the remaining work is not knowable in advance and any bar would be a fiction.
    """
    now = time.time()

    def age(path: Path) -> float | None:
        try:
            return now - path.stat().st_mtime
        except OSError:
            return None

    pairs, graph = workdir / "pairs.rst", workdir / "product_graph.rst"
    graph_age, pairs_age = age(graph), age(pairs)

    if graph_age is not None:
        if graph_age < WRITING_WINDOW:
            phase, detail = "building the product graph", "Comparing the two structures."
        elif pairs_age is not None and pairs_age < WRITING_WINDOW:
            phase, detail = "writing results", "The search is done; collecting the output."
        else:
            phase = "searching for the maximum clique"
            detail = ("This is the long step and it reports nothing while it runs. It is "
                      "branch-and-bound, so the time left cannot be estimated -- a pair of "
                      "binding sites finishes instantly, two whole chains took three hours.")
    elif pairs_age is not None:
        phase, detail = "enumerating residue pairs", "Finding candidate correspondences."
    elif any(workdir.glob("*-cf.pdb")):
        phase, detail = "starting the scorer", "Chemical features are assigned."
    else:
        phase, detail = "assigning chemical features", "Reading the two structures."

    progress: dict = {"phase": phase, "detail": detail}

    # The checkpoint is more informative than the file-mtime guesswork above whenever it exists,
    # because the scorer says outright which stage and iteration it is on.  The guesswork stays
    # as the fallback for the window before the first iteration completes.
    checkpoint = read_checkpoint(workdir)
    if checkpoint:
        progress["best_score"] = checkpoint.get("best_score")
        progress["iteration"] = checkpoint.get("iteration")
        progress["iterations"] = checkpoint.get("iterations")
        progress["stage"] = checkpoint.get("stage")
        progress["has_preview"] = "matrix" in checkpoint

    # Just the count: the page uses it to decide whether the plot is worth re-fetching, so
    # sending the rows themselves on every poll would be wasted payload.
    rows = read_convergence(workdir)
    if rows:
        progress["points"] = len(rows)

    scratch = sum(p.stat().st_size for p in (pairs, graph) if p.exists())
    if scratch:
        progress["scratch_mb"] = round(scratch / 1024 / 1024, 1)

    try:
        stats = _read_proc_stat(int((workdir / PID_FILE).read_text().strip()))
    except (OSError, ValueError):
        stats = None
    if stats:
        progress.update(stats)

    return progress


def normalise_pdb(src: Path, dst: Path) -> Path:
    """Copy `src` to `dst` keeping only coordinate records, and guarantee a trailing TER.

    G-LoSA's parsers stop at the first `TER` and require one to be present; a file without it is
    read as truncated and yields an empty structure rather than an error.  Uploaded files also
    routinely carry CRLF line endings and multi-model NMR ensembles, both of which break the
    fixed-column parsing, so only the first model's ATOM/HETATM records survive this pass.
    """
    kept: list[str] = []
    with src.open("r", errors="replace") as handle:
        for line in handle:
            record = line[:6]
            if record.startswith(("ATOM", "HETATM")):
                kept.append(line.rstrip("\r\n"))
            elif record.startswith("ENDMDL"):
                break  # first model only

    if not kept:
        raise GlosaError(f"{src.name} contains no ATOM or HETATM records")

    dst.write_text("\n".join(kept) + "\nTER\n")
    return dst


def assign_chemical_features(pdb: Path, workdir: Path, timeout: int, label: str = "") -> Path:
    """Produce the `<stem>-cf.pdb` feature file that `glosa` needs for a structure.

    AssignChemicalFeatures writes beside its input and derives the name from it, so the caller
    gets the path back rather than choosing it.  `label` is the name to use when complaining
    about the input: the working copies are called s1.pdb and s2.pdb, which would mean nothing
    to whoever uploaded the file.
    """
    _run([JAVA_BIN, "-cp", str(CLASSES_ACF), "AssignChemicalFeatures", pdb.name], workdir, timeout)
    produced = workdir / f"{pdb.stem}-cf.pdb"
    if not produced.exists():
        raise GlosaError(f"no chemical-feature file was produced for {label or pdb.name}")
    _require_features(produced, label or pdb.name)
    return produced


def _require_features(feature_file: Path, source_name: str) -> None:
    """Reject a feature file with no points in it, before glosa gets a chance to crash on it.

    AssignChemicalFeatures only recognises the twenty standard residues, so a ligand-only or
    nucleic-acid PDB yields a file of three bytes rather than an error.  glosa does not check,
    and segfaults.  The structures that hit this are usually a ligand picked for the reference
    or mobile slot, which is what the message points at.
    """
    points = sum(1 for line in feature_file.read_text().splitlines()
                 if line[:6].startswith(("ATOM", "HETATM")))
    if points == 0:
        raise GlosaError(
            f"{source_name} yielded no chemical feature points, so there is nothing to score. "
            f"G-LoSA reads features from standard amino-acid residues only -- if this is a "
            f"ligand, put it in the 'carried along' slot rather than the reference or mobile one."
        )


def parse_matrix(path: Path) -> dict[str, list[float]] | None:
    """Pull the translation and rotation out of glosa's `matrix.txt`.

    The three numbered rows carry `t(i) u(i,1) u(i,2) u(i,3)`; everything else in the file is the
    header and the restatement of X/Y/Z as formulae.
    """
    if not path.exists():
        return None
    translation: list[float] = []
    rotation: list[list[float]] = []
    for line in path.read_text().splitlines():
        fields = line.split()
        if len(fields) == 5 and fields[0] in {"1", "2", "3"}:
            try:
                values = [float(f) for f in fields[1:]]
            except ValueError:
                continue
            translation.append(values[0])
            rotation.append(values[1:])
    if len(translation) != 3:
        return None
    return {"translation": translation, "rotation": rotation}


def _chain_ids(pdb_text: str) -> set[str]:
    return {line[21] for line in pdb_text.splitlines() if line[:6].startswith(("ATOM", "HETATM"))}


def _relabel_chains(pdb_text: str, taken: set[str]) -> tuple[str, dict[str, str]]:
    """Move every chain in `pdb_text` onto a chain ID that is not already in `taken`."""
    alphabet = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789"
    mapping: dict[str, str] = {}
    for original in sorted(_chain_ids(pdb_text)):
        free = next((c for c in alphabet if c not in taken and c not in mapping.values()), "Z")
        mapping[original] = free
        taken.add(free)

    out = []
    for line in pdb_text.splitlines():
        if line[:6].startswith(("ATOM", "HETATM")) and len(line) > 21:
            line = line[:21] + mapping.get(line[21], line[21]) + line[22:]
        out.append(line)
    return "\n".join(out), mapping


def write_combined(reference: Path, aligned: Path, dest: Path, ligand: Path | None = None) -> dict:
    """Write one PDB holding the reference and the aligned structure on distinct chains.

    A multi-model file would be the other option, but most viewers show only the first model
    unless told otherwise, so the point of the overlay would be lost.  Renaming the moved
    structure's chains instead means the file opens with both bodies visible and colourable by
    chain.
    """
    ref_text = reference.read_text()
    reference_chains = sorted(_chain_ids(ref_text))
    taken = set(reference_chains)

    aligned_text, aligned_map = _relabel_chains(aligned.read_text(), taken)

    ligand_text, ligand_map = "", {}
    if ligand is not None and ligand.exists():
        ligand_text, ligand_map = _relabel_chains(ligand.read_text(), taken)

    # The reference keeps the chains it arrived with, so which letters the other two end up on
    # depends on the input -- a reference on chain H leaves the aligned copy at A, not B.  The
    # mapping is recorded here and returned, rather than left for the reader to infer, because
    # the viewer colours by role and would otherwise have to guess.
    roles = {
        "reference": reference_chains,
        "aligned": sorted(aligned_map.values()),
        "transferred": sorted(ligand_map.values()),
    }

    parts = ["REMARK   1 G-LoSA overlay: reference + aligned mobile structure"]
    parts += [f"REMARK   1 {role + ' chains':<20}{','.join(chains) or '-'}"
              for role, chains in roles.items()]
    parts.append(ref_text.rstrip("\n"))
    parts.append(aligned_text.rstrip("\n"))
    if ligand_text:
        parts.append(ligand_text.rstrip("\n"))

    dest.write_text("\n".join(parts) + "\nEND\n")
    return roles


def _feature_file(structure: Path, supplied: Path | None, workdir: Path, timeout: int,
                  label: str) -> Path:
    """Use a caller-supplied feature file if there is one, otherwise derive it from the PDB.

    The upstream instructions have the user run AssignChemicalFeatures by hand and pass the
    result to `-s1cf`/`-s2cf`, but nothing in it needs a human decision -- the features follow
    from the residue and atom names alone.  Deriving it here means a structure is the only thing
    anyone has to provide; `supplied` exists for the case where someone has hand-edited a feature
    set and wants that used verbatim.
    """
    if supplied is None:
        return assign_chemical_features(structure, workdir, timeout, label)
    destination = workdir / f"{structure.stem}-cf.pdb"
    shutil.copyfile(supplied, destination)
    _require_features(destination, supplied.name)
    return destination


def align(
    structure1: Path,
    structure2: Path,
    outdir: Path,
    transfer: Path | None = None,
    features1: Path | None = None,
    features2: Path | None = None,
    extra_args: list[str] | None = None,
    timeout: int = DEFAULT_TIMEOUT,
    keep_scratch: bool = False,
) -> AlignResult:
    """Align `structure2` onto `structure1` and collect every artefact in `outdir`.

    `transfer` is glosa's `-s2w`: a structure that is not scored but is moved by the resulting
    matrix, which is how a ligand is carried along with its binding site.

    `features1` and `features2` are optional.  Left unset -- the normal case -- the chemical
    feature files are computed from the two structures.
    """
    if not GLOSA_BIN.exists():
        raise GlosaError(f"the glosa binary is missing at {GLOSA_BIN}; compile it first")

    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    # glosa names its outputs after nothing at all, so the inputs are copied in under stable
    # names and the whole run happens here.
    s1 = normalise_pdb(structure1, outdir / "s1.pdb")
    s2 = normalise_pdb(structure2, outdir / "s2.pdb")
    s2w = normalise_pdb(transfer, outdir / "s2w.pdb") if transfer is not None else None

    cf1 = _feature_file(s1, features1, outdir, timeout, structure1.name)
    cf2 = _feature_file(s2, features2, outdir, timeout, structure2.name)

    cmd = [GLOSA_BIN, "-s1", s1.name, "-s1cf", cf1.name, "-s2", s2.name, "-s2cf", cf2.name]
    if s2w is not None:
        cmd += ["-s2w", s2w.name]
    cmd += extra_args or []

    try:
        proc = _run(cmd, outdir, timeout, pid_file=outdir / PID_FILE)
    except GlosaError as exc:
        # Salvage before the `finally` below removes the checkpoint this reads.
        exc.partial = salvage_checkpoint(outdir)
        raise
    finally:
        if not keep_scratch:
            for name in SCRATCH_FILES:
                (outdir / name).unlink(missing_ok=True)

    match = _GA_SCORE.search(proc.stdout)
    if not match:
        raise GlosaError(f"glosa produced no GA-score; output was:\n{proc.stdout.strip()[-2000:]}")

    aligned = outdir / "ali_struct.pdb"
    transferred = outdir / "ali_struct_with.pdb"
    outputs: dict[str, Path] = {}
    chain_roles: dict[str, list[str]] = {}
    if aligned.exists():
        outputs["aligned"] = aligned
        combined = outdir / "combined.pdb"
        chain_roles = write_combined(
            s1, aligned, combined, transferred if transferred.exists() else None
        )
        outputs["combined"] = combined
    if transferred.exists():
        outputs["transferred"] = transferred
    matrix_path = outdir / "matrix.txt"
    if matrix_path.exists():
        outputs["matrix"] = matrix_path

    return AlignResult(
        ga_score=float(match.group(1)),
        stdout=proc.stdout,
        workdir=outdir,
        matrix=parse_matrix(matrix_path),
        outputs=outputs,
        chain_roles=chain_roles,
    )


def align_temp(structure1: Path, structure2: Path, **kwargs) -> AlignResult:
    """`align()` into a fresh temporary directory; the caller owns and must remove it."""
    return align(structure1, structure2, Path(tempfile.mkdtemp(prefix="glosa-")), **kwargs)
