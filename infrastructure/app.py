"""Web front end for G-LoSA: upload two structures, get the alignment back.

Alignments are run on a background thread rather than inside the request.  A binding-site pair
finishes in well under a second, but the runtime scales with the square of the residue count and
a whole-protein pair can run for hours -- far past any sensible HTTP timeout.  The request
therefore only creates the job directory and returns; the status page polls.

All job state lives in files under `JOBS_DIR`, never in module globals, so it does not matter
which gunicorn worker serves the poll after another worker accepted the upload.
"""

from __future__ import annotations

import json
import shutil
import threading
import time
import traceback
import uuid
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
from pathlib import Path

from flask import (
    Flask,
    abort,
    jsonify,
    redirect,
    render_template,
    request,
    send_file,
    url_for,
)
from werkzeug.utils import secure_filename

import glosa_runner
from glosa_runner import GlosaError, align, probe_progress

HERE = Path(__file__).resolve().parent
PACKAGE_ROOT = HERE.parent
EXAMPLE_DIR = Path(glosa_runner.os.environ.get("GLOSA_EXAMPLES") or PACKAGE_ROOT / "example")
JOBS_DIR = Path(glosa_runner.os.environ.get("GLOSA_JOBS_DIR") or PACKAGE_ROOT / "runs" / "jobs")

MAX_UPLOAD_BYTES = int(glosa_runner.os.environ.get("GLOSA_MAX_UPLOAD_MB", "32")) * 1024 * 1024
JOB_RETENTION_HOURS = float(glosa_runner.os.environ.get("GLOSA_JOB_RETENTION_HOURS", "72"))

# One alignment at a time per worker.  The clique search is single-threaded but memory-hungry --
# a whole-protein product graph runs to gigabytes -- so letting several start at once is how the
# container gets itself OOM-killed rather than how it goes faster.
EXECUTOR = ThreadPoolExecutor(max_workers=1, thread_name_prefix="glosa")

# Files offered for download on the results page, in the order they are shown.
DOWNLOADS = [
    ("combined", "combined.pdb", "Overlay: reference + aligned structure in one file"),
    ("aligned", "ali_struct.pdb", "Mobile structure moved onto the reference"),
    ("transferred", "ali_struct_with.pdb", "Carried-along structure after the same transform"),
    ("matrix", "matrix.txt", "Rotation and translation"),
    ("convergence", "convergence.tsv", "Score at each search iteration"),
]

app = Flask(__name__)
app.config["MAX_CONTENT_LENGTH"] = MAX_UPLOAD_BYTES


# --------------------------------------------------------------------------- job state on disk


def job_dir(job_id: str) -> Path:
    """Resolve a job ID to its directory, refusing anything that is not a bare UUID.

    The ID reaches this from the URL, so treating it as a path fragment without checking would
    let `..` walk out of JOBS_DIR.
    """
    try:
        uuid.UUID(job_id)
    except ValueError:
        abort(404)
    return JOBS_DIR / job_id


def read_status(job_id: str) -> dict:
    path = job_dir(job_id) / "status.json"
    if not path.exists():
        abort(404)
    return decorate(json.loads(path.read_text()), job_dir(job_id))


def decorate(state: dict, directory: Path) -> dict:
    """Add the live fields -- elapsed time and current phase -- to a job's stored state.

    These are deliberately not written into status.json.  They change every second while a job
    runs, and persisting them would mean the worker thread and every polling request writing to
    the same file.  Computing them when asked costs a couple of stat() calls.
    """
    state = dict(state)
    if state.get("status") in ("queued", "running"):
        started = state.get("started") or state.get("submitted")
        if started:
            began = datetime.fromisoformat(started)
            state["elapsed"] = round((datetime.now(timezone.utc) - began).total_seconds(), 1)
        if state.get("status") == "running":
            state["progress"] = probe_progress(directory / "out")
    return state


def write_status(directory: Path, **fields) -> dict:
    path = directory / "status.json"
    state = json.loads(path.read_text()) if path.exists() else {}
    state.update(fields)
    path.write_text(json.dumps(state, indent=2, default=str))
    return state


def purge_old_jobs() -> None:
    """Drop job directories past the retention window.

    Each whole-protein run can leave tens of megabytes behind, so without this the volume fills.
    Called opportunistically on the index page rather than from a timer; the traffic here is
    occasional and a cron inside the container would be more machinery than the problem needs.
    """
    if not JOBS_DIR.exists():
        return
    cutoff = time.time() - JOB_RETENTION_HOURS * 3600
    for entry in JOBS_DIR.iterdir():
        try:
            if entry.is_dir() and entry.stat().st_mtime < cutoff:
                shutil.rmtree(entry, ignore_errors=True)
        except OSError:
            continue


# --------------------------------------------------------------------------- running a job


def save_upload(field: str, directory: Path) -> Path | None:
    """Persist an uploaded structure, or fall back to the named example if one was chosen.

    The form offers both an upload and an example dropdown for each slot so the service can be
    tried without hunting for a PDB file; the upload wins when both are set.
    """
    uploaded = request.files.get(field)
    if uploaded and uploaded.filename:
        name = secure_filename(uploaded.filename) or f"{field}.pdb"
        destination = directory / name
        uploaded.save(destination)
        return destination

    chosen = (request.form.get(f"{field}_example") or "").strip()
    if chosen:
        source = EXAMPLE_DIR / secure_filename(chosen)
        if not source.exists():
            raise GlosaError(f"unknown example structure {chosen!r}")
        destination = directory / source.name
        shutil.copyfile(source, destination)
        return destination

    return None


def cancel_marker(directory: Path) -> Path:
    """A file, not a flag in memory: the stop request and the worker thread can be in different
    gunicorn workers, and the marker is what lets the worker tell "I was stopped" apart from
    "I crashed" when the kill lands as a signal."""
    return directory / "cancelled"


def run_job(job_id: str, paths: dict[str, Path | None], extra: list[str]) -> None:
    """Execute one alignment and record the outcome in `status.json`."""
    directory = JOBS_DIR / job_id
    if cancel_marker(directory).exists():
        # Stopped while still queued, so the alignment never started.
        write_status(directory, status="cancelled", seconds=0)
        return
    write_status(directory, status="running", started=datetime.now(timezone.utc).isoformat())
    began = time.monotonic()
    try:
        result = align(
            structure1=paths["s1"],
            structure2=paths["s2"],
            transfer=paths["s2w"],
            features1=paths["s1cf"],
            features2=paths["s2cf"],
            outdir=directory / "out",
            extra_args=extra,
        )
        write_status(
            directory,
            status="done",
            ga_score=result.ga_score,
            matrix=result.matrix,
            chain_roles=result.chain_roles,
            stdout=result.stdout,
            outputs={k: v.name for k, v in result.outputs.items()},
            seconds=round(time.monotonic() - began, 2),
            finished=datetime.now(timezone.utc).isoformat(),
        )
    except GlosaError as exc:
        # A stop arrives as a SIGTERM to the scorer, which surfaces here as an ordinary
        # failure, so the marker is what distinguishes the two.
        if cancel_marker(directory).exists():
            # exc.partial is whatever the last checkpoint could be turned into, so stopping
            # early still yields a usable alignment rather than nothing.
            write_status(directory, status="cancelled",
                         seconds=round(time.monotonic() - began, 2),
                         **exc.partial)
        else:
            write_status(directory, status="failed", error=str(exc),
                         seconds=round(time.monotonic() - began, 2))
    except Exception:  # noqa: BLE001 - the traceback is the useful artefact for an unexpected fault
        write_status(directory, status="failed", error=traceback.format_exc()[-4000:],
                     seconds=round(time.monotonic() - began, 2))


# --------------------------------------------------------------------------- routes


@app.get("/")
def index():
    purge_old_jobs()
    examples = sorted(p.name for p in EXAMPLE_DIR.glob("*.pdb")
                      if not p.name.startswith("._") and "-cf" not in p.name)
    return render_template("index.html", examples=examples, recent=recent_jobs())


def recent_jobs(limit: int | None = 10) -> list[dict]:
    """Every job on disk, newest first, with the live fields filled in.

    Unfinished jobs are floated to the top regardless of age: the board exists to answer "is my
    submission still going", and a running job buried under newer finished ones would defeat it.
    """
    if not JOBS_DIR.exists():
        return []
    entries = []
    for directory in JOBS_DIR.iterdir():
        status_file = directory / "status.json"
        if not status_file.is_file():
            continue
        try:
            state = json.loads(status_file.read_text())
        except (OSError, json.JSONDecodeError):
            continue
        state = decorate(state, directory)
        state["id"] = directory.name
        state["mtime"] = status_file.stat().st_mtime
        entries.append(state)
    entries.sort(key=lambda s: (s.get("status") in ("queued", "running"), s["mtime"]),
                 reverse=True)
    return entries if limit is None else entries[:limit]


@app.post("/align")
def start_alignment():
    job_id = str(uuid.uuid4())
    directory = JOBS_DIR / job_id
    inputs = directory / "in"
    inputs.mkdir(parents=True, exist_ok=True)

    try:
        # s1cf/s2cf are overrides.  Left empty -- the normal case -- the feature files are
        # derived from the structures, so a PDB each is all the form actually needs.
        paths = {field: save_upload(field, inputs)
                 for field in ("s1", "s2", "s2w", "s1cf", "s2cf")}
    except GlosaError as exc:
        shutil.rmtree(directory, ignore_errors=True)
        return render_template("error.html", message=str(exc)), 400

    if paths["s1"] is None or paths["s2"] is None:
        shutil.rmtree(directory, ignore_errors=True)
        return render_template(
            "error.html",
            message="Two structures are required: pick an example or upload a PDB for each slot.",
        ), 400

    extra: list[str] = []
    for flag, field in [("-iter", "iter"), ("-itercf", "itercf"), ("-n", "norm")]:
        value = (request.form.get(field) or "").strip()
        if value:
            extra += [flag, value]

    write_status(
        directory,
        status="queued",
        id=job_id,
        s1=paths["s1"].name,
        s2=paths["s2"].name,
        s2w=paths["s2w"].name if paths["s2w"] else None,
        features="supplied" if (paths["s1cf"] or paths["s2cf"]) else "derived from the structures",
        extra=" ".join(extra),
        submitted=datetime.now(timezone.utc).isoformat(),
    )
    EXECUTOR.submit(run_job, job_id, paths, extra)
    return redirect(url_for("show_job", job_id=job_id))


@app.get("/jobs")
def jobs_board():
    """Every job and its current state, refreshed in place while anything is running."""
    jobs = recent_jobs(limit=None)
    return render_template(
        "jobs.html",
        jobs=jobs,
        active=sum(1 for j in jobs if j.get("status") in ("queued", "running")),
    )


@app.get("/jobs.json")
def jobs_feed():
    """Backs the board's auto-refresh, and is the endpoint to poll from a script."""
    return jsonify(jobs=recent_jobs(limit=None))


@app.get("/job/<job_id>")
def show_job(job_id: str):
    state = read_status(job_id)
    out = job_dir(job_id) / "out"
    downloads = [
        (label, name, blurb)
        for label, name, blurb in DOWNLOADS
        if (out / name).exists()
    ]
    # Laid out server-side and drawn as inline SVG, so the page pulls in no plotting library.
    plot = glosa_runner.convergence_plot(glosa_runner.read_convergence(out))
    return render_template("job.html", job=state, job_id=job_id, downloads=downloads, plot=plot)


@app.get("/job/<job_id>/status.json")
def job_status(job_id: str):
    """Polled by the status page while the alignment runs."""
    return jsonify(read_status(job_id))


@app.post("/job/<job_id>/stop")
def stop_job(job_id: str):
    """Stop a running or queued job.

    The marker is written before the kill, not after: the worker thread checks it to decide
    whether the scorer's death was a cancellation or a crash, and writing it second would leave
    a race in which the job is recorded as failed.
    """
    directory = job_dir(job_id)
    state = json.loads((directory / "status.json").read_text())
    if state.get("status") not in ("queued", "running"):
        return redirect(url_for("show_job", job_id=job_id))

    cancel_marker(directory).touch()
    killed = glosa_runner.kill_scorer(directory / "out")
    app.logger.warning("stop requested for job %s%s", job_id,
                       f" (killed scorer pid {killed})" if killed else " (nothing running yet)")
    # A queued job has no scorer to kill, and one between steps may have none either, so the
    # status is set here rather than relying on the worker's exception path alone.
    if not killed:
        write_status(directory, status="cancelled")
    return redirect(url_for("show_job", job_id=job_id))


@app.get("/job/<job_id>/convergence.svg")
def job_convergence(job_id: str):
    """The convergence plot as a standalone SVG fragment.

    A running job's page swaps this in as new iterations land.  Rendering it here rather than
    redrawing in JavaScript means the layout maths lives in one place.
    """
    out = job_dir(job_id) / "out"
    plot = glosa_runner.convergence_plot(glosa_runner.read_convergence(out))
    return render_template("_convergence.html", plot=plot)


@app.get("/job/<job_id>/preview.pdb")
def job_preview(job_id: str):
    """The alignment as it currently stands, built from the latest checkpoint.

    The mobile structure is transformed by the best matrix found so far and combined with the
    reference, which gives the same kind of overlay as the finished `combined.pdb` but while
    the search is still running.
    """
    directory = job_dir(job_id)
    out = directory / "out"
    checkpoint = glosa_runner.read_checkpoint(out)
    if not checkpoint or "matrix" not in checkpoint:
        abort(404)

    reference, mobile = out / "s1.pdb", out / "s2.pdb"
    if not (reference.is_file() and mobile.is_file()):
        abort(404)

    moved = out / "preview_mobile.pdb"
    moved.write_text(glosa_runner.transform_pdb(mobile.read_text(), checkpoint["matrix"]))
    preview = out / "preview.pdb"
    glosa_runner.write_combined(reference, moved, preview)
    return send_file(preview, mimetype="text/plain", max_age=0)


@app.get("/job/<job_id>/file/<name>")
def job_file(job_id: str, name: str):
    safe = secure_filename(name)
    path = job_dir(job_id) / "out" / safe
    if not path.is_file():
        abort(404)
    # as_attachment is off for .pdb so the viewer can fetch the same URL.
    return send_file(path, as_attachment=not safe.endswith(".pdb"), download_name=safe,
                     mimetype="text/plain")


@app.template_filter("duration")
def format_duration(seconds: float | None) -> str:
    """Seconds as something readable, since a job here spans sub-second to several hours."""
    if seconds is None:
        return ""
    seconds = int(seconds)
    if seconds < 60:
        return f"{seconds} s"
    if seconds < 3600:
        return f"{seconds // 60} min {seconds % 60:02d} s"
    return f"{seconds // 3600} h {(seconds % 3600) // 60:02d} min"


@app.get("/healthz")
def healthz():
    """Liveness probe that also proves the two binaries the service depends on are present."""
    ok = glosa_runner.GLOSA_BIN.exists() and glosa_runner.CLASSES_ACF.exists()
    return jsonify(ok=ok, glosa=str(glosa_runner.GLOSA_BIN), java=glosa_runner.JAVA_BIN), (
        200 if ok else 503
    )


@app.errorhandler(413)
def too_large(_):
    limit = MAX_UPLOAD_BYTES // (1024 * 1024)
    return render_template("error.html", message=f"Upload exceeds the {limit} MB limit."), 413


def reap_orphans() -> None:
    """Fail any job still marked queued or running when the process starts.

    The executor lives in the process, so a restart -- a redeploy, a crash, an OOM kill during a
    large alignment -- abandons whatever was in flight while leaving `status.json` saying
    "running".  The results page polls on that status, so without this the job's page would
    refresh forever against work that nothing is doing.
    """
    if not JOBS_DIR.exists():
        return
    for directory in JOBS_DIR.iterdir():
        status_file = directory / "status.json"
        if not status_file.is_file():
            continue
        try:
            state = json.loads(status_file.read_text())
        except (OSError, json.JSONDecodeError):
            continue
        if state.get("status") in ("queued", "running"):
            # The scorer does not die with the process that spawned it, and it is CPU-bound, so
            # an abandoned whole-protein job will otherwise sit on a core for hours against
            # work nobody is waiting for.  Its scratch has to go too: align()'s `finally` never
            # ran, so tens of megabytes of product graph are still on the volume.
            killed = glosa_runner.kill_scorer(directory / "out")
            for name in glosa_runner.SCRATCH_FILES:
                (directory / "out" / name).unlink(missing_ok=True)
            # warning, not info: the default level would swallow it, and a restart that
            # abandoned somebody's multi-hour alignment is exactly what you want in the log.
            app.logger.warning("reaped abandoned job %s%s", directory.name,
                               f" (killed scorer pid {killed})" if killed else "")
            write_status(
                directory,
                status="failed",
                error="The service restarted while this job was running, so it was abandoned. "
                      "Submit it again.",
            )


JOBS_DIR.mkdir(parents=True, exist_ok=True)
reap_orphans()
