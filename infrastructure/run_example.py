"""Smoke-test the pipeline on the bundled binding-site example, outside the container.

    make example                     # from the package root
    python infrastructure/run_example.py

Takes no arguments.  The bundled 1IA1 and 2BL9 binding sites are the cheap pair -- a few hundred
milliseconds -- so this checks the wiring without waiting on a whole-protein clique search.  It
should score 0.896194; G-LoSA is deterministic, so anything else means the binary was built
differently rather than that the inputs changed.
"""

from __future__ import annotations

from pathlib import Path

import glosa_runner
from glosa_runner import align

ROOT = Path(__file__).resolve().parent.parent
EXAMPLE = ROOT / "example"
OUTDIR = ROOT / "runs" / "1IA1_vs_2BL9"


def main() -> None:
    print(f"glosa binary : {glosa_runner.GLOSA_BIN}")
    print(f"java         : {glosa_runner.JAVA_BIN}")

    result = align(
        structure1=EXAMPLE / "1IA1-bs.pdb",
        structure2=EXAMPLE / "2BL9-bs.pdb",
        transfer=EXAMPLE / "2BL9-lig.pdb",
        outdir=OUTDIR,
    )

    print(f"\nGA-score     : {result.ga_score:.6f}")
    if result.matrix:
        print(f"translation  : {[round(v, 3) for v in result.matrix['translation']]}")
    print("\noutputs:")
    for label, path in sorted(result.outputs.items()):
        print(f"  {label:<12} {path.relative_to(ROOT)}  ({path.stat().st_size} bytes)")

    combined = result.outputs["combined"]
    chains = {line[21] for line in combined.read_text().splitlines()
              if line[:6].startswith(("ATOM", "HETATM"))}
    print(f"\ncombined.pdb chains: {sorted(chains)}")


if __name__ == "__main__":
    main()
