"""Prove the three G-LoSA programs work together inside the image.

`make verify` only shows that gunicorn answers HTTP.  This runs a real alignment of the bundled
1IA1 and 2BL9 binding sites, which exercises the jlink runtime, the compiled feature assigner and
the C++ scorer in sequence, and checks the score against the value the upstream package produces.

Takes no arguments so it can be driven straight from `docker exec`.
"""

from __future__ import annotations

import sys
import tempfile
from pathlib import Path

import glosa_runner
from glosa_runner import GlosaError, align

EXAMPLES = Path(glosa_runner.os.environ.get("GLOSA_EXAMPLES", "/opt/glosa/example"))

# G-LoSA is deterministic, so the bundled pair always scores the same; a drift here means the
# binary was built differently, not that the inputs changed.
EXPECTED = 0.896194
TOLERANCE = 1e-4


def main() -> int:
    print(f"glosa : {glosa_runner.GLOSA_BIN}")
    print(f"java  : {glosa_runner.JAVA_BIN}")

    with tempfile.TemporaryDirectory(prefix="glosa-smoke-") as tmp:
        try:
            result = align(
                structure1=EXAMPLES / "1IA1-bs.pdb",
                structure2=EXAMPLES / "2BL9-bs.pdb",
                transfer=EXAMPLES / "2BL9-lig.pdb",
                outdir=Path(tmp),
            )
        except GlosaError as exc:
            print(f"FAIL: {exc}")
            return 1

        print(f"GA-score: {result.ga_score:.6f} (expected {EXPECTED:.6f})")
        for label, path in sorted(result.outputs.items()):
            print(f"  {label:<12} {path.name} ({path.stat().st_size} bytes)")

        if abs(result.ga_score - EXPECTED) > TOLERANCE:
            print(f"FAIL: GA-score differs from the reference value by more than {TOLERANCE}")
            return 1
        if "combined" not in result.outputs:
            print("FAIL: no combined overlay was written")
            return 1

        # The convergence log comes from the patch to glosa.cpp, so this doubles as a check that
        # the binary in the image is the patched one and not a stock rebuild.
        rows = glosa_runner.read_convergence(Path(tmp))
        print(f"  {'convergence':<12} {len(rows)} iterations logged")
        if not rows:
            print("FAIL: no convergence history; is the scorer built from the patched source?")
            return 1
        if abs(rows[-1]["best_score"] - result.ga_score) > TOLERANCE:
            print(f"FAIL: last logged best {rows[-1]['best_score']:.6f} is not the final score")
            return 1

    print("OK")
    return 0


if __name__ == "__main__":
    sys.exit(main())
