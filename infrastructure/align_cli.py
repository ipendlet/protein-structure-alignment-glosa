"""Command-line front end to the same pipeline the web service runs.

The web service caps a job at `GLOSA_TIMEOUT` because a whole-protein pair can occupy the shared
endpoint for hours and spill a multi-gigabyte product graph.  Those runs belong here instead,
where the timeout defaults to off and the scratch files can be kept for inspection.

    python infrastructure/align_cli.py -s1 example/1IA1-bs.pdb -s2 example/2BL9-bs.pdb \
        -o runs/1IA1_vs_2BL9 --transfer example/2BL9-lig.pdb
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from glosa_runner import GlosaError, align


def main() -> int:
    parser = argparse.ArgumentParser(description="Align two structures with G-LoSA.")
    parser.add_argument("-s1", required=True, type=Path, help="reference structure (PDB)")
    parser.add_argument("-s2", required=True, type=Path, help="mobile structure (PDB)")
    parser.add_argument("-o", "--outdir", required=True, type=Path, help="output directory")
    parser.add_argument("--transfer", type=Path,
                        help="structure carried along by the matrix but not scored (-s2w)")
    parser.add_argument("--s1cf", type=Path,
                        help="override the reference chemical feature file; derived from -s1 "
                             "when omitted")
    parser.add_argument("--s2cf", type=Path,
                        help="override the mobile chemical feature file; derived from -s2 "
                             "when omitted")
    parser.add_argument("--timeout", type=int, default=0,
                        help="seconds before the run is killed; 0 means no limit")
    parser.add_argument("--keep-scratch", action="store_true",
                        help="keep pairs.rst and product_graph.rst, which can be very large")
    parser.add_argument("extra", nargs="*",
                        help="further glosa flags, e.g. -iter 1 -n 3")
    args = parser.parse_args()

    try:
        result = align(
            structure1=args.s1,
            structure2=args.s2,
            transfer=args.transfer,
            features1=args.s1cf,
            features2=args.s2cf,
            outdir=args.outdir,
            extra_args=args.extra,
            timeout=args.timeout or 2**31 - 1,
            keep_scratch=args.keep_scratch,
        )
    except GlosaError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1

    print(f"GA-score: {result.ga_score:.6f}")
    for label, path in sorted(result.outputs.items()):
        print(f"  {label:<12} {path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
