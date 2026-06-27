"""Step 02: run predictors and merge final epitope evidence.

This step should run predictor adapters on the de-duplicated task table, then
merge results back to epitope windows and mutation evidence.
"""

from __future__ import annotations

import argparse
from typing import Optional


def build_parser() -> argparse.ArgumentParser:
    """Create the command-line parser for prediction and final merge."""

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("-s", "--sample", required=True)
    parser.add_argument("--tasks", required=True)
    parser.add_argument("--epitope-windows", required=True)
    parser.add_argument("--variant-annotation", required=True)
    parser.add_argument("-o", "--outdir", required=True)
    parser.add_argument("--threads", type=int, default=1)
    return parser


def main(argv: Optional[list[str]] = None) -> int:
    """CLI entry point."""

    _ = build_parser().parse_args(argv)
    raise NotImplementedError("Step 02 skeleton only; implementation pending.")


if __name__ == "__main__":
    raise SystemExit(main())

