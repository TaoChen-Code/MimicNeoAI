# coding=utf-8
from __future__ import annotations
import argparse
import sys
from typing import List, Callable


def _call_main(func: Callable, rest: List[str]) -> int:
    """
    Call a module-level main with argv passthrough if supported.
    Falls back to sys.argv patching if the signature is main().
    """
    try:
        # Prefer main(argv: Optional[List[str]] = None) style
        return int(func(rest) or 0)
    except TypeError:
        # Fallback: emulate CLI argv
        old_argv = sys.argv[:]
        try:
            sys.argv = [sys.argv[0]] + rest
            return int(func() or 0)
        finally:
            sys.argv = old_argv


def _run_download_database(rest: List[str]) -> int:
    from mimicneoai import download_database
    return _call_main(download_database.main, rest)


def _run_cryptic(rest: List[str]) -> int:
    from mimicneoai.cryptic_pipeline.cryptic import main as cryptic_main
    return _call_main(cryptic_main, rest)


def _run_microbial(rest: List[str]) -> int:
    from mimicneoai.microbial_pipeline.microbial import main as microbial_main
    return _call_main(microbial_main, rest)


def _run_mutation_derived(rest: List[str]) -> int:
    # Now that the module is mutation_derived.py, a clean import is enough
    from mimicneoai.mutation_derived_pipeline.mutation_derived import main as md_main
    return _call_main(md_main, rest)


def main(argv: List[str] | None = None) -> int:
    argv = argv if argv is not None else sys.argv[1:]

    parser = argparse.ArgumentParser(prog="mimicneoai", description="MimicNeoAI unified CLI")
    subparsers = parser.add_subparsers(dest="command", required=True)

    sp_dl = subparsers.add_parser("download_database",
                                  help="Download and extract the reference database")
    sp_dl.add_argument("rest", nargs=argparse.REMAINDER,
                       help="Args passed to the download module (e.g., --target-dir PATH)")

    sp_cryptic = subparsers.add_parser("cryptic",
                                       help="Run the cryptic pipeline")
    sp_cryptic.add_argument("rest", nargs=argparse.REMAINDER,
                            help="Args passed to cryptic (e.g., -c ... )")

    sp_micro = subparsers.add_parser("microbial",
                                     help="Run the microbial pipeline")
    sp_micro.add_argument("rest", nargs=argparse.REMAINDER,
                          help="Args passed to microbial (e.g., -c ... )")

    sp_md = subparsers.add_parser("mutation-derived",
                                  help="Run the mutation-derived pipeline")
    sp_md.add_argument("rest", nargs=argparse.REMAINDER,
                       help="Args passed to mutation-derived (e.g., -c ... )")

    args = parser.parse_args(argv)

    if args.command == "download_database":
        return _run_download_database(args.rest)
    if args.command == "cryptic":
        return _run_cryptic(args.rest)
    if args.command == "microbial":
        return _run_microbial(args.rest)
    if args.command == "mutation-derived":
        return _run_mutation_derived(args.rest)

    parser.print_help()
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
