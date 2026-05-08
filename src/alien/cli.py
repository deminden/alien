from __future__ import annotations

import argparse
from pathlib import Path

from .build import build


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="alien", description="ALIEN: Audited Library Integration for External Namespaces")
    subparsers = parser.add_subparsers(dest="command", required=True)

    build_parser = subparsers.add_parser("build", help="Build audited combined GMT files")
    build_parser.add_argument("--config", type=Path, required=True, help="YAML configuration file")
    build_parser.add_argument("--outdir", type=Path, default=None, help="Override project.outdir from the config")
    build_parser.add_argument("--workers", type=_positive_int, default=None, help="Number of local worker processes")
    build_parser.add_argument("--dry-run", action="store_true", help="Validate configuration and report planned output")
    build_parser.add_argument("--force-download", action="store_true", help="Refresh mapping-resource downloads and caches")

    args = parser.parse_args(argv)
    if args.command == "build":
        result = build(
            config=args.config,
            outdir=args.outdir,
            workers=args.workers,
            dry_run=args.dry_run,
            force_download=args.force_download,
        )
        print(f"Wrote ALIEN outputs under {result.outdir}")
        if result.namespaces:
            print("Namespaces: " + ", ".join(result.namespaces))
        return 0
    parser.error(f"Unknown command: {args.command}")
    return 2


def _positive_int(value: str) -> int:
    try:
        parsed = int(value)
    except ValueError as exc:
        raise argparse.ArgumentTypeError("must be a positive integer") from exc
    if parsed < 1:
        raise argparse.ArgumentTypeError("must be a positive integer")
    return parsed


if __name__ == "__main__":
    raise SystemExit(main())
