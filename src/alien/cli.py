from __future__ import annotations

import argparse
from pathlib import Path

from .build import build
from .presets import is_preset_name, list_presets


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="alien", description="ALIEN: Audited Library Integration for External Namespaces")
    subparsers = parser.add_subparsers(dest="command", required=True)

    subparsers.add_parser("presets", help="List bundled ALIEN build presets")

    build_parser = subparsers.add_parser("build", help="Build audited combined GMT files")
    build_parser.add_argument("config_or_preset", nargs="?", help="Bundled preset name/alias or YAML configuration file")
    build_parser.add_argument("--config", type=Path, default=None, help="YAML configuration file")
    build_parser.add_argument("--override", type=Path, action="append", default=[], help="YAML overlay merged after the preset/config; repeatable")
    build_parser.add_argument("--source-dir", type=Path, default=None, help="Override project.source_dir from the config")
    build_parser.add_argument("--outdir", type=Path, default=None, help="Override project.outdir from the config")
    build_parser.add_argument("--workers", type=_positive_int, default=None, help="Number of local worker processes")
    build_parser.add_argument("--output-mode", choices=["full", "minimal", "gmt"], default="full", help="Which outputs to keep after the build")
    build_parser.add_argument("--dry-run", action="store_true", help="Validate configuration and report planned output")
    build_parser.add_argument("--force-download", action="store_true", help="Refresh mapping-resource downloads and caches")

    args = parser.parse_args(argv)
    if args.command == "presets":
        _print_presets()
        return 0

    if args.command == "build":
        if args.config is not None and args.config_or_preset is not None:
            build_parser.error("use either positional CONFIG_OR_PRESET or --config PATH, not both")
        if args.config is None and args.config_or_preset is None:
            build_parser.error("build requires a preset name, config path, or --config PATH")
        config = args.config if args.config is not None else _resolve_config_or_preset(args.config_or_preset)
        result = build(
            config=config,
            overrides=args.override,
            source_dir=args.source_dir,
            outdir=args.outdir,
            workers=args.workers,
            dry_run=args.dry_run,
            force_download=args.force_download,
            output_mode=args.output_mode,
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


def _resolve_config_or_preset(value: str | None) -> str | Path:
    if value is None:
        raise ValueError("Missing config or preset.")
    if is_preset_name(value):
        return value
    return Path(value)


def _print_presets() -> None:
    print("name\taliases\tdescription\tpath")
    for preset in list_presets():
        aliases = ", ".join(preset.aliases) if preset.aliases else "-"
        print(f"{preset.name}\t{aliases}\t{preset.description}\t{preset.path}")


if __name__ == "__main__":
    raise SystemExit(main())
