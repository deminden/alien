"""Command-line interface. Argument parsing and output only; every decision lives in the library."""

from __future__ import annotations

import argparse
import json
import logging
import sys
from collections.abc import Callable, Sequence
from pathlib import Path
from typing import Any

from .config import AlienConfig, ConfigInput, is_preset_name, list_presets
from .errors import AlienError
from .outputs import verify_output_dir
from .pipeline import build, effective_config_text, fetch, resolve_build_config
from .version import __version__

LOGGER = logging.getLogger("alien")
LOG_FORMAT = "%(asctime)s %(levelname)s %(message)s"


def _positive_int(value: str) -> int:
    try:
        parsed = int(value)
    except ValueError as exc:
        raise argparse.ArgumentTypeError("must be a positive integer") from exc
    if parsed < 1:
        raise argparse.ArgumentTypeError("must be a positive integer")
    return parsed


def _add_config_arguments(parser: argparse.ArgumentParser, *, with_config_flag: bool = True) -> None:
    parser.add_argument("config_or_preset", nargs="?", help="Bundled preset name/alias or YAML configuration file")
    if with_config_flag:
        parser.add_argument("--config", type=Path, default=None, help="YAML configuration file")
    parser.add_argument(
        "--override", type=Path, action="append", default=[], help="YAML overlay merged after the config; repeatable"
    )
    parser.add_argument("--source-dir", type=Path, default=None, help="Override project.source_dir")


def _add_cache_arguments(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "--lock", type=Path, default=None, help="Lock file (default: alien.lock.json beside the config)"
    )
    parser.add_argument("--force-download", action="store_true", help="Refresh downloads and derived caches")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="alien", description="ALIEN: Audited Library Integration for External Namespaces"
    )
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    parser.add_argument("-v", "--verbose", action="store_true", help="Debug logging and full tracebacks")
    parser.add_argument("-q", "--quiet", action="store_true", help="Only log warnings and errors")
    parser.add_argument("--log-file", type=Path, default=None, help="Also write the log to this file")
    sub = parser.add_subparsers(dest="command", required=True)

    sub.add_parser("presets", help="List bundled ALIEN build presets")

    build_cmd = sub.add_parser("build", help="Build audited combined GMT files")
    _add_config_arguments(build_cmd)
    _add_cache_arguments(build_cmd)
    build_cmd.add_argument("--outdir", type=Path, default=None, help="Override project.outdir")
    build_cmd.add_argument("--workers", type=_positive_int, default=None, help="Number of local worker processes")
    build_cmd.add_argument(
        "--output-mode", choices=["full", "minimal", "gmt"], default=None, help="Override outputs.mode"
    )
    build_cmd.add_argument("--dry-run", action="store_true", help="Validate inputs and print the build plan")
    build_cmd.add_argument("--force", action="store_true", help="Replace --outdir even if it is not an ALIEN build")
    build_cmd.add_argument("--offline", action="store_true", help="Never use the network; everything must be cached")
    build_cmd.add_argument(
        "--frozen", action="store_true", help="Offline, and every cached resource must match the lock"
    )

    fetch_cmd = sub.add_parser("fetch", help="Download every resource a config needs and write the lock file")
    _add_config_arguments(fetch_cmd)
    _add_cache_arguments(fetch_cmd)

    validate_cmd = sub.add_parser("validate", help="Validate a configuration without building")
    _add_config_arguments(validate_cmd)

    verify_cmd = sub.add_parser("verify", help="Re-hash an output directory against its MANIFEST.json")
    verify_cmd.add_argument("outdir", type=Path)

    config_cmd = sub.add_parser("config", help="Inspect configuration")
    config_sub = config_cmd.add_subparsers(dest="config_command", required=True)
    show = config_sub.add_parser("show", help="Print the effective configuration as YAML")
    _add_config_arguments(show, with_config_flag=False)
    config_sub.add_parser("schema", help="Print the configuration JSON Schema")
    return parser


def _config_source(args: argparse.Namespace, parser: argparse.ArgumentParser) -> ConfigInput | None:
    flag = getattr(args, "config", None)
    if flag is not None and args.config_or_preset is not None:
        parser.error("use either positional CONFIG_OR_PRESET or --config PATH, not both")
    if flag is not None:
        return Path(flag)
    value = args.config_or_preset
    if value is None:
        return None
    return value if is_preset_name(value) else Path(value)


def _require_config(args: argparse.Namespace, parser: argparse.ArgumentParser) -> ConfigInput:
    config = _config_source(args, parser)
    if config is None:
        parser.error(f"{args.command} requires a preset name, config path, or --config PATH")
    return config


def _common(args: argparse.Namespace) -> dict[str, Any]:
    return {"overrides": args.override, "source_dir": args.source_dir}


def _cmd_build(args: argparse.Namespace, parser: argparse.ArgumentParser) -> int:
    result = build(
        _require_config(args, parser),
        **_common(args),
        outdir=args.outdir,
        workers=args.workers,
        output_mode=args.output_mode,
        dry_run=args.dry_run,
        force=args.force,
        force_download=args.force_download,
        offline=args.offline,
        frozen=args.frozen,
        lock=args.lock,
    )
    if args.dry_run:
        print("\n".join(result.plan))
        return 0
    print(f"Wrote ALIEN outputs under {result.outdir}")
    if result.namespaces:
        print("Namespaces: " + ", ".join(result.namespaces))
    return 0


def _cmd_fetch(args: argparse.Namespace, parser: argparse.ArgumentParser) -> int:
    lock = fetch(_require_config(args, parser), **_common(args), lock=args.lock, force_download=args.force_download)
    print(f"Wrote {lock}")
    return 0


def _cmd_validate(args: argparse.Namespace, parser: argparse.ArgumentParser) -> int:
    model = resolve_build_config(_require_config(args, parser), **_common(args))
    print(f"Configuration is valid: {len(model.sources)} source(s); namespaces: {', '.join(model.namespaces())}")
    return 0


def _cmd_verify(args: argparse.Namespace, parser: argparse.ArgumentParser) -> int:
    problems = verify_output_dir(args.outdir)
    for problem in problems:
        print(problem)
    if problems:
        return 3
    print(f"{args.outdir}: every file matches MANIFEST.json")
    return 0


def _cmd_config(args: argparse.Namespace, parser: argparse.ArgumentParser) -> int:
    if args.config_command == "schema":
        print(json.dumps(AlienConfig.model_json_schema(), indent=2, sort_keys=True))
        return 0
    model = resolve_build_config(_config_source(args, parser), **_common(args))
    print(effective_config_text(model, keep_outdir=True), end="")
    return 0


def _cmd_presets(args: argparse.Namespace, parser: argparse.ArgumentParser) -> int:
    print("name\taliases\tdescription\tpath")
    for preset in list_presets():
        aliases = ", ".join(preset.aliases) if preset.aliases else "-"
        print(f"{preset.name}\t{aliases}\t{preset.description}\t{preset.path}")
    return 0


COMMANDS: dict[str, Callable[[argparse.Namespace, argparse.ArgumentParser], int]] = {
    "build": _cmd_build,
    "fetch": _cmd_fetch,
    "validate": _cmd_validate,
    "verify": _cmd_verify,
    "config": _cmd_config,
    "presets": _cmd_presets,
}


def configure_logging(verbose: bool, quiet: bool, log_file: Path | None) -> None:
    level = logging.DEBUG if verbose else logging.WARNING if quiet else logging.INFO
    handlers: list[logging.Handler] = [logging.StreamHandler(sys.stderr)]
    if log_file is not None:
        handlers.append(logging.FileHandler(log_file, encoding="utf-8"))
    for handler in handlers:
        handler.setFormatter(logging.Formatter(LOG_FORMAT, datefmt="%H:%M:%S"))
    LOGGER.handlers = [h for h in LOGGER.handlers if isinstance(h, logging.NullHandler)] + handlers
    LOGGER.setLevel(level)
    LOGGER.propagate = False


def main(argv: Sequence[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    configure_logging(args.verbose, args.quiet, args.log_file)
    try:
        return COMMANDS[args.command](args, parser)
    except AlienError as exc:
        if args.verbose:
            raise
        print(f"alien: error: {exc}", file=sys.stderr)
        return exc.exit_code
    except KeyboardInterrupt:
        print("alien: interrupted", file=sys.stderr)
        return 130
    except Exception as exc:
        if args.verbose:
            raise
        print(
            f"alien: unexpected error: {type(exc).__name__}: {exc} (rerun with --verbose for a traceback)",
            file=sys.stderr,
        )
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
