"""Load, layer and validate ALIEN configuration.

Layering, lowest first: model defaults, the ``extends`` chain of the config file, the config file,
then each override in order. Mappings deep-merge; lists replace. After merging, the list-extension
keys ``gene_mapping.extra_non_gene_tokens``, ``gene_mapping.extra_non_gene_token_regex``,
``filtering.extra_broad_disease_regex`` and ``source_priority.extra.<family>`` append to their lists.

Relative paths in a YAML layer resolve against that file's directory (bundled presets excepted: their
paths stay relative to the working directory). Dict layers are taken as given.
"""

from __future__ import annotations

import os
from collections.abc import Callable, Sequence
from copy import deepcopy
from pathlib import Path
from typing import Any

import yaml
from pydantic import ValidationError

from ..errors import ConfigError
from .schema import SOURCE_TYPES, AlienConfig, default_layer

ConfigInput = str | Path | dict[str, Any]
LIST_EXTENSIONS = (
    ("gene_mapping", "extra_non_gene_tokens", "non_gene_tokens"),
    ("gene_mapping", "extra_non_gene_token_regex", "non_gene_token_regex"),
    ("filtering", "extra_broad_disease_regex", "broad_disease_regex"),
)

#: Set by alien.config.presets so ``extends: <preset>`` works without a circular import.
PRESET_RESOLVER: Callable[[str], Path | None] = lambda name: None  # noqa: E731


def load_config(
    path: ConfigInput | None = None,
    overrides: Sequence[ConfigInput] | None = None,
    *,
    resolve_relative_paths: bool = True,
) -> AlienConfig:
    """Return the validated config built from defaults, ``path`` and ``overrides``."""
    layers: list[dict[str, Any]] = []
    if path is not None:
        layers.append(read_layer(path, resolve_relative_paths=resolve_relative_paths))
    layers.extend(read_layer(override) for override in overrides or [])
    return validate_layers(layers)


def validate_layers(layers: Sequence[dict[str, Any]]) -> AlienConfig:
    merged = default_layer()
    for layer in layers:
        merged = merge_config(merged, layer)
    return validate_config(finalize_config(merged))


def validate_config(data: dict[str, Any]) -> AlienConfig:
    try:
        return AlienConfig.model_validate(data)
    except ValidationError as exc:
        raise ConfigError(format_validation_error(exc)) from None


def format_validation_error(exc: ValidationError) -> str:
    lines = []
    for error in exc.errors():
        location = _yaml_path(error["loc"])
        message = str(error["msg"]).removeprefix("Value error, ")
        if error["type"] == "extra_forbidden":
            message = "unknown option (check the spelling)"
        lines.append(f"{location or '<root>'}: {message}")
    return "Invalid ALIEN configuration:\n  " + "\n  ".join(lines)


def _yaml_path(loc: tuple[int | str, ...]) -> str:
    parts: list[str] = []
    for index, item in enumerate(loc):
        if isinstance(item, int):
            parts.append(f"[{item}]")
            continue
        previous = loc[index - 1] if index else None
        if isinstance(previous, int) and item in SOURCE_TYPES:
            continue  # discriminated-union tag, not a key the user wrote
        if item.startswith(("function-", "tagged-union")):
            continue
        parts.append(f".{item}" if parts else item)
    return "".join(parts)


def read_layer(
    source: ConfigInput,
    *,
    resolve_relative_paths: bool = True,
    _seen: frozenset[Path] = frozenset(),
) -> dict[str, Any]:
    """One config layer as a plain mapping, with its ``extends`` chain merged underneath."""
    if isinstance(source, dict):
        data = deepcopy(source)
        base_dir: Path | None = None
    else:
        path = Path(source).resolve()
        if path in _seen:
            raise ConfigError(f"Config extends chain loops back to {path}.")
        _seen = _seen | {path}
        data = read_config_file(path, resolve_relative_paths=resolve_relative_paths)
        base_dir = path.parent
    parent = data.pop("extends", None)
    if parent is None:
        return data
    if not isinstance(parent, str) or not parent.strip():
        raise ConfigError("extends must name a preset or a config file path.")
    preset = PRESET_RESOLVER(parent)
    if preset is not None:
        base = read_layer(preset, resolve_relative_paths=False, _seen=_seen)
    else:
        parent_path = Path(parent).expanduser()
        if not parent_path.is_absolute():
            parent_path = (base_dir or Path.cwd()) / parent_path
        if not parent_path.exists():
            raise ConfigError(f"extends target {parent!r} is neither a bundled preset nor an existing file.")
        base = read_layer(parent_path, _seen=_seen)
    return merge_config(base, data)


def read_config_file(path: str | Path, *, resolve_relative_paths: bool = True) -> dict[str, Any]:
    path = Path(path)
    try:
        with path.open("r", encoding="utf-8") as handle:
            data = yaml.safe_load(handle) or {}
    except FileNotFoundError:
        raise ConfigError(f"Config file {path} does not exist.") from None
    except yaml.YAMLError as exc:
        raise ConfigError(f"Config file {path} is not valid YAML: {exc}") from None
    if not isinstance(data, dict):
        raise ConfigError(f"Config file {path} must contain a YAML mapping.")
    if resolve_relative_paths:
        resolve_config_paths(data, path.resolve().parent)
    return data


def finalize_config(cfg: dict[str, Any]) -> dict[str, Any]:
    cfg = deepcopy(cfg)
    for section, source_key, target_key in LIST_EXTENSIONS:
        block = cfg.get(section)
        if isinstance(block, dict) and source_key in block:
            extras = block.pop(source_key) or []
            block[target_key] = [*(block.get(target_key) or []), *extras]
    priority = cfg.get("source_priority")
    if isinstance(priority, dict) and isinstance(priority.get("extra"), dict):
        for family, extras in priority.pop("extra").items():
            priority[family] = [*(priority.get(family) or []), *extras]
    return cfg


def merge_config(base: dict[str, Any], override: dict[str, Any]) -> dict[str, Any]:
    """Deep-merge mappings; any other value (lists included) replaces the base value."""
    merged = deepcopy(base)
    for key, value in override.items():
        if isinstance(value, dict) and isinstance(merged.get(key), dict):
            merged[key] = merge_config(merged[key], value)
        else:
            merged[key] = deepcopy(value)
    return merged


def _resolve_path_value(value: object, base_dir: Path) -> object:
    if not isinstance(value, str) or value.strip() in {"", "null"} or "://" in value:
        return value
    path = Path(value).expanduser()
    return os.path.normpath(path if path.is_absolute() else base_dir / path)


def _resolve_keys(mapping: object, keys: Sequence[str], base_dir: Path) -> None:
    if not isinstance(mapping, dict):
        return
    for key in keys:
        if key in mapping:
            value = mapping[key]
            if isinstance(value, list):
                mapping[key] = [_resolve_path_value(item, base_dir) for item in value]
            else:
                mapping[key] = _resolve_path_value(value, base_dir)


def resolve_config_paths(data: dict[str, Any], base_dir: Path) -> None:
    """Rewrite every path-typed field of one config layer relative to ``base_dir`` (in place)."""
    _resolve_keys(data.get("project"), ["source_dir", "outdir"], base_dir)
    for spec in data.get("sources") or []:
        _resolve_keys(spec, ["path", "paths", "cache_dir"], base_dir)
    for target in data.get("targets") or []:
        if isinstance(target, dict):
            _resolve_target_paths(target, base_dir)


def _resolve_target_paths(target: dict[str, Any], base_dir: Path) -> None:
    _resolve_keys(target, ["annotation_gtf"], base_dir)
    for key in ["output_genes", "gene_filter"]:
        if isinstance(target.get(key), str):
            target[key] = _resolve_path_value(target[key], base_dir)
        _resolve_keys(target.get(key), ["path", "file"], base_dir)
    annotation = target.get("annotation")
    _resolve_keys(annotation, ["path"], base_dir)
    if isinstance(annotation, dict):
        fallbacks = annotation.get("metadata_fallbacks") or []
        for fallback in fallbacks if isinstance(fallbacks, list) else [fallbacks]:
            _resolve_keys(fallback, ["path"], base_dir)
