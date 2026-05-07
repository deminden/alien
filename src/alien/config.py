from __future__ import annotations

from copy import deepcopy
from pathlib import Path
from typing import Any

import yaml

from .defaults import default_config


REQUIRED_SECTIONS = [
    "sources",
    "targets",
    "filtering",
    "redundancy",
    "source_priority",
    "gene_mapping",
    "ncbi_gene",
]


def load_config(path: str | Path | None = None) -> dict[str, Any]:
    cfg = default_config()
    if path is not None:
        with Path(path).open("r", encoding="utf-8") as handle:
            override = yaml.safe_load(handle) or {}
        cfg = merge_config(cfg, override)
    _apply_list_extensions(cfg)
    _validate_config(cfg)
    return cfg


def merge_config(base: dict[str, Any], override: dict[str, Any]) -> dict[str, Any]:
    merged = deepcopy(base)
    _deep_update(merged, override)
    return merged


def _deep_update(target: dict[str, Any], override: dict[str, Any]) -> None:
    for key, value in override.items():
        if isinstance(value, dict) and isinstance(target.get(key), dict):
            _deep_update(target[key], value)
        else:
            target[key] = deepcopy(value)


def _apply_list_extensions(cfg: dict[str, Any]) -> None:
    gene_mapping = cfg.get("gene_mapping", {})
    for source_key, target_key in [
        ("extra_non_gene_tokens", "non_gene_tokens"),
        ("extra_non_gene_token_regex", "non_gene_token_regex"),
    ]:
        extras = gene_mapping.pop(source_key, None)
        if extras:
            gene_mapping.setdefault(target_key, [])
            gene_mapping[target_key].extend(extras)

    filtering = cfg.get("filtering", {})
    extras = filtering.pop("extra_broad_disease_regex", None)
    if extras:
        filtering.setdefault("broad_disease_regex", [])
        filtering["broad_disease_regex"].extend(extras)

    source_priority = cfg.get("source_priority", {})
    extras_by_family = source_priority.pop("extra", None)
    if isinstance(extras_by_family, dict):
        for family, extras in extras_by_family.items():
            source_priority.setdefault(family, [])
            source_priority[family].extend(extras)


def _validate_config(cfg: dict[str, Any]) -> None:
    missing = [section for section in REQUIRED_SECTIONS if section not in cfg]
    if missing:
        raise ValueError(f"Config is missing required sections: {', '.join(missing)}")
    if not cfg.get("sources"):
        raise ValueError("Config must define at least one source.")
    include_symbols = bool(cfg.get("outputs", {}).get("include_symbols", True))
    if not include_symbols and not cfg.get("targets"):
        raise ValueError("Config must enable symbols output or define at least one target namespace.")
    for target in cfg.get("targets", []):
        if not target.get("name"):
            raise ValueError("Every target must define a non-empty name.")
        target_type = str(target.get("type", "ensembl_gtf"))
        if target_type != "ensembl_gtf":
            raise ValueError(f"Target {target.get('name', '<unnamed>')} has unsupported type {target_type!r}.")
        annotation = target.get("annotation", {})
        has_structured_annotation = isinstance(annotation, dict) and (
            annotation.get("path") or annotation.get("url") or annotation.get("version")
        )
        has_legacy_annotation = target.get("annotation_gtf") or target.get("gencode_version")
        if not has_structured_annotation and not has_legacy_annotation:
            raise ValueError(
                f"Target {target.get('name', '<unnamed>')} must define annotation.path, annotation.url, annotation.version, "
                "annotation_gtf, or gencode_version."
            )
