from __future__ import annotations

from collections.abc import Sequence
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


REMOVED_OPTIONS: dict[tuple[str, str], str] = {
    ("redundancy", "apply_after_output_gene_projection"): "redundancy always runs after target projection",
    ("redundancy", "nested_overlap_filter"): "the nested-overlap filter was never implemented; remove the block",
    ("filtering", "remove_empty_after_output_gene_projection"): "empty terms are always dropped; use min_size_default",
    ("gene_mapping", "strip_ensembl_version"): "Ensembl versions are always stripped",
    ("gene_mapping", "allow_ambiguous_alias_mapping"): "ambiguous aliases are always dropped and audited",
    ("gene_mapping", "prefer_gencode_exact_gene_name"): "annotation gene names are always preferred",
    ("gene_mapping", "prefer_current_hgnc_symbol"): "current HGNC symbols are always preferred",
    ("gene_mapping", "drop_unmapped_from_ensg_gmt"): "unmapped genes are never written to Ensembl GMTs",
    ("ensembl_archive", "unresolved_example_threshold"): "use ensembl_archive.unresolved_print_limit",
}


def load_config(
    path: str | Path | None = None, overrides: Sequence[str | Path | dict[str, Any]] | None = None
) -> dict[str, Any]:
    cfg = default_config()
    if path is not None:
        cfg = merge_config(cfg, read_config_file(path))
    for override in overrides or []:
        cfg = merge_config(cfg, _coerce_override(override))
    return finalize_config(cfg)


def read_config_file(path: str | Path) -> dict[str, Any]:
    with Path(path).open("r", encoding="utf-8") as handle:
        data = yaml.safe_load(handle) or {}
    if not isinstance(data, dict):
        raise ValueError(f"Config file {path} must contain a YAML mapping.")
    return data


def finalize_config(cfg: dict[str, Any]) -> dict[str, Any]:
    cfg = deepcopy(cfg)
    _apply_list_extensions(cfg)
    _validate_config(cfg)
    return cfg


def merge_config(base: dict[str, Any], override: dict[str, Any]) -> dict[str, Any]:
    merged = deepcopy(base)
    _deep_update(merged, override)
    return merged


def _coerce_override(override: str | Path | dict[str, Any]) -> dict[str, Any]:
    if isinstance(override, dict):
        return override
    return read_config_file(override)


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
    for (section, key), replacement in REMOVED_OPTIONS.items():
        if isinstance(cfg.get(section), dict) and key in cfg[section]:
            raise ValueError(f"Config option {section}.{key} was removed: {replacement}.")
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
