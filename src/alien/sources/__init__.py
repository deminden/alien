"""Gene-set sources: a registry of loaders that each return the canonical membership table.

Built-in loaders register themselves with :func:`register`. A project can add its own source without
forking ALIEN: publish a callable ``load(options: dict, cache: ResourceCache) -> DataFrame`` under the
``alien.sources`` entry-point group and configure it as ``{type: plugin, plugin: <entry-point name>,
options: {...}}``.
"""

from __future__ import annotations

import logging
from importlib.metadata import entry_points
from pathlib import Path
from typing import Any

import pandas as pd

from ..cache import ResourceCache
from ..config.schema import AlienConfig, PluginSource
from ..errors import ConfigError, SourceError, TermIdCollisionError
from ..io.tables import write_table
from ..schema import CANONICAL_COLUMNS, TERM_ID_COLLISION_COLUMNS, TERM_ID_IDENTITY_COLUMNS
from . import enrichr, local, msigdb  # noqa: F401  (importing them registers the built-in loaders)
from .frames import concat_memberships, empty_memberships
from .registry import loader_for, register, registered_types

LOGGER = logging.getLogger("alien")
ENTRY_POINT_GROUP = "alien.sources"


@register("plugin")
def load_plugin(spec: PluginSource, cache: ResourceCache) -> pd.DataFrame:
    matches = [ep for ep in entry_points(group=ENTRY_POINT_GROUP) if ep.name == spec.plugin]
    if not matches:
        raise ConfigError(
            f"No source plugin named {spec.plugin!r} is installed (entry-point group {ENTRY_POINT_GROUP})."
        )
    frame = matches[0].load()(dict(spec.options), cache)
    if not isinstance(frame, pd.DataFrame):
        raise SourceError(f"Source plugin {spec.plugin!r} did not return a DataFrame.")
    missing = [column for column in ("term_id", "gene_symbol") if column not in frame]
    if missing:
        raise SourceError(f"Source plugin {spec.plugin!r} returned no {', '.join(missing)} column.")
    return concat_memberships([frame])


def load_source_memberships(config: AlienConfig, cache: ResourceCache) -> pd.DataFrame:
    """Load every enabled source into one canonical membership table."""
    memberships = concat_memberships([loader_for(spec.type)(spec, cache) for spec in config.sources if spec.enabled])
    if memberships.empty:
        raise SourceError("No source gene-set memberships were loaded.")
    return memberships


SAMPLE_SYMBOLS = 12


def _distinct(values: pd.Series) -> list[str]:
    return sorted({str(value) for value in values if str(value)})


def _n_distinct(values: pd.Series) -> int:
    return len(_distinct(values))


def _sample(values: pd.Series) -> str:
    return ";".join(_distinct(values)[:SAMPLE_SYMBOLS])


def _identity_counts(working: pd.DataFrame, columns: list[str]) -> pd.DataFrame:
    identity = working[["term_id", *columns]].drop_duplicates()
    return identity.groupby("term_id", sort=False).size().rename("n_identities_for_term").reset_index()


def term_id_collision_audit(memberships: pd.DataFrame) -> pd.DataFrame:
    """Rows describing each distinct identity behind a term_id that more than one identity uses."""
    columns = [column for column in TERM_ID_IDENTITY_COLUMNS if column in memberships.columns]
    audit_columns = [
        c for c in TERM_ID_COLLISION_COLUMNS if c != "action" and (c in columns or c not in TERM_ID_IDENTITY_COLUMNS)
    ]
    working = (
        memberships.reindex(columns=list(dict.fromkeys(["term_id", "gene_symbol", *columns]))).fillna("").astype(str)
    )
    counts = _identity_counts(working, columns)
    collided = set(counts.loc[counts["n_identities_for_term"] > 1, "term_id"])
    if not collided:
        return pd.DataFrame(columns=audit_columns)
    grouped = (
        working[working["term_id"].isin(collided)]
        .groupby(["term_id", *columns], dropna=False, sort=False)
        .agg(
            n_rows=("gene_symbol", "size"),
            n_gene_symbols=("gene_symbol", _n_distinct),
            sample_gene_symbols=("gene_symbol", _sample),
        )
        .reset_index()
        .merge(counts, on="term_id", how="left")
        .sort_values(["term_id", *columns], kind="mergesort")
        .reset_index(drop=True)
    )
    grouped["collision_identity_id"] = grouped.groupby("term_id").cumcount() + 1
    return grouped[audit_columns]


def check_term_id_collisions(
    memberships: pd.DataFrame, action: str, metadata_dir: Path, report_metadata_dir: Path, warnings: list[str]
) -> None:
    """Write term_id_collisions.tsv when term IDs collide; raise for action 'error', else warn."""
    audit = term_id_collision_audit(memberships)
    if audit.empty:
        return
    audit = audit.assign(action=action)
    name = "term_id_collisions.tsv"
    write_table(audit, metadata_dir / name)
    message = (
        f"Detected {len(audit)} conflicting term/source identities across {audit['term_id'].nunique()} "
        "term_id value(s)."
    )
    if action == "error":
        raise TermIdCollisionError(
            f"{message} See {report_metadata_dir / name}. Rename or prefix the colliding term IDs, "
            "or set term_id_collisions.action explicitly."
        )
    message = f"{message} See metadata/{name}."
    warnings.append(f"{message} Continuing because term_id_collisions.action is {action!r}.")
    LOGGER.warning("%s Continuing because term_id_collisions.action is %r.", message, action)


_MOVED = {
    "build_hgnc_maps": ("alien.resources.hgnc", "build_hgnc_maps"),
    "read_output_genes": ("alien.resources.output_genes", "read_gene_id_list"),
    "gencode_annotation_url": ("alien.resources.annotation", "gencode_annotation_url"),
    "resolve_enrichr_library": ("alien.sources.enrichr", "resolve_enrichr_library"),
    "read_symbol_gmt_source": ("alien.sources.local", "read_symbol_gmt"),
    "read_canonical_memberships": ("alien.sources.local", "read_canonical_memberships"),
    "read_msigdb_like_table": ("alien.sources.msigdb", "read_msigdb_like_table"),
}


def __getattr__(name: str) -> Any:
    """Deprecated 0.1.x import paths; removed in 0.3."""
    if name in _MOVED:
        import importlib  # noqa: PLC0415
        import warnings  # noqa: PLC0415

        module, attribute = _MOVED[name]
        warnings.warn(f"alien.sources.{name} moved to {module}.{attribute}.", DeprecationWarning, stacklevel=2)
        return getattr(importlib.import_module(module), attribute)
    raise AttributeError(f"module 'alien.sources' has no attribute {name!r}")


__all__ = [
    "CANONICAL_COLUMNS",
    "check_term_id_collisions",
    "empty_memberships",
    "load_source_memberships",
    "register",
    "registered_types",
    "term_id_collision_audit",
]
