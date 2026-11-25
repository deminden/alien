"""Mapping: symbol resolution, the per-key decision chain, projection onto namespaces, and audit tables."""

from __future__ import annotations

import warnings
from typing import Any

from .project import MappingResult, Memberships, map_namespaces, problematic_source_ids, target_context
from .resolver import MappingContext, MappingOutcome, resolve
from .symbols import SymbolResolution, resolve_symbol

__all__ = [
    "MappingContext",
    "MappingOutcome",
    "MappingResult",
    "Memberships",
    "SymbolResolution",
    "map_namespaces",
    "problematic_source_ids",
    "resolve",
    "resolve_symbol",
    "target_context",
]

_REMOVED = {
    "map_ensembl": "alien.mapping.map_namespaces (or resolve() for a single key)",
    "map_symbols": "alien.mapping.map_namespaces",
}


def resolve_current_symbol(
    symbol: str, hgnc_maps: Any, manual_repairs: dict[str, str] | None = None, *, use_alias: bool = True
) -> tuple[str, str, list[str]]:
    """Deprecated 0.1.x API: ``(current symbol, status, candidates)``. Use :func:`resolve_symbol`."""
    warnings.warn(
        "resolve_current_symbol is deprecated; use alien.mapping.resolve_symbol.", DeprecationWarning, stacklevel=2
    )
    result = resolve_symbol(symbol, hgnc_maps, manual_repairs, use_alias=use_alias)
    return result.symbol, result.status, list(result.candidates)


def __getattr__(name: str) -> Any:
    if name == "EnsemblArchiveResolver":
        warnings.warn(
            "alien.mapping.EnsemblArchiveResolver moved to alien.resources.ensembl_archive.",
            DeprecationWarning,
            stacklevel=2,
        )
        from ..resources.ensembl_archive import EnsemblArchiveResolver  # noqa: PLC0415

        return EnsemblArchiveResolver
    if name in _REMOVED:
        raise AttributeError(f"alien.mapping.{name} was removed in 0.2; use {_REMOVED[name]}.")
    raise AttributeError(f"module 'alien.mapping' has no attribute {name!r}")
