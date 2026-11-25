"""Symbol resolution against HGNC, shared by the symbols namespace and every target namespace."""

from __future__ import annotations

import re
from collections.abc import Callable, Mapping
from dataclasses import dataclass
from typing import Final

from ..config.schema import GeneMappingConfig
from ..resources.hgnc import HgncMaps

#: ambiguous resolution status -> mapping_stage written to ambiguous_gene_mappings
AMBIGUOUS_STAGES: Final = {"ambiguous_alias": "hgnc_alias_symbol", "ambiguous_previous": "hgnc_previous_symbol"}
#: resolution status -> MappingStatus value when the symbol maps to a target gene
SYMBOL_STATUSES: Final = {
    "exact_current_symbol": "hgnc_current_symbol",
    "hgnc_previous_symbol": "hgnc_previous_symbol",
    "hgnc_alias_symbol": "hgnc_alias_symbol",
    "manual_symbol_repair": "manual_symbol_repair",
}


@dataclass(frozen=True, slots=True)
class SymbolResolution:
    """``symbol`` is the current HGNC symbol, or "" when unmapped or ambiguous."""

    symbol: str
    status: str
    candidates: tuple[str, ...]

    @property
    def is_ambiguous(self) -> bool:
        return self.status in AMBIGUOUS_STAGES


def _history_tables(hgnc: HgncMaps, use_alias: bool) -> list[tuple[str, Mapping[str, str | list[str]]]]:
    tables: list[tuple[str, Mapping[str, str | list[str]]]] = [
        ("hgnc_previous_symbol", hgnc.previous),
        ("ambiguous_previous", hgnc.ambiguous_previous),
    ]
    if use_alias:
        tables += [("ambiguous_alias", hgnc.ambiguous_alias), ("hgnc_alias_symbol", hgnc.alias)]
    return tables


def resolve_symbol(
    raw: str, hgnc: HgncMaps, repairs: Mapping[str, str] | None = None, *, use_alias: bool = True
) -> SymbolResolution:
    """Exact current symbol, then a manual repair, previous symbol, ambiguous previous, ambiguous alias, alias.

    A manual repair outranks every history or alias lookup but never an exact current symbol.
    """
    symbol = str(raw).strip()
    if symbol and symbol in hgnc.current:
        return SymbolResolution(symbol, "exact_current_symbol", (symbol,))
    if symbol and repairs and symbol in repairs:
        repaired = str(repairs[symbol]).strip()
        known = repaired in hgnc.current
        return SymbolResolution(repaired if known else "", "manual_symbol_repair" if known else "unmapped", (repaired,))
    for status, table in _history_tables(hgnc, use_alias) if symbol else []:
        if symbol in table:
            match = table[symbol]
            if isinstance(match, list):
                return SymbolResolution("", status, tuple(match))
            return SymbolResolution(match, status, (match,))
    return SymbolResolution("", "unmapped", ())


def non_gene_matcher(config: GeneMappingConfig) -> Callable[[object], bool]:
    """Case-insensitive exact tokens or regex prefixes that are not genes (e.g. "disease")."""
    exact = {token.casefold() for token in config.non_gene_tokens}
    patterns = [re.compile(pattern, re.IGNORECASE) for pattern in config.non_gene_token_regex]

    def is_non_gene(value: object) -> bool:
        text = str(value).strip()
        if not text:
            return False
        return text.casefold() in exact or any(pattern.match(text) for pattern in patterns)

    return is_non_gene
