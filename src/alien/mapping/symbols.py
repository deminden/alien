"""Symbol resolution against HGNC, shared by the symbols namespace and every target namespace."""

from __future__ import annotations

import re
from collections.abc import Callable, Mapping, Sequence
from dataclasses import dataclass
from typing import Final

from ..config.schema import GeneMappingConfig
from ..resources.hgnc import HgncMaps

#: ambiguous resolution status -> mapping_stage written to ambiguous_gene_mappings
AMBIGUOUS_STAGES: Final = {
    "ambiguous_alias": "hgnc_alias_symbol",
    "ambiguous_previous": "hgnc_previous_symbol",
    "ambiguous_case_insensitive": "hgnc_case_insensitive_symbol",
    "ambiguous_excel_date": "excel_date_symbol",
}
#: resolution status -> MappingStatus value when the symbol maps to a target gene
SYMBOL_STATUSES: Final = {
    "exact_current_symbol": "hgnc_current_symbol",
    "hgnc_previous_symbol": "hgnc_previous_symbol",
    "hgnc_alias_symbol": "hgnc_alias_symbol",
    "manual_symbol_repair": "manual_symbol_repair",
    "hgnc_case_insensitive_symbol": "hgnc_case_insensitive_symbol",
    "mitochondrial_legacy_symbol": "mitochondrial_legacy_symbol",
    "excel_date_symbol": "excel_date_symbol",
}
#: Resolution status for an alias-only match shorter than ``min_alias_length`` (treated as unmapped).
SHORT_ALIAS: Final = "short_alias_rejected"
#: Legacy mitochondrial-genome gene names (tRNA one-letter codes, ND/COX/ATP/CYTB/RNR) -> current HGNC symbols.
#: COX1 and COX2 are left out on purpose: in disease text they usually mean PTGS1/PTGS2.
MITOCHONDRIAL_LEGACY: Final = {
    **{f"TRN{code}": f"MT-T{code}" for code in "ARNDCEQGHIKMFPTWYV"},
    "TRNL1": "MT-TL1",
    "TRNL2": "MT-TL2",
    "TRNS1": "MT-TS1",
    "TRNS2": "MT-TS2",
    **{f"ND{n}": f"MT-ND{n}" for n in range(1, 7)},
    "ND4L": "MT-ND4L",
    "COX3": "MT-CO3",
    "CYTB": "MT-CYB",
    "ATP6": "MT-ATP6",
    "ATP8": "MT-ATP8",
    "RNR1": "MT-RNR1",
    "RNR2": "MT-RNR2",
}
_EXCEL_DATE = re.compile(r"^(?:(\d{1,2})-(SEPT?|MAR)|(SEPT?|MAR)-(\d{1,2}))$", re.IGNORECASE)
#: Month -> gene-family prefix; numbers listed in ``_EXCEL_DATE_AMBIGUOUS`` also name another family.
_EXCEL_DATE_PREFIX: Final = {"SEP": ("SEPTIN",), "SEPT": ("SEPTIN",), "MAR": ("MARCHF",)}
_EXCEL_DATE_AMBIGUOUS: Final = {("MAR", 1): "MTARC1", ("MAR", 2): "MTARC2"}


@dataclass(frozen=True, slots=True)
class SymbolOptions:
    """How far symbol resolution may go beyond exact HGNC matches."""

    use_alias: bool = True
    case_insensitive: bool = False
    min_alias_length: int = 1
    mitochondrial_legacy: bool = False
    excel_dates: bool = False

    @classmethod
    def from_config(cls, config: GeneMappingConfig, *, short_alias_guard: bool = True) -> SymbolOptions:
        """Options for one source; ``short_alias_guard`` says whether it is listed in short_alias_source_tags."""
        return cls(
            use_alias=config.use_hgnc_alias_mapping,
            case_insensitive=config.case_insensitive_symbols,
            min_alias_length=config.min_alias_length if short_alias_guard else 1,
            mitochondrial_legacy=config.mitochondrial_legacy_symbols,
            excel_dates=config.excel_date_symbols,
        )


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


def _resolve_folded(symbol: str, hgnc: HgncMaps, options: SymbolOptions) -> SymbolResolution | None:
    """The same order ignoring letter case: current, previous, then alias symbols; several genes are ambiguous."""
    folded = symbol.casefold()
    for kind in ("current", "previous", "alias") if options.use_alias else ("current", "previous"):
        matches = hgnc.folded[kind].get(folded)
        if not matches:
            continue
        if len(matches) == 1:
            if kind == "alias" and len(symbol) < options.min_alias_length:
                return SymbolResolution("", SHORT_ALIAS, matches)
            return SymbolResolution(matches[0], "hgnc_case_insensitive_symbol", matches)
        return SymbolResolution("", "ambiguous_case_insensitive", matches)
    return None


def _resolve_excel_date(symbol: str, hgnc: HgncMaps) -> SymbolResolution | None:
    """``1-Sep`` -> SEPTIN1, ``MAR-6`` -> MARCHF6: gene names a spreadsheet converted to dates."""
    match = _EXCEL_DATE.match(symbol)
    if match is None:
        return None
    number_text = match.group(1) or match.group(4)
    month = (match.group(2) or match.group(3)).upper()
    number = int(number_text)
    candidates = [f"{prefix}{number}" for prefix in _EXCEL_DATE_PREFIX[month]]
    other = _EXCEL_DATE_AMBIGUOUS.get((month[:3], number))
    if other:
        candidates.append(other)
    known = tuple(sorted(c for c in candidates if c in hgnc.current))
    if not known:
        return None
    if len(known) > 1:
        return SymbolResolution("", "ambiguous_excel_date", known)
    return SymbolResolution(known[0], "excel_date_symbol", known)


def resolve_symbol(
    raw: str,
    hgnc: HgncMaps,
    repairs: Mapping[str, str] | None = None,
    *,
    use_alias: bool = True,
    case_insensitive: bool = False,
    options: SymbolOptions | None = None,
) -> SymbolResolution:
    """Exact current symbol, manual repair, mitochondrial legacy name, previous symbol, ambiguous previous,
    ambiguous alias, alias, then (when enabled) a case-insensitive lookup and an Excel-date repair.

    A manual repair outranks every history or alias lookup but never an exact current symbol. With
    ``case_insensitive``, a symbol that none of these match is looked up again ignoring letter case
    (``C1ORF112`` finds ``C1orf112``); an exact match of any kind always wins over a case-insensitive one.
    Alias-only matches shorter than ``min_alias_length`` are rejected (status ``short_alias_rejected``).
    """
    if options is None:
        options = SymbolOptions(use_alias=use_alias, case_insensitive=case_insensitive)
    symbol = str(raw).strip()
    if not symbol:
        return SymbolResolution("", "unmapped", ())
    steps: tuple[Callable[[], SymbolResolution | None], ...] = (
        lambda: _resolve_exact(symbol, hgnc, repairs, options),
        lambda: _resolve_history(symbol, hgnc, options),
        lambda: _resolve_folded(symbol, hgnc, options) if options.case_insensitive else None,
        lambda: _resolve_excel_date(symbol, hgnc) if options.excel_dates else None,
    )
    for step in steps:
        if (resolution := step()) is not None:
            return resolution
    return SymbolResolution("", "unmapped", ())


def _resolve_exact(
    symbol: str, hgnc: HgncMaps, repairs: Mapping[str, str] | None, options: SymbolOptions
) -> SymbolResolution | None:
    """Exact current symbol, then a manual repair, then a legacy mitochondrial name."""
    if symbol in hgnc.current:
        return SymbolResolution(symbol, "exact_current_symbol", (symbol,))
    if repairs and symbol in repairs:
        repaired = str(repairs[symbol]).strip()
        known = repaired in hgnc.current
        return SymbolResolution(repaired if known else "", "manual_symbol_repair" if known else "unmapped", (repaired,))
    legacy = MITOCHONDRIAL_LEGACY.get(symbol, "") if options.mitochondrial_legacy else ""
    if legacy in hgnc.current:
        return SymbolResolution(legacy, "mitochondrial_legacy_symbol", (legacy,))
    return None


def _resolve_history(symbol: str, hgnc: HgncMaps, options: SymbolOptions) -> SymbolResolution | None:
    """Previous, ambiguous previous, ambiguous alias, then alias symbols (short aliases rejected)."""
    for status, table in _history_tables(hgnc, options.use_alias):
        if symbol not in table:
            continue
        match = table[symbol]
        if isinstance(match, list):
            return SymbolResolution("", status, tuple(match))
        if status == "hgnc_alias_symbol" and len(symbol) < options.min_alias_length:
            return SymbolResolution("", SHORT_ALIAS, (match,))
        return SymbolResolution(match, status, (match,))
    return None


def short_alias_guarded(config: GeneMappingConfig, source_tags: Sequence[str]) -> list[bool]:
    """Per source tag: whether the short-alias rule applies to it."""
    listed = set(config.short_alias_source_tags)
    return ["*" in listed or tag in listed for tag in source_tags]


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
