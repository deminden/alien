"""Core domain types shared by every stage, and the single implementation of ID normalization.

Human-specific assumptions live in exactly three places: :func:`ensembl_stable_id` (the ``ENSG`` prefix),
``alien.resources.hgnc`` and ``alien.resources.ncbi`` (taxon 9606).
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from enum import StrEnum
from pathlib import Path
from typing import Any, Final

import pandas as pd

ENSEMBL_GENE_PREFIX: Final = "ENSG"
_ENSEMBL_STABLE = re.compile(r"^(ENSG[0-9]+)(?:\.[0-9]+)?(?:_PAR_Y)?$")
_MISSING_TOKENS: Final = frozenset({"NA", "NAN", "NULL"})


class MappingStatus(StrEnum):
    """What happened to one gene membership in one target namespace."""

    DIRECT_SOURCE_ENSEMBL = "direct_source_ensembl"
    ARCHIVE_ENSEMBL_RESCUE = "archive_ensembl_rescue"
    HGNC_CURRENT_SYMBOL = "hgnc_current_symbol"
    HGNC_PREVIOUS_SYMBOL = "hgnc_previous_symbol"
    HGNC_ALIAS_SYMBOL = "hgnc_alias_symbol"
    HGNC_ENSEMBL_TIEBREAK = "hgnc_ensembl_tiebreak"
    MANUAL_SYMBOL_REPAIR = "manual_symbol_repair"
    NCBI_GENE_HISTORY = "ncbi_gene_history"
    NCBI_GENE_INFO_SYNONYM = "ncbi_gene_info_synonym"
    AMBIGUOUS_DROP = "ambiguous_drop"
    UNMAPPED_DROP = "unmapped_drop"
    NOT_IN_OUTPUT_GENES = "not_in_output_genes"
    NON_GENE_TOKEN_DROP = "non_gene_token_drop"

    @property
    def is_mapped(self) -> bool:
        """The membership contributes a gene to the output GMT."""
        return self in _MAPPED

    @property
    def is_rescue(self) -> bool:
        """Mapped through the Ensembl archive or NCBI rather than a direct or HGNC match."""
        return self in _RESCUES


_MAPPED: Final = frozenset(
    {
        MappingStatus.DIRECT_SOURCE_ENSEMBL,
        MappingStatus.ARCHIVE_ENSEMBL_RESCUE,
        MappingStatus.HGNC_CURRENT_SYMBOL,
        MappingStatus.HGNC_PREVIOUS_SYMBOL,
        MappingStatus.HGNC_ALIAS_SYMBOL,
        MappingStatus.HGNC_ENSEMBL_TIEBREAK,
        MappingStatus.MANUAL_SYMBOL_REPAIR,
        MappingStatus.NCBI_GENE_HISTORY,
        MappingStatus.NCBI_GENE_INFO_SYNONYM,
    }
)
_RESCUES: Final = frozenset(
    {MappingStatus.ARCHIVE_ENSEMBL_RESCUE, MappingStatus.NCBI_GENE_HISTORY, MappingStatus.NCBI_GENE_INFO_SYNONYM}
)


class RemovalReason(StrEnum):
    BROAD_DISEASE_TERM = "broad_disease_term"
    TOO_SMALL = "too_small"
    TOO_LARGE = "too_large"
    EXACT_DUPLICATE = "exact_duplicate"
    JACCARD_DUPLICATE = "jaccard_duplicate"
    NO_MAPPED_GENES = "no_mapped_genes"


class Family(StrEnum):
    """Term families with built-in behaviour; configs may use other family names too."""

    BIOLOGY_PROCESS_PATHWAY = "biology_process_pathway"
    BIOLOGY_FUNCTION_LOCATION = "biology_function_location"
    DISEASE_PHENOTYPE = "disease_phenotype"
    CANCER_DEPENDENCY_STATE = "cancer_dependency_state"


SYMBOLS_NAMESPACE: Final = "symbols"


def is_missing(value: Any) -> bool:
    if value is None:
        return True
    try:
        return bool(pd.isna(value))
    except (TypeError, ValueError):
        return False


def ensembl_stable_id(value: Any) -> str:
    """``ENSG00000141510.18`` -> ``ENSG00000141510``; non-Ensembl text is returned stripped."""
    if is_missing(value):
        return ""
    text = str(value).strip()
    if not text or text.upper() in _MISSING_TOKENS:
        return ""
    match = _ENSEMBL_STABLE.match(text)
    return match.group(1) if match else text


def is_ensembl_gene_id(value: Any) -> bool:
    return ensembl_stable_id(value).startswith(ENSEMBL_GENE_PREFIX)


def clean_field(value: Any) -> str:
    """Single-line text: tabs and newlines become spaces."""
    text = "" if is_missing(value) else str(value)
    return re.sub(r"[\t\r\n]+", " ", text).strip()


def sanitize_term_id(value: Any) -> str:
    """A GMT-safe identifier: letters, digits, ``_``, ``.`` and ``-`` only."""
    text = clean_field(value)
    text = re.sub(r"\s+", "_", text)
    text = re.sub(r"[^A-Za-z0-9_.-]+", "_", text)
    text = re.sub(r"_+", "_", text)
    return text.strip("_") or "UNNAMED"


@dataclass(slots=True)
class TermRecord:
    """One gene set in one namespace: its identity metadata and member genes."""

    term_id: str
    description: str
    meta: dict[str, str]
    genes: set[str] = field(default_factory=set)
    original_genes: set[str] = field(default_factory=set)

    @property
    def family(self) -> str:
        return self.meta.get("family", "")

    @property
    def source_tag(self) -> str:
        return self.meta.get("source_tag", "")


TermSet = dict[str, TermRecord]


@dataclass(frozen=True, slots=True)
class AnnotationCandidate:
    """A target gene that a symbol could map to."""

    ensembl_gene_id: str
    gene_biotype: str = ""
    is_in_output_genes: bool = True
    annotation_source: str = "annotation_gtf"


@dataclass(frozen=True)
class GeneIdList:
    """A resolved output_genes or gene_filter list for one target."""

    ids: frozenset[str] | None
    table: pd.DataFrame
    id_type: str
    source_path: Path | None = None
    id_column: str | None = None
    symbol_column: str | None = None


@dataclass(frozen=True)
class Target:
    """A prepared target namespace: its annotation helper and the IDs it may emit."""

    name: str
    annotation: pd.DataFrame
    annotation_label: str
    annotation_path: Path
    output_genes: GeneIdList
    gene_filter: GeneIdList
    primary_annotation_ids: frozenset[str]
    helper_ids: frozenset[str]  # GTF plus metadata fallbacks, before output-gene symbol rows are added
    metadata_fallbacks: tuple[dict[str, Any], ...] = ()

    @property
    def annotation_ids(self) -> frozenset[str]:
        """Stable IDs in the annotation used for mapping (GTF, fallbacks and output-gene symbol rows)."""
        return frozenset(self.annotation["ensembl_gene_id"].astype(str)) if not self.annotation.empty else frozenset()

    @property
    def output_gene_ids(self) -> frozenset[str]:
        return self.output_genes.ids or frozenset()

    @property
    def filter_ids(self) -> frozenset[str]:
        return self.gene_filter.ids or frozenset()

    @property
    def output_gene_source(self) -> str:
        if self.output_gene_ids:
            return "output_genes"
        if self.filter_ids:
            return "annotation_filtered"
        return "annotation"

    @property
    def effective_ids(self) -> frozenset[str]:
        """The IDs this target may write: output_genes, else annotation & gene_filter, else annotation."""
        if self.output_gene_ids:
            return self.output_gene_ids
        if self.filter_ids:
            return self.helper_ids & self.filter_ids
        return self.helper_ids
