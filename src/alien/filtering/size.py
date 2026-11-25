"""Size limits per family and removal of broad disease terms."""

from __future__ import annotations

import re
from dataclasses import dataclass

import pandas as pd

from ..config.schema import FilteringConfig
from ..model import Family, RemovalReason, TermRecord, TermSet
from ..schema import SIZE_REMOVED_COLUMNS


@dataclass
class SizeResult:
    kept: TermSet
    removed: pd.DataFrame


def max_size_for_family(config: FilteringConfig, family: str) -> int:
    if family == Family.DISEASE_PHENOTYPE:
        return config.max_size_disease
    if family == Family.CANCER_DEPENDENCY_STATE:
        return config.max_size_cancer
    return config.max_size_default


def removal_reason(record: TermRecord, config: FilteringConfig, broad: list[re.Pattern[str]]) -> str:
    if record.family == Family.DISEASE_PHENOTYPE and config.drop_broad_disease_terms:
        names = [record.meta.get("original_name", ""), record.meta.get("display_name", "")]
        if any(pattern.match(name.strip()) for pattern in broad for name in names):
            return RemovalReason.BROAD_DISEASE_TERM
    size = len(record.genes)
    if size < config.min_size_default:
        return RemovalReason.TOO_SMALL
    if size > max_size_for_family(config, record.family):
        return RemovalReason.TOO_LARGE
    return ""


def apply_size_filters(terms: TermSet, namespace: str, config: FilteringConfig) -> SizeResult:
    broad = [re.compile(pattern, re.IGNORECASE) for pattern in config.broad_disease_regex]
    kept: TermSet = {}
    removed: list[dict[str, object]] = []
    for term_id, record in sorted(terms.items()):
        reason = removal_reason(record, config, broad)
        if not reason:
            kept[term_id] = record
            continue
        removed.append(
            {
                "target_namespace": namespace,
                "family": record.family,
                "term_id": term_id,
                "original_name": record.meta.get("original_name", ""),
                "source_tag": record.source_tag,
                "original_size_symbols": len(record.original_genes),
                "mapped_size": len(record.genes),
                "reason": str(reason),
            }
        )
    return SizeResult(kept, pd.DataFrame(removed, columns=list(SIZE_REMOVED_COLUMNS)))
