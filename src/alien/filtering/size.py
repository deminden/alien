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


@dataclass(frozen=True)
class ExclusionRule:
    name: str
    source_tags: frozenset[str]
    patterns: tuple[re.Pattern[str], ...]

    def matches(self, record: TermRecord, names: list[str]) -> bool:
        if self.source_tags and record.source_tag not in self.source_tags:
            return False
        return any(pattern.search(name) for pattern in self.patterns for name in names)


def compile_exclusions(config: FilteringConfig) -> list[ExclusionRule]:
    return [
        ExclusionRule(
            name,
            frozenset(rule.source_tags),
            tuple(re.compile(pattern, re.IGNORECASE) for pattern in rule.patterns),
        )
        for name, rule in sorted(config.exclude_terms.items())
        if rule.enabled and rule.patterns
    ]


def term_names(record: TermRecord) -> list[str]:
    """The source and display names, each also with underscores read as spaces (MSigDB-style names)."""
    names: list[str] = []
    for key in ("original_name", "display_name"):
        name = str(record.meta.get(key, "") or "").strip()
        for variant in (name, re.sub(r"[_\s]+", " ", name).strip()):
            if variant and variant not in names:
                names.append(variant)
    return names


def removal_rule(
    record: TermRecord, config: FilteringConfig, broad: list[re.Pattern[str]], exclusions: list[ExclusionRule]
) -> tuple[str, str]:
    """(removal reason, exclusion rule name) for one term; ("", "") when it is kept."""
    names = term_names(record)
    for rule in exclusions:
        if rule.matches(record, names):
            return RemovalReason.EXCLUDED_TERM, rule.name
    return removal_reason(record, config, broad, names), ""


def removal_reason(
    record: TermRecord, config: FilteringConfig, broad: list[re.Pattern[str]], names: list[str] | None = None
) -> str:
    if record.family == Family.DISEASE_PHENOTYPE and config.drop_broad_disease_terms:
        names = term_names(record) if names is None else names
        if any(pattern.match(name) for pattern in broad for name in names):
            return RemovalReason.BROAD_DISEASE_TERM
    size = len(record.genes)
    if size < config.min_size_default:
        return RemovalReason.TOO_SMALL
    if size > max_size_for_family(config, record.family):
        return RemovalReason.TOO_LARGE
    return ""


def apply_size_filters(terms: TermSet, namespace: str, config: FilteringConfig) -> SizeResult:
    broad = [re.compile(pattern, re.IGNORECASE) for pattern in config.broad_disease_regex]
    exclusions = compile_exclusions(config)
    kept: TermSet = {}
    removed: list[dict[str, object]] = []
    for term_id, record in sorted(terms.items()):
        reason, rule = removal_rule(record, config, broad, exclusions)
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
                "exclusion_rule": rule,
            }
        )
    return SizeResult(kept, pd.DataFrame(removed, columns=list(SIZE_REMOVED_COLUMNS)))
