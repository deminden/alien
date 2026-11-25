"""Redundancy removal: exact duplicates, then Jaccard similarity above ``jaccard_cutoff``.

Candidate pairs come from a prefix-filter similarity join (AllPairs): genes are ranked by ascending
document frequency, and two sets can only have Jaccard similarity >= t if their first
``|S| - ceil(t*|S|) + 1`` genes overlap. Only those pairs are verified, so hub genes no longer cost
k^2 comparisons. The comparison is strict (``J > cutoff``), exactly as before.

Linkage:

* ``single`` (default): similar terms are joined transitively (union-find), so A~B and B~C put A and C
  in one cluster even when A and C are dissimilar. One representative per cluster is kept.
* ``representative``: terms are visited in representative order; a term is removed only when it is
  similar to a term already kept, so every removed term is above the cutoff against its representative.

Representative order: source priority (``source_priority`` for the family; unlisted sources last), then
larger mapped size, then names not starting with a generic prefix, then shorter name, then term ID.
"""

from __future__ import annotations

import logging
import math
from collections import Counter
from collections.abc import Callable, Iterable, Sequence
from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass
from typing import Final

import pandas as pd

from ..config.schema import AlienConfig, RedundancyConfig
from ..model import RemovalReason, TermSet
from ..parallel import process_pool_context
from ..schema import REDUNDANCY_REMOVED_COLUMNS, REDUNDANCY_SUMMARY_COLUMNS, SIZE_REMOVED_COLUMNS
from .size import apply_size_filters

LOGGER = logging.getLogger("alien")
ALL_FAMILIES: Final = "*"
GENERIC_PREFIXES: Final = ("regulation of", "positive regulation of", "negative regulation of", "process")
#: Tolerance that keeps the prefix length conservative when t*|S| lands on an integer in floating point.
PREFIX_EPSILON: Final = 1e-9


@dataclass
class RedundancyResult:
    kept: TermSet
    removed: pd.DataFrame
    summary: pd.DataFrame


def representative_key(terms: TermSet, priority: Sequence[str]) -> Callable[[str], tuple[object, ...]]:
    rank = {source: index for index, source in enumerate(priority)}
    unlisted = len(priority)

    def key(term_id: str) -> tuple[object, ...]:
        record = terms[term_id]
        name = str(record.meta.get("display_name") or record.meta.get("original_name") or term_id)
        generic = int(name.lower().startswith(GENERIC_PREFIXES))
        return (rank.get(record.source_tag, unlisted), -len(record.genes), generic, len(name), term_id)

    return key


def _removed(namespace: str, family: str, term_id: str, representative: str, kind: str, size: int) -> dict[str, object]:
    return {
        "target_namespace": namespace,
        "family": family,
        "term_id": term_id,
        "representative_term_id": representative,
        "redundancy_cluster_id": f"{namespace}:{family}:{kind}:{representative}",
        "reason": str(RemovalReason.EXACT_DUPLICATE if kind == "exact" else RemovalReason.JACCARD_DUPLICATE),
        "mapped_size": size,
    }


def remove_exact_duplicates(
    terms: TermSet, namespace: str, family: str, priority: Sequence[str]
) -> tuple[TermSet, list[dict[str, object]]]:
    groups: dict[frozenset[str], list[str]] = {}
    for term_id, record in terms.items():
        groups.setdefault(frozenset(record.genes), []).append(term_id)
    key = representative_key(terms, priority)
    kept: TermSet = {}
    removed: list[dict[str, object]] = []
    for genes, term_ids in groups.items():
        representative = min(term_ids, key=key)
        kept[representative] = terms[representative]
        removed.extend(
            _removed(namespace, family, term_id, representative, "exact", len(genes))
            for term_id in sorted(set(term_ids) - {representative})
        )
    return kept, removed


def _prefix_length(size: int, cutoff: float) -> int:
    return max(1, min(size, size - math.ceil(cutoff * size - PREFIX_EPSILON) + 1))


def _is_similar(a: frozenset[str], b: frozenset[str], cutoff: float) -> bool:
    """Jaccard(a, b) > cutoff, evaluated exactly as in 0.1.x (size-ratio bound first, then float ratio)."""
    size_a, size_b = len(a), len(b)
    if min(size_a, size_b) / max(size_a, size_b) <= cutoff:
        return False
    inter = len(a & b)
    union = size_a + size_b - inter
    return bool(union) and inter / union > cutoff


def _prefixes(gene_sets: Sequence[frozenset[str]], cutoff: float) -> list[list[int]]:
    """Each set's prefix of gene ranks, genes ordered by ascending document frequency."""
    frequency = Counter(gene for genes in gene_sets for gene in genes)
    rank = {gene: index for index, gene in enumerate(sorted(frequency, key=lambda g: (frequency[g], g)))}
    return [sorted(rank[gene] for gene in genes)[: _prefix_length(len(genes), cutoff)] for genes in gene_sets]


def similar_pairs(gene_sets: Sequence[frozenset[str]], cutoff: float) -> list[tuple[int, int]]:
    """Index pairs (i < j) with Jaccard similarity strictly above ``cutoff``."""
    prefixes = _prefixes(gene_sets, cutoff)
    index: dict[int, list[int]] = {}
    pairs: list[tuple[int, int]] = []
    for x in sorted(range(len(gene_sets)), key=lambda i: (len(gene_sets[i]), i)):
        candidates = {y for token in prefixes[x] for y in index.get(token, ())} if gene_sets[x] else set()
        pairs += [(min(x, y), max(x, y)) for y in candidates if _is_similar(gene_sets[x], gene_sets[y], cutoff)]
        for token in prefixes[x] if gene_sets[x] else ():
            index.setdefault(token, []).append(x)
    return sorted(pairs)


def _find(parent: list[int], item: int) -> int:
    root = item
    while parent[root] != root:
        root = parent[root]
    while parent[item] != root:
        parent[item], item = root, parent[item]
    return root


def single_linkage_clusters(n: int, pairs: Iterable[tuple[int, int]]) -> list[list[int]]:
    """Connected components, ordered by their first member."""
    parent = list(range(n))
    for a, b in pairs:
        root_a, root_b = _find(parent, a), _find(parent, b)
        if root_a != root_b:
            parent[max(root_a, root_b)] = min(root_a, root_b)
    components: dict[int, list[int]] = {}
    for item in range(n):
        components.setdefault(_find(parent, item), []).append(item)
    return list(components.values())


def _representative_linkage(
    terms: TermSet, term_ids: list[str], pairs: list[tuple[int, int]], key: Callable[[str], tuple[object, ...]]
) -> tuple[TermSet, list[tuple[str, str]]]:
    """Visit terms in representative order; remove a term only if it is similar to an already kept term."""
    neighbours: dict[int, list[int]] = {}
    for a, b in pairs:
        neighbours.setdefault(a, []).append(b)
        neighbours.setdefault(b, []).append(a)
    position = {term_id: i for i, term_id in enumerate(term_ids)}
    kept: TermSet = {}
    removed: list[tuple[str, str]] = []
    for term_id in sorted(term_ids, key=key):
        kept_neighbours = [term_ids[n] for n in neighbours.get(position[term_id], []) if term_ids[n] in kept]
        if kept_neighbours:
            removed.append((term_id, min(kept_neighbours, key=key)))
        else:
            kept[term_id] = terms[term_id]
    return kept, removed


def _single_linkage(
    terms: TermSet, term_ids: list[str], pairs: list[tuple[int, int]], key: Callable[[str], tuple[object, ...]]
) -> tuple[TermSet, list[tuple[str, str]]]:
    kept: TermSet = {}
    removed: list[tuple[str, str]] = []
    for component in single_linkage_clusters(len(term_ids), pairs):
        members = sorted(term_ids[i] for i in component)
        representative = min(members, key=key)
        kept[representative] = terms[representative]
        removed += [(term_id, representative) for term_id in members if term_id != representative]
    return kept, removed


def remove_jaccard(
    terms: TermSet, namespace: str, family: str, priority: Sequence[str], config: RedundancyConfig
) -> tuple[TermSet, list[dict[str, object]]]:
    term_ids = list(terms)
    pairs = similar_pairs([frozenset(terms[t].genes) for t in term_ids], config.jaccard_cutoff)
    linkage = _representative_linkage if config.linkage == "representative" else _single_linkage
    kept, removed = linkage(terms, term_ids, pairs, representative_key(terms, priority))
    rows = [_removed(namespace, family, t, rep, "jaccard", len(terms[t].genes)) for t, rep in removed]
    return kept, rows


def group_priority(terms: TermSet, priority: dict[str, list[str]], family: str) -> list[str]:
    if family != ALL_FAMILIES:
        return list(priority.get(family, []))
    merged: list[str] = []
    for name in sorted({record.family for record in terms.values()}):
        merged.extend(tag for tag in priority.get(name, []) if tag not in merged)
    return merged


def remove_redundancy(
    terms: TermSet, namespace: str, family: str, config: RedundancyConfig, priority: dict[str, list[str]]
) -> RedundancyResult:
    """Redundancy removal for one (namespace, family) group; ``family`` may be ``*`` for all families."""
    order = group_priority(terms, priority, family)
    if config.exact_duplicate_removal:
        after_exact, exact_removed = remove_exact_duplicates(terms, namespace, family, order)
    else:
        after_exact, exact_removed = dict(terms), []
    kept, jaccard_removed = remove_jaccard(after_exact, namespace, family, order, config)
    summary = {
        "target_namespace": namespace,
        "family": family,
        "n_input": len(terms),
        "n_after_exact": len(after_exact),
        "n_after_jaccard": len(kept),
        "n_removed_exact": len(exact_removed),
        "n_removed_jaccard": len(jaccard_removed),
        "jaccard_cutoff": config.jaccard_cutoff,
    }
    return RedundancyResult(
        kept,
        pd.DataFrame(exact_removed + jaccard_removed, columns=list(REDUNDANCY_REMOVED_COLUMNS)),
        pd.DataFrame([summary], columns=list(REDUNDANCY_SUMMARY_COLUMNS)),
    )


def redundancy_groups(terms: TermSet, per_family: bool) -> list[tuple[str, TermSet]]:
    if not terms:
        return []
    if not per_family:
        return [(ALL_FAMILIES, dict(terms))]
    families = sorted({record.family for record in terms.values()})
    return [(family, {t: r for t, r in terms.items() if r.family == family}) for family in families]


@dataclass(frozen=True)
class _Task:
    namespace: str
    family: str
    terms: TermSet
    config: RedundancyConfig
    priority: dict[str, list[str]]


def _run_task(task: _Task) -> RedundancyResult:
    return remove_redundancy(task.terms, task.namespace, task.family, task.config, task.priority)


def run_tasks[T, R](function: Callable[[T], R], tasks: Sequence[T], workers: int) -> list[R]:
    """Run ``function`` over ``tasks`` in order, in a process pool when ``workers > 1``."""
    if workers <= 1 or len(tasks) <= 1:
        return [function(task) for task in tasks]
    with ProcessPoolExecutor(max_workers=min(workers, len(tasks)), mp_context=process_pool_context()) as pool:
        return list(pool.map(function, tasks))


@dataclass
class FilterResult:
    size_filtered: dict[str, TermSet]
    filtered: dict[str, TermSet]
    removed_size: pd.DataFrame
    removed_redundancy: pd.DataFrame
    redundancy_summary: pd.DataFrame


def _concat(frames: Sequence[pd.DataFrame], columns: Sequence[str]) -> pd.DataFrame:
    kept = [frame for frame in frames if not frame.empty]
    return pd.concat(kept, ignore_index=True)[list(columns)] if kept else pd.DataFrame(columns=list(columns))


def filter_namespaces(mapped: dict[str, TermSet], config: AlienConfig, workers: int) -> FilterResult:
    """Size filters per namespace, then redundancy per (namespace, family) group, serially or in a pool."""
    size_filtered: dict[str, TermSet] = {}
    removed_size: list[pd.DataFrame] = []
    tasks: list[_Task] = []
    for namespace, terms in mapped.items():
        LOGGER.info("Filtering namespace %s with %d projected terms", namespace, len(terms))
        size = apply_size_filters(terms, namespace, config.filtering)
        size_filtered[namespace] = size.kept
        removed_size.append(size.removed)
        for family, subset in redundancy_groups(size.kept, config.redundancy.apply_per_family):
            tasks.append(_Task(namespace, family, subset, config.redundancy, dict(config.source_priority)))
    if tasks:
        LOGGER.info("Removing redundancy across %d namespace/family groups with %d worker(s)", len(tasks), workers)
    results = run_tasks(_run_task, tasks, workers)
    filtered: dict[str, TermSet] = {namespace: {} for namespace in mapped}
    for task, result in zip(tasks, results, strict=True):
        filtered[task.namespace].update(result.kept)
    return FilterResult(
        size_filtered,
        filtered,
        _concat(removed_size, SIZE_REMOVED_COLUMNS),
        _concat([r.removed for r in results], REDUNDANCY_REMOVED_COLUMNS),
        _concat([r.summary for r in results], REDUNDANCY_SUMMARY_COLUMNS),
    )
