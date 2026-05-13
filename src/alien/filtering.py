from __future__ import annotations

import re
from typing import Any

import pandas as pd


def apply_size_and_broad_filters(
    term_genes: dict[str, dict[str, Any]],
    namespace: str,
    cfg: dict,
) -> tuple[dict[str, dict[str, Any]], pd.DataFrame]:
    # Remove terms outside configured size limits
    kept: dict[str, dict[str, Any]] = {}
    removed: list[dict[str, object]] = []
    min_size = int(cfg["filtering"]["min_size_default"])
    broad_patterns = [re.compile(pattern, re.I) for pattern in cfg["filtering"].get("broad_disease_regex", [])]

    for term_id, record in sorted(term_genes.items()):
        meta = record["meta"]
        family = meta.get("family", "")
        max_size = _max_size_for_family(cfg, family)
        size = len(record["genes"])
        reason = ""
        if family == "disease_phenotype" and cfg["filtering"].get("drop_broad_disease_terms", True):
            names = [str(meta.get("original_name", "")), str(meta.get("display_name", ""))]
            if any(pattern.match(name.strip()) for pattern in broad_patterns for name in names):
                reason = "broad_disease_term"
        if not reason and size < min_size:
            reason = "too_small"
        if not reason and size > max_size:
            reason = "too_large"
        if reason:
            removed.append(_removed_row(namespace, term_id, record, reason))
        else:
            kept[term_id] = record
    return kept, pd.DataFrame(removed)


def _max_size_for_family(cfg: dict, family: str) -> int:
    if family == "disease_phenotype":
        return int(cfg["filtering"].get("max_size_disease", cfg["filtering"]["max_size_default"]))
    if family == "cancer_dependency_state":
        return int(cfg["filtering"].get("max_size_cancer", cfg["filtering"]["max_size_default"]))
    return int(cfg["filtering"]["max_size_default"])


def _removed_row(namespace: str, term_id: str, record: dict[str, Any], reason: str) -> dict[str, object]:
    meta = record["meta"]
    return {
        "target_namespace": namespace,
        "family": meta.get("family", ""),
        "term_id": term_id,
        "original_name": meta.get("original_name", ""),
        "source_tag": meta.get("source_tag", ""),
        "original_size_symbols": len(record.get("original_genes", [])),
        "mapped_size_before_universe": len(record.get("genes", [])),
        "mapped_size_after_universe": len(record.get("genes", [])),
        "reason": reason,
    }


GENERIC_PREFIXES = ("regulation of", "positive regulation of", "negative regulation of", "process")


def remove_redundancy(
    term_genes: dict[str, dict[str, Any]],
    namespace: str,
    cfg: dict,
) -> tuple[dict[str, dict[str, Any]], pd.DataFrame, pd.DataFrame]:
    # Remove duplicate terms separately inside each family
    jaccard = float(cfg["redundancy"]["jaccard_cutoff"])
    priority = cfg.get("source_priority", {})
    kept: dict[str, dict[str, Any]] = {}
    removed_rows: list[dict[str, object]] = []
    summary_rows: list[dict[str, object]] = []

    for family in sorted({record["meta"].get("family", "") for record in term_genes.values()}):
        subset = {tid: rec for tid, rec in term_genes.items() if rec["meta"].get("family", "") == family}
        after_exact, exact_removed = _remove_exact_duplicates(subset, namespace, family, priority.get(family, []))
        after_jaccard, jaccard_removed = _remove_jaccard(after_exact, namespace, family, priority.get(family, []), jaccard)
        kept.update(after_jaccard)
        removed_rows.extend(exact_removed + jaccard_removed)
        summary_rows.append(
            {
                "target_namespace": namespace,
                "family": family,
                "n_input": len(subset),
                "n_after_exact": len(after_exact),
                "n_after_jaccard": len(after_jaccard),
                "n_removed_exact": len(exact_removed),
                "n_removed_jaccard": len(jaccard_removed),
                "jaccard_cutoff": jaccard,
            }
        )
    return kept, pd.DataFrame(removed_rows), pd.DataFrame(summary_rows)


def _remove_exact_duplicates(
    terms: dict[str, dict[str, Any]],
    namespace: str,
    family: str,
    priority: list[str],
) -> tuple[dict[str, dict[str, Any]], list[dict[str, object]]]:
    # Collapse terms with identical gene membership
    by_genes: dict[tuple[str, ...], list[str]] = {}
    for term_id, record in terms.items():
        by_genes.setdefault(tuple(sorted(record["genes"])), []).append(term_id)
    kept: dict[str, dict[str, Any]] = {}
    removed: list[dict[str, object]] = []
    for gene_key, term_ids in by_genes.items():
        representative = _choose_representative(term_ids, terms, priority)
        kept[representative] = terms[representative]
        cluster = f"{namespace}:{family}:exact:{representative}"
        for term_id in sorted(set(term_ids) - {representative}):
            removed.append(_redundancy_removed(namespace, family, term_id, representative, cluster, "exact_duplicate", len(gene_key)))
    return kept, removed


def _remove_jaccard(
    terms: dict[str, dict[str, Any]],
    namespace: str,
    family: str,
    priority: list[str],
    cutoff: float,
) -> tuple[dict[str, dict[str, Any]], list[dict[str, object]]]:
    # Cluster highly overlapping terms by Jaccard similarity
    parent = {term_id: term_id for term_id in terms}
    rank = {term_id: 0 for term_id in terms}
    gene_sets = {term_id: record["genes"] for term_id, record in terms.items()}
    sizes = {term_id: len(record["genes"]) for term_id, record in terms.items()}
    inverted: dict[str, list[str]] = {}
    for term_id, record in terms.items():
        for gene in record["genes"]:
            inverted.setdefault(gene, []).append(term_id)

    compared: set[tuple[str, str]] = set()
    for term_ids in inverted.values():
        term_ids = sorted(set(term_ids))
        for i, a in enumerate(term_ids):
            for b in term_ids[i + 1 :]:
                pair = (a, b)
                if pair in compared:
                    continue
                compared.add(pair)
                size_a = sizes[a]
                size_b = sizes[b]
                if min(size_a, size_b) / max(size_a, size_b) <= cutoff:
                    continue
                if size_a < size_b:
                    inter = sum(1 for gene in gene_sets[a] if gene in gene_sets[b])
                else:
                    inter = sum(1 for gene in gene_sets[b] if gene in gene_sets[a])
                union = size_a + size_b - inter
                if union and inter / union > cutoff:
                    _union(parent, rank, a, b)

    components: dict[str, list[str]] = {}
    for term_id in terms:
        components.setdefault(_find(parent, term_id), []).append(term_id)

    kept: dict[str, dict[str, Any]] = {}
    removed: list[dict[str, object]] = []
    for component in components.values():
        members = sorted(component)
        representative = _choose_representative(members, terms, priority)
        kept[representative] = terms[representative]
        cluster = f"{namespace}:{family}:jaccard:{representative}"
        for term_id in sorted(set(members) - {representative}):
            removed.append(_redundancy_removed(namespace, family, term_id, representative, cluster, "jaccard_duplicate", len(terms[term_id]["genes"])))
    return kept, removed


def _find(parent: dict[str, str], item: str) -> str:
    root = item
    while parent[root] != root:
        root = parent[root]
    while parent[item] != item:
        item, parent[item] = parent[item], root
    return root


def _union(parent: dict[str, str], rank: dict[str, int], a: str, b: str) -> None:
    root_a = _find(parent, a)
    root_b = _find(parent, b)
    if root_a == root_b:
        return
    if rank[root_a] < rank[root_b]:
        root_a, root_b = root_b, root_a
    parent[root_b] = root_a
    if rank[root_a] == rank[root_b]:
        rank[root_a] += 1


def _choose_representative(term_ids: list[str], terms: dict[str, dict[str, Any]], priority: list[str]) -> str:
    priority_index = {source: i for i, source in enumerate(priority)}

    def key(term_id: str) -> tuple[object, ...]:
        record = terms[term_id]
        meta = record["meta"]
        name = str(meta.get("display_name") or meta.get("original_name") or term_id)
        lower = name.lower()
        generic = int(any(lower.startswith(prefix) for prefix in GENERIC_PREFIXES))
        return (
            priority_index.get(meta.get("source_tag", ""), 999),
            -len(record["genes"]),
            generic,
            len(name),
            term_id,
        )

    return sorted(term_ids, key=key)[0]


def _redundancy_removed(namespace: str, family: str, term_id: str, representative: str, cluster: str, reason: str, size: int) -> dict[str, object]:
    return {
        "target_namespace": namespace,
        "family": family,
        "term_id": term_id,
        "representative_term_id": representative,
        "redundancy_cluster_id": cluster,
        "reason": reason,
        "mapped_size_after_universe": size,
    }
