"""QC summaries: collection, mapping, targets, archive, NCBI and non-gene tokens."""

from __future__ import annotations

import logging
from collections.abc import Sequence

import pandas as pd

from ..model import SYMBOLS_NAMESPACE, MappingStatus, Target, TermSet, ensembl_stable_id
from ..schema import (
    ARCHIVE_SUMMARY_COLUMNS,
    COLLECTION_SUMMARY_COLUMNS,
    MAPPING_SUMMARY_COLUMNS,
    METADATA_FALLBACK_COLUMNS,
    NCBI_SUMMARY_COLUMNS,
    NON_GENE_SUMMARY_COLUMNS,
    TARGET_GENE_ID_COLUMNS,
    TARGET_NAMESPACE_SUMMARY_COLUMNS,
)

LOGGER = logging.getLogger("alien")
GroupKey = tuple[str, str]


def group_terms(terms: TermSet) -> dict[GroupKey, TermSet]:
    grouped: dict[GroupKey, TermSet] = {}
    for term_id, record in terms.items():
        grouped.setdefault((record.family, record.source_tag), {})[term_id] = record
    return grouped


def collection_summary(
    pre: dict[str, TermSet], post_size: dict[str, TermSet], final: dict[str, TermSet]
) -> pd.DataFrame:
    rows = []
    for namespace, terms in final.items():
        final_groups, size_groups, pre_groups = (
            group_terms(terms),
            group_terms(post_size.get(namespace, {})),
            group_terms(pre.get(namespace, {})),
        )
        for key in sorted(set(pre_groups) | set(size_groups) | set(final_groups)):
            subset = final_groups.get(key, {})
            sizes = [len(record.genes) for record in subset.values()]
            rows.append(
                {
                    "target_namespace": namespace,
                    "family": key[0],
                    "source_tag": key[1],
                    "n_terms_before_filter": len(pre_groups.get(key, {})),
                    "n_terms_after_size_filter": len(size_groups.get(key, {})),
                    "n_terms_after_redundancy": len(subset),
                    "median_size": float(pd.Series(sizes).median()) if sizes else 0,
                    "min_size": min(sizes) if sizes else 0,
                    "max_size": max(sizes) if sizes else 0,
                    "n_unique_genes": len({gene for record in subset.values() for gene in record.genes}),
                }
            )
    return pd.DataFrame(rows, columns=list(COLLECTION_SUMMARY_COLUMNS))


def _key(key: object) -> tuple[str, ...]:
    return tuple(str(part) for part in key) if isinstance(key, tuple) else (str(key),)


def _counts(frame: pd.DataFrame, keys: Sequence[str]) -> dict[tuple[str, ...], int]:
    if frame.empty:
        return {}
    return {_key(key): int(n) for key, n in frame.groupby(list(keys)).size().items()}


def mapping_summary(
    memberships: pd.DataFrame,
    mapped: dict[str, TermSet],
    tables: dict[str, pd.DataFrame],
    removed_size: pd.DataFrame,
    status: pd.DataFrame,
) -> pd.DataFrame:
    """Raw mapping events and unique genes per namespace, family and source tag."""
    identity = ["target_namespace", "family", "source_tag"]
    original = _counts(memberships, ["family", "source_tag"])
    unmapped, ambiguous, non_gene = (_counts(tables[name], identity) for name in ("unmapped", "ambiguous", "non_gene"))
    removed = _counts(removed_size, [*identity, "reason"])
    mapped_status = (
        status[status["mapping_status"].map(lambda value: MappingStatus(value).is_mapped)]
        if not status.empty
        else status
    )
    mapped_raw = {
        _key(key): int(n)
        for key, n in (
            mapped_status.groupby(identity)["n_gene_members"].sum().items() if not mapped_status.empty else []
        )
    }
    rows = []
    for namespace, terms in mapped.items():
        grouped = group_terms(terms)
        for family, source_tag in sorted(original):
            key = (namespace, family, source_tag)
            n_original = original[(family, source_tag)]
            n_non_gene = non_gene.get(key, 0)
            n_considered = max(0, n_original - n_non_gene)
            n_raw = n_considered if namespace == SYMBOLS_NAMESPACE else mapped_raw.get(key, 0)
            n_unique = sum(len(record.genes) for record in grouped.get((family, source_tag), {}).values())
            rows.append(
                {
                    "target_namespace": namespace,
                    "family": family,
                    "source_tag": source_tag,
                    "n_gene_members_original": n_original,
                    "n_non_gene_token_dropped": n_non_gene,
                    "n_gene_members_considered_for_mapping": n_considered,
                    "n_gene_members_mapped_raw": n_raw,
                    "n_gene_members_mapped_unique": n_unique,
                    "n_gene_members_collapsed_after_mapping": max(0, n_raw - n_unique),
                    "n_gene_members_unmapped": unmapped.get(key, 0),
                    "n_ambiguous": ambiguous.get(key, 0),
                    "mapping_rate_raw": round(n_raw / n_considered, 6) if n_considered else 0,
                    "mapping_rate_unique": round(n_unique / n_considered, 6) if n_considered else 0,
                    "n_terms_lost_because_too_small": removed.get((*key, "too_small"), 0),
                    "n_terms_lost_because_too_large": removed.get((*key, "too_large"), 0),
                }
            )
    return pd.DataFrame(rows, columns=list(MAPPING_SUMMARY_COLUMNS))


def _grouped_count(frame: pd.DataFrame, columns: Sequence[str]) -> pd.DataFrame:
    if frame.empty:
        return pd.DataFrame(columns=list(columns))
    return frame.groupby(list(columns[:-1])).size().reset_index(name=columns[-1])


def archive_summary(audit: pd.DataFrame) -> pd.DataFrame:
    return _grouped_count(audit, ARCHIVE_SUMMARY_COLUMNS)


def non_gene_summary(non_gene: pd.DataFrame) -> pd.DataFrame:
    return _grouped_count(non_gene, NON_GENE_SUMMARY_COLUMNS)


def ncbi_summary(audit: pd.DataFrame) -> pd.DataFrame:
    return _grouped_count(audit, NCBI_SUMMARY_COLUMNS)


def log_archive_audit(audit: pd.DataFrame, print_limit: int) -> None:
    if audit.empty:
        LOGGER.info("Ensembl archive audit: no source Ensembl IDs needed archive fallback.")
        return
    LOGGER.info("Ensembl archive fallback/audit summary:\n%s", archive_summary(audit).to_string(index=False))
    resolved = {
        status.value for status in MappingStatus if status.is_mapped and status != MappingStatus.DIRECT_SOURCE_ENSEMBL
    }
    unresolved = audit[~audit["rescue_status"].isin(resolved)]
    if unresolved.empty:
        LOGGER.info("Ensembl archive audit: all %d problematic source Ensembl mappings were resolved.", len(audit))
        return
    LOGGER.warning("Ensembl archive audit: %d problematic source Ensembl mappings were not resolved.", len(unresolved))
    for row in unresolved.head(print_limit).to_dict("records"):
        LOGGER.warning(
            "Unrescued Ensembl ID example: namespace=%s term=%s source=%s symbol=%s ensembl=%s status=%s",
            row["target_namespace"],
            row["term_id"],
            row["source_tag"],
            row["input_symbol"],
            row["source_ensembl_stable"],
            row["rescue_status"],
        )


def target_namespace_summary(targets: Sequence[Target], memberships: pd.DataFrame) -> pd.DataFrame:
    source_ids = {
        stable for value in pd.unique(memberships["gene_ensembl_from_source"]) if (stable := ensembl_stable_id(value))
    }
    rows = []
    for target in targets:
        helper, target_ids = target.helper_ids, target.effective_ids
        annotated = target_ids & helper
        missing = target_ids - annotated
        rows.append(
            {
                "target_namespace": target.name,
                "output_gene_source": target.output_gene_source,
                "output_genes_size": len(target.output_gene_ids),
                "gene_filter_size": len(target.filter_ids),
                "effective_target_size": len(target_ids),
                "annotation_helper_size": len(helper),
                "primary_annotation_helper_size": len(target.primary_annotation_ids or helper),
                "metadata_fallback_helper_size": max(0, len(helper) - len(target.primary_annotation_ids or helper)),
                "output_genes_with_annotation_metadata": len(annotated),
                "annotation_metadata_coverage": round(len(annotated) / len(target_ids), 6) if target_ids else 0,
                "target_ids_missing_annotation_metadata": len(missing),
                "annotation_ids_absent_from_target": len(helper - target_ids),
                "source_ensembl_ids_kept_without_annotation": len(source_ids & missing),
            }
        )
    return pd.DataFrame(rows, columns=list(TARGET_NAMESPACE_SUMMARY_COLUMNS))


def metadata_fallback_table(targets: Sequence[Target]) -> pd.DataFrame:
    rows = [record for target in targets for record in target.metadata_fallbacks]
    return pd.DataFrame(rows, columns=list(METADATA_FALLBACK_COLUMNS))


def target_gene_id_table(targets: Sequence[Target], kind: str) -> pd.DataFrame:
    """The configured output_genes (``kind="output_genes"``) or gene_filter IDs of every target."""
    rows = []
    for target in targets:
        gene_list = target.output_genes if kind == "output_genes" else target.gene_filter
        table = gene_list.table
        if table.empty and gene_list.ids:
            ids = sorted(gene_list.ids)
            table = pd.DataFrame(
                {"input_gene_id": ids, "ensembl_gene_id": ids, "gene_symbol": "", "id_type": "ensembl_stable"}
            )
        for row in table.to_dict("records"):
            stable = ensembl_stable_id(row.get("ensembl_gene_id", ""))
            if not stable:
                continue
            rows.append(
                {
                    "target_namespace": target.name,
                    "source_path": str(gene_list.source_path or ""),
                    "source_id_column": str(gene_list.id_column or ""),
                    "source_symbol_column": str(gene_list.symbol_column or "") if kind == "output_genes" else "",
                    "input_gene_id": row.get("input_gene_id", ""),
                    "ensembl_gene_id": stable,
                    "gene_symbol": row.get("gene_symbol", ""),
                    "id_type": row.get("id_type", ""),
                    "has_annotation_metadata": stable in target.helper_ids,
                }
            )
    return pd.DataFrame(rows, columns=list(TARGET_GENE_ID_COLUMNS)).drop_duplicates()
