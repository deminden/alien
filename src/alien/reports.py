from __future__ import annotations

import platform
import subprocess
from pathlib import Path
from typing import Any

import pandas as pd

from .utils import ensure_dirs, utc_now, write_json


def build_source_provenance(
    cfg: dict[str, Any],
    source_dir: Path,
    targets: list[dict[str, Any]],
    warnings: list[str],
) -> dict[str, Any]:
    # Record source versions and local cache paths.
    msigdb_dir = source_dir / "msigdb"
    db_version_path = msigdb_dir / "db_version.txt"
    db_version = db_version_path.read_text(encoding="utf-8").strip() if db_version_path.exists() else ""
    target_records = [
        {
            "name": target.get("name", ""),
            "type": target.get("type", ""),
            "annotation": target.get("annotation_label", ""),
            "annotation_path": str(target.get("annotation_path") or target.get("annotation_gtf", "")),
            "restrict_to": str(target.get("source_path") or ""),
        }
        for target in targets
    ]
    return {
        "created_at": utc_now(),
        "tool": "ALIEN",
        "python_version": platform.python_version(),
        "r_version": _r_version(),
        "packages": _python_packages(),
        "configured_sources": cfg.get("sources", []),
        "configured_targets": target_records,
        "sources": {
            "MSigDB": {
                "target_release": "2026.1.Hs",
                "db_version": db_version,
                "source_url": "https://www.gsea-msigdb.org/gsea/msigdb",
                "local_raw_dir": str(msigdb_dir),
                "license_note": "MSigDB license and attribution terms apply.",
            },
            "HGNC": {
                "source_url": "https://ftp.ebi.ac.uk/pub/databases/genenames/hgnc/tsv/hgnc_complete_set.txt",
                "local_raw_file": str(source_dir / "hgnc" / "hgnc_complete_set.txt"),
                "license_note": "HGNC data license and attribution terms apply.",
            },
            "NCBI_Gene": {
                "gene_info_url": "https://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/Homo_sapiens.gene_info.gz",
                "gene_history_url": "https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene_history.gz",
                "local_cache_file": str(source_dir / "ncbi_gene" / "ncbi_gene_symbol_rescue.tsv.gz"),
                "license_note": "NCBI Gene history and gene info are used only as a lower-confidence audited rescue layer.",
            },
            "GENCODE_v47": {
                "source_url": "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_47/gencode.v47.annotation.gtf.gz",
                "license_note": "GENCODE release 47; used as a symbol/metadata helper while configured target IDs define the GMT namespace.",
            },
            "GENCODE_v29": {
                "source_url": "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_29/gencode.v29.annotation.gtf.gz",
                "license_note": "GENCODE release 29; used as a symbol/metadata helper while configured target IDs define the GMT namespace.",
            },
            "Ensembl_archive": {
                "source_url": "https://rest.ensembl.org/archive/id",
                "local_cache_file": str(source_dir / "ensembl_archive" / "archive_id_cache.json"),
                "license_note": "Ensembl REST archive stable-ID lookup; cached for reproducible reruns.",
            },
        },
        "reference_facts": {
            "msigdb_target_release": "2026.1.Hs",
            "min_msigdbr_version": "26.1.0",
            "official_v0_1_scope": "human gene sets with HGNC symbols and Ensembl target namespaces",
        },
        "warnings": warnings,
    }


def write_provenance(obj: dict[str, Any], path: Path) -> None:
    write_json(obj, path)


def build_term_manifest(
    source_terms: pd.DataFrame,
    final_terms: dict[str, dict[str, dict[str, Any]]],
    removed_size: pd.DataFrame,
    removed_redundancy: pd.DataFrame,
) -> pd.DataFrame:
    # Summarize term status across all output namespaces
    rows: list[dict[str, object]] = []
    if source_terms.empty:
        return pd.DataFrame()
    term_meta = source_terms.drop_duplicates("term_id").set_index("term_id")
    original_sizes = source_terms.groupby("term_id")["gene_symbol"].nunique().to_dict()
    final_sizes = {
        namespace: {term_id: len(record["genes"]) for term_id, record in terms.items()}
        for namespace, terms in final_terms.items()
    }
    removed_reason = _removed_reason_map(removed_size, removed_redundancy)
    for term_id in sorted(term_meta.index):
        row = term_meta.loc[term_id].to_dict()
        representative = removed_reason.get(term_id, {}).get("representative_term_id", term_id)
        manifest_row = {
            "term_id": term_id,
            "original_name": row.get("original_name", ""),
            "display_name": row.get("display_name", ""),
            "family": row.get("family", ""),
            "aspect": row.get("aspect", ""),
            "source": row.get("source", ""),
            "source_tag": row.get("source_tag", ""),
            "collection": row.get("collection", ""),
            "subcollection": row.get("subcollection", ""),
            "db_version": row.get("db_version", ""),
            "source_url": row.get("source_url", ""),
            "source_license_note": row.get("source_license_note", ""),
            "original_size_symbols": original_sizes.get(term_id, 0),
            "representative_term_id": representative,
            "is_representative": representative == term_id,
            "redundancy_cluster_id": removed_reason.get(term_id, {}).get("redundancy_cluster_id", ""),
            "removed_reason": removed_reason.get(term_id, {}).get("reason", ""),
            "created_at": utc_now(),
        }
        for namespace, sizes in final_sizes.items():
            safe_namespace = str(namespace).replace("-", "_").replace(".", "_")
            manifest_row[f"{safe_namespace}_size_after_mapping"] = sizes.get(term_id, 0)
            manifest_row[f"kept_{safe_namespace}"] = term_id in sizes
        rows.append(manifest_row)
    return pd.DataFrame(rows)


def collection_summary(
    pre_filter: dict[str, dict[str, dict[str, Any]]],
    post_size: dict[str, dict[str, dict[str, Any]]],
    final: dict[str, dict[str, dict[str, Any]]],
) -> pd.DataFrame:
    # Summarize term counts after each filtering stage
    rows: list[dict[str, object]] = []
    for namespace, terms in final.items():
        final_groups = _group_terms(terms)
        size_groups = _group_terms(post_size.get(namespace, {}))
        pre_groups = _group_terms(pre_filter.get(namespace, {}))
        keys = sorted(set(pre_groups) | set(size_groups) | set(final_groups))
        for family, source_tag in keys:
            subset = final_groups.get((family, source_tag), {})
            sizes = [len(record["genes"]) for record in subset.values()]
            genes = sorted({gene for record in subset.values() for gene in record["genes"]})
            rows.append(
                {
                    "target_namespace": namespace,
                    "family": family,
                    "source_tag": source_tag,
                    "n_terms_before_filter": len(pre_groups.get((family, source_tag), {})),
                    "n_terms_after_size_filter": len(size_groups.get((family, source_tag), {})),
                    "n_terms_after_redundancy": len(subset),
                    "median_size": float(pd.Series(sizes).median()) if sizes else 0,
                    "min_size": min(sizes) if sizes else 0,
                    "max_size": max(sizes) if sizes else 0,
                    "n_unique_genes": len(genes),
                }
            )
    return pd.DataFrame(rows)


def mapping_summary(
    original: pd.DataFrame,
    mapped_pre_filter: dict[str, dict[str, dict[str, Any]]],
    unmapped: pd.DataFrame,
    ambiguous: pd.DataFrame,
    non_gene: pd.DataFrame,
    removed_size: pd.DataFrame,
    mapping_status: pd.DataFrame | None = None,
) -> pd.DataFrame:
    # Summarize raw mapping events and collapsed unique genes
    rows: list[dict[str, object]] = []
    original_counts = original.groupby(["family", "source_tag"]).size().to_dict() if not original.empty else {}
    for namespace, terms in mapped_pre_filter.items():
        grouped = _group_terms(terms)
        for family, source_tag in sorted(original_counts):
            n_original = int(original_counts.get((family, source_tag), 0))
            n_mapped_unique = sum(len(record["genes"]) for record in grouped.get((family, source_tag), {}).values())
            n_mapped_raw = _mapped_status_count(mapping_status, namespace, family, source_tag)
            n_unmapped = _count(unmapped, namespace, family, source_tag)
            n_ambiguous = _count(ambiguous, namespace, family, source_tag)
            n_non_gene = _count(non_gene, namespace, family, source_tag)
            n_considered = max(0, n_original - n_non_gene)
            if namespace == "symbols":
                n_mapped_raw = n_considered
            too_small = _removed_count(removed_size, namespace, family, "too_small")
            too_large = _removed_count(removed_size, namespace, family, "too_large")
            rows.append(
                {
                    "target_namespace": namespace,
                    "family": family,
                    "source_tag": source_tag,
                    "n_gene_members_original": n_original,
                    "n_non_gene_token_dropped": n_non_gene,
                    "n_gene_members_considered_for_mapping": n_considered,
                    "n_gene_members_mapped_raw": n_mapped_raw,
                    "n_gene_members_mapped_unique": n_mapped_unique,
                    "n_gene_members_collapsed_after_mapping": max(0, n_mapped_raw - n_mapped_unique),
                    "n_gene_members_unmapped": n_unmapped,
                    "n_ambiguous": n_ambiguous,
                    "mapping_rate_raw": round(n_mapped_raw / n_considered, 6) if n_considered else 0,
                    "mapping_rate_unique": round(n_mapped_unique / n_considered, 6) if n_considered else 0,
                    "n_terms_lost_because_too_small": too_small,
                    "n_terms_lost_because_too_large": too_large,
                }
            )
    return pd.DataFrame(rows)


def write_qc(
    qc_dir: Path,
    collection: pd.DataFrame,
    mapping: pd.DataFrame,
    redundancy: pd.DataFrame,
    warnings: list[str],
) -> None:
    # Write QC tables and warnings
    ensure_dirs(qc_dir)
    collection.to_csv(qc_dir / "collection_summary.tsv", sep="\t", index=False)
    write_json(collection.to_dict(orient="records"), qc_dir / "collection_summary.json")
    mapping.to_csv(qc_dir / "mapping_summary.tsv", sep="\t", index=False)
    redundancy.to_csv(qc_dir / "redundancy_summary.tsv", sep="\t", index=False)
    (qc_dir / "warnings.txt").write_text("\n".join(warnings) + ("\n" if warnings else ""), encoding="utf-8")


def _removed_reason_map(*frames: pd.DataFrame) -> dict[str, dict[str, object]]:
    out: dict[str, dict[str, object]] = {}
    for frame in frames:
        if frame is None or frame.empty:
            continue
        for _, row in frame.iterrows():
            out.setdefault(str(row["term_id"]), row.to_dict())
    return out


def _group_terms(terms: dict[str, dict[str, Any]]) -> dict[tuple[str, str], dict[str, dict[str, Any]]]:
    grouped: dict[tuple[str, str], dict[str, dict[str, Any]]] = {}
    for term_id, record in terms.items():
        meta = record["meta"]
        key = (str(meta.get("family", "")), str(meta.get("source_tag", "")))
        grouped.setdefault(key, {})[term_id] = record
    return grouped


def _count(df: pd.DataFrame, namespace: str, family: str, source_tag: str) -> int:
    if df.empty:
        return 0
    mask = df["target_namespace"].eq(namespace) & df["family"].eq(family) & df["source_tag"].eq(source_tag)
    return int(mask.sum())


def _mapped_status_count(df: pd.DataFrame | None, namespace: str, family: str, source_tag: str) -> int:
    if df is None or df.empty:
        return 0
    mapped_statuses = {
        "direct_source_ensembl",
        "archive_ensembl_rescue",
        "hgnc_current_symbol",
        "hgnc_previous_symbol",
        "hgnc_alias_symbol",
        "hgnc_ensembl_tiebreak",
        "manual_symbol_repair",
        "ncbi_gene_history",
        "ncbi_gene_info_synonym",
    }
    mask = (
        df["target_namespace"].eq(namespace)
        & df["family"].eq(family)
        & df["source_tag"].eq(source_tag)
        & df["mapping_status"].isin(mapped_statuses)
    )
    return int(df.loc[mask, "n_gene_members"].sum())


def _removed_count(df: pd.DataFrame, namespace: str, family: str, reason: str) -> int:
    if df.empty:
        return 0
    mask = df["target_namespace"].eq(namespace) & df["family"].eq(family) & df["reason"].eq(reason)
    return int(mask.sum())


def _r_version() -> str:
    try:
        result = subprocess.run(["Rscript", "--version"], check=False, capture_output=True, text=True)
        return (result.stderr or result.stdout).strip()
    except FileNotFoundError:
        return "Rscript unavailable"


def _python_packages() -> dict[str, str]:
    packages = {}
    for package in ["pandas", "pyarrow", "yaml", "requests", "networkx"]:
        try:
            mod = __import__(package)
            packages[package] = getattr(mod, "__version__", "unknown")
        except Exception:
            packages[package] = "unavailable"
    return packages
