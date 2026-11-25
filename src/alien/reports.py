from __future__ import annotations

import hashlib
import json
import platform
import subprocess
from importlib import metadata
from pathlib import Path
from typing import Any

import pandas as pd

from .dryrun import LOCAL_SOURCE_TYPES
from .outputs import sha256_file
from .sources import (
    ENRICHR_BASE_URL,
    HGNC_FALLBACK_URLS,
    HGNC_URL,
    NCBI_GENE_HISTORY_URL,
    NCBI_GENE_INFO_URL,
    _msigdb_release_info,
    _source_paths,
)
from .utils import ensure_dirs, utc_now, write_json
from .version import __version__

RUNTIME_PACKAGES = ("pandas", "pyarrow", "PyYAML", "rdata", "requests", "urllib3", "pydantic", "filelock")


def build_source_provenance(
    cfg: dict[str, Any],
    source_dir: Path,
    targets: list[dict[str, Any]],
    warnings: list[str],
    *,
    workers: int = 1,
    effective_config_text: str = "",
) -> dict[str, Any]:
    """Record what this build actually read and how, assembled from the run rather than constants."""
    return {
        "tool": "ALIEN",
        "alien_version": __version__,
        "git_describe": _git_describe(),
        "output_schema_version": OUTPUT_SCHEMA_VERSION,
        "created_at": utc_now(),
        "python_version": platform.python_version(),
        "platform": platform.platform(),
        "packages": _python_packages(),
        "config_sha256": hashlib.sha256(effective_config_text.encode("utf-8")).hexdigest(),
        "configured_sources": cfg.get("sources", []),
        "configured_targets": [_target_record(target) for target in targets],
        "input_files": _input_files(cfg, source_dir, targets),
        "resources": _resources(cfg, source_dir),
        "redundancy": {
            "jaccard_cutoff": cfg.get("redundancy", {}).get("jaccard_cutoff"),
            "linkage": cfg.get("redundancy", {}).get("linkage", "single"),
            "exact_duplicate_removal": cfg.get("redundancy", {}).get("exact_duplicate_removal", True),
            "apply_per_family": cfg.get("redundancy", {}).get("apply_per_family", True),
        },
        "workers": workers,
        "warnings": list(warnings),
    }


def _target_record(target: dict[str, Any]) -> dict[str, Any]:
    return {
        "name": target.get("name", ""),
        "type": target.get("type", ""),
        "annotation": target.get("annotation_label", ""),
        "annotation_path": str(target.get("annotation_path") or target.get("annotation_gtf", "")),
        "metadata_fallbacks": target.get("metadata_fallbacks", []),
        "output_genes_path": str(target.get("output_source_path") or ""),
        "output_genes_id_column": str(target.get("output_id_column") or ""),
        "output_genes_symbol_column": str(target.get("output_symbol_column") or ""),
        "gene_filter_path": str(target.get("filter_source_path") or ""),
        "gene_filter_column": str(target.get("filter_source_column") or ""),
    }


def _file_record(role: str, path: Path) -> dict[str, Any]:
    record: dict[str, Any] = {"role": role, "path": str(path)}
    if path.is_file():
        record.update({"sha256": sha256_file(path), "size": path.stat().st_size})
    return record


def _input_files(cfg: dict[str, Any], source_dir: Path, targets: list[dict[str, Any]]) -> list[dict[str, Any]]:
    records: list[dict[str, Any]] = []
    for index, spec in enumerate(cfg.get("sources", [])):
        if spec.get("enabled", True) and spec.get("type") in LOCAL_SOURCE_TYPES:
            records += [_file_record(f"sources[{index}]", path) for path in _source_paths(spec, source_dir)]
    for target in targets:
        name = str(target.get("name", ""))
        if target.get("annotation_path"):
            records.append(_file_record(f"target:{name}:annotation", Path(target["annotation_path"])))
        for fallback in target.get("metadata_fallbacks", []):
            role = f"target:{name}:metadata_fallback_{fallback.get('fallback_index', '')}"
            records.append(_file_record(role, Path(str(fallback.get("annotation_path", "")))))
        for key, role in [("output_source_path", "output_genes"), ("filter_source_path", "gene_filter")]:
            if target.get(key):
                records.append(_file_record(f"target:{name}:{role}", Path(target[key])))
    records.append(_file_record("HGNC", source_dir / "hgnc" / "hgnc_complete_set.txt"))
    if cfg.get("ncbi_gene", {}).get("enabled", False):
        records += [
            _file_record("NCBI_Gene", source_dir / "ncbi_gene" / name)
            for name in ["Homo_sapiens.gene_info.gz", "gene_history.gz"]
        ]
    return records


def _resources(cfg: dict[str, Any], source_dir: Path) -> dict[str, Any]:
    resources: dict[str, Any] = {
        "HGNC": {
            "urls": [HGNC_URL, *HGNC_FALLBACK_URLS],
            "local_file": str(source_dir / "hgnc" / "hgnc_complete_set.txt"),
        }
    }
    ncbi = cfg.get("ncbi_gene", {})
    if ncbi.get("enabled", False):
        resources["NCBI_Gene"] = {
            "gene_info_url": ncbi.get("gene_info_url", NCBI_GENE_INFO_URL),
            "gene_history_url": ncbi.get("gene_history_url", NCBI_GENE_HISTORY_URL),
            "allowed_source_tags": ncbi.get("allowed_source_tags", []),
            "use_gene_info_synonyms": ncbi.get("use_gene_info_synonyms", True),
        }
    archive = cfg.get("ensembl_archive", {})
    if archive.get("enabled", True):
        resources["Ensembl_archive"] = {
            "endpoint": archive.get("endpoint", ""),
            "local_cache_file": str(source_dir / "ensembl_archive" / "archive_id_cache.json"),
        }
    sources = [spec for spec in cfg.get("sources", []) if spec.get("enabled", True)]
    msigdb = [spec for spec in sources if spec.get("type") == "msigdb_remote"]
    if msigdb:
        resources["MSigDB"] = [
            {
                "version": str(spec.get("version", "")),
                "db_species": str(spec.get("db_species", "HS")),
                "collection": str(spec.get("collection", "")),
                "subcollection": str(spec.get("subcollection", "")),
                **{key: _msigdb_release_info(spec)[key] for key in ("zip_url", "zip_md5", "zip_name")},
            }
            for spec in msigdb
        ]
    if any(spec.get("type") == "enrichr_remote" for spec in sources):
        resources["Enrichr"] = {"source_url": ENRICHR_BASE_URL}
    return resources


def _git_describe() -> str:
    root = Path(__file__).resolve().parents[2]
    if not (root / ".git").exists():
        return ""
    try:
        result = subprocess.run(
            ["git", "describe", "--always", "--dirty", "--tags"],
            cwd=root,
            capture_output=True,
            text=True,
            check=False,
            timeout=5,
        )
    except (OSError, subprocess.SubprocessError):
        return ""
    return result.stdout.strip()


def write_provenance(obj: dict[str, Any], path: Path) -> None:
    write_json(obj, path)


# Bumped whenever an output file or column changes incompatibly.
OUTPUT_SCHEMA_VERSION = 2

TERM_MANIFEST_COLUMNS = [
    "term_id",
    "target_namespace",
    "original_name",
    "display_name",
    "family",
    "aspect",
    "source",
    "source_tag",
    "collection",
    "subcollection",
    "db_version",
    "source_url",
    "source_license_note",
    "original_size_symbols",
    "size_after_mapping",
    "kept",
    "removed_reason",
    "representative_term_id",
    "redundancy_cluster_id",
]


def build_term_manifest(
    source_terms: pd.DataFrame,
    mapped_pre_filter: dict[str, dict[str, dict[str, Any]]],
    final_terms: dict[str, dict[str, dict[str, Any]]],
    removed_size: pd.DataFrame,
    removed_redundancy: pd.DataFrame,
) -> pd.DataFrame:
    """One row per (term, target namespace) with the namespace-specific outcome."""
    if source_terms.empty:
        return pd.DataFrame(columns=TERM_MANIFEST_COLUMNS)
    term_meta = source_terms.drop_duplicates("term_id").set_index("term_id")
    original_sizes = source_terms.groupby("term_id")["gene_symbol"].nunique().to_dict()
    removed = _removed_reason_map(removed_size, removed_redundancy)
    rows: list[dict[str, object]] = []
    for term_id in sorted(term_meta.index):
        meta = term_meta.loc[term_id].to_dict()
        for namespace, mapped in mapped_pre_filter.items():
            removal = removed.get((namespace, term_id), {})
            kept = term_id in final_terms.get(namespace, {})
            reason = str(removal.get("reason", ""))
            if not kept and not reason:
                reason = "no_mapped_genes"
            rows.append(
                {
                    "term_id": term_id,
                    "target_namespace": namespace,
                    **{column: meta.get(column, "") for column in TERM_MANIFEST_COLUMNS[2:13]},
                    "original_size_symbols": original_sizes.get(term_id, 0),
                    "size_after_mapping": len(mapped[term_id]["genes"]) if term_id in mapped else 0,
                    "kept": kept,
                    "removed_reason": reason,
                    "representative_term_id": removal.get("representative_term_id", "") or term_id,
                    "redundancy_cluster_id": removal.get("redundancy_cluster_id", ""),
                }
            )
    return pd.DataFrame(rows, columns=TERM_MANIFEST_COLUMNS)


def build_source_manifest(source_terms: pd.DataFrame) -> pd.DataFrame:
    """Summarize the resolved source collections used in a build."""
    columns = [
        "source",
        "source_tag",
        "collection",
        "subcollection",
        "family",
        "aspect",
        "db_version",
        "selected_library",
        "configured_name",
        "match",
        "match_method",
        "candidate_libraries",
        "n_terms",
        "n_memberships",
        "n_unique_symbols",
        "source_url",
        "source_license_note",
    ]
    if source_terms.empty:
        return pd.DataFrame(columns=columns)

    frame = source_terms.copy()
    group_columns = ["source", "source_tag", "collection", "subcollection", "family", "aspect", "db_version"]
    for column in group_columns + ["source_url", "source_license_note", "metadata_json", "term_id", "gene_symbol"]:
        if column not in frame:
            frame[column] = ""

    rows: list[dict[str, object]] = []
    grouped = frame.groupby(group_columns, dropna=False, sort=True)
    for key, group in grouped:
        record = dict(zip(group_columns, key))
        metadata = _source_group_metadata(group)
        rows.append(
            {
                **record,
                "selected_library": metadata.get("selected_library", ""),
                "configured_name": metadata.get("configured_name", ""),
                "match": metadata.get("match", ""),
                "match_method": metadata.get("match_method", ""),
                "candidate_libraries": metadata.get("candidate_libraries", ""),
                "n_terms": int(group["term_id"].nunique()),
                "n_memberships": len(group),
                "n_unique_symbols": int(group["gene_symbol"].nunique()),
                "source_url": _single_value(group["source_url"]),
                "source_license_note": _single_value(group["source_license_note"]),
            }
        )
    return pd.DataFrame(rows, columns=columns)


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
            too_small = _removed_count(removed_size, namespace, family, source_tag, "too_small")
            too_large = _removed_count(removed_size, namespace, family, source_tag, "too_large")
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


def _removed_reason_map(*frames: pd.DataFrame) -> dict[tuple[str, str], dict[str, object]]:
    out: dict[tuple[str, str], dict[str, object]] = {}
    for frame in frames:
        if frame is None or frame.empty:
            continue
        for row in frame.to_dict("records"):
            out.setdefault((str(row["target_namespace"]), str(row["term_id"])), row)
    return out


def _group_terms(terms: dict[str, dict[str, Any]]) -> dict[tuple[str, str], dict[str, dict[str, Any]]]:
    grouped: dict[tuple[str, str], dict[str, dict[str, Any]]] = {}
    for term_id, record in terms.items():
        meta = record["meta"]
        key = (str(meta.get("family", "")), str(meta.get("source_tag", "")))
        grouped.setdefault(key, {})[term_id] = record
    return grouped


def _source_group_metadata(group: pd.DataFrame) -> dict[str, str]:
    parsed = [
        _parse_metadata_json(value) for value in group.get("metadata_json", pd.Series(dtype=str)).drop_duplicates()
    ]
    return {
        "selected_library": _metadata_value(parsed, "selected_library"),
        "configured_name": _metadata_value(parsed, "configured_name"),
        "match": _metadata_value(parsed, "match"),
        "match_method": _metadata_value(parsed, "match_method"),
        "candidate_libraries": _metadata_list_value(parsed, "candidate_libraries"),
    }


def _parse_metadata_json(value: object) -> dict[str, Any]:
    text = str(value or "").strip()
    if not text:
        return {}
    try:
        parsed = json.loads(text)
    except json.JSONDecodeError:
        return {}
    return parsed if isinstance(parsed, dict) else {}


def _metadata_value(records: list[dict[str, Any]], key: str) -> str:
    values = sorted({str(record.get(key, "")).strip() for record in records if str(record.get(key, "")).strip()})
    return ";".join(values)


def _metadata_list_value(records: list[dict[str, Any]], key: str) -> str:
    values: set[str] = set()
    for record in records:
        raw = record.get(key, [])
        if isinstance(raw, list):
            values.update(str(value).strip() for value in raw if str(value).strip())
        elif str(raw or "").strip():
            values.add(str(raw).strip())
    return json.dumps(sorted(values), ensure_ascii=True) if values else ""


def _single_value(series: pd.Series) -> str:
    values = sorted({str(value).strip() for value in series if str(value).strip()})
    return values[0] if len(values) == 1 else ";".join(values)


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


def _removed_count(df: pd.DataFrame, namespace: str, family: str, source_tag: str, reason: str) -> int:
    if df.empty:
        return 0
    mask = (
        df["target_namespace"].eq(namespace)
        & df["family"].eq(family)
        & df["source_tag"].eq(source_tag)
        & df["reason"].eq(reason)
    )
    return int(mask.sum())


def _python_packages() -> dict[str, str]:
    packages: dict[str, str] = {}
    for package in RUNTIME_PACKAGES:
        try:
            packages[package] = metadata.version(package)
        except metadata.PackageNotFoundError:
            packages[package] = "not installed"
    return packages
