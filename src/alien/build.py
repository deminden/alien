from __future__ import annotations

from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import pandas as pd

from .config import load_config, merge_config
from .filtering import apply_size_and_broad_filters, remove_redundancy
from .mapping import map_namespaces
from .reports import (
    build_source_provenance,
    build_term_manifest,
    collection_summary,
    mapping_summary,
    write_provenance,
    write_qc,
)
from .sources import (
    load_ensembl_gtf_target,
    load_hgnc,
    load_ncbi_gene_maps,
    load_source_memberships,
    mark_universe,
    read_gene_universe,
    write_gene_symbol_master,
    write_mapping,
)
from .utils import LOGGER, ensure_dirs, process_pool_context, setup_logging, strip_ensembl_version, write_gmt, write_tsv_gz


DEFAULT_OUTDIR = Path("data/alien_gmt")
TERM_ID_IDENTITY_COLUMNS = [
    "original_name",
    "display_name",
    "description",
    "source",
    "source_tag",
    "collection",
    "subcollection",
    "family",
    "aspect",
    "db_version",
    "source_url",
    "source_license_note",
    "metadata_json",
]


@dataclass(frozen=True)
class BuildResult:
    outdir: Path
    namespaces: tuple[str, ...]
    n_source_memberships: int
    n_source_terms: int
    warnings: tuple[str, ...]


def build(
    config: str | Path | dict[str, Any] | None = None,
    outdir: str | Path | None = None,
    workers: int | None = None,
    dry_run: bool = False,
    force_download: bool = False,
) -> BuildResult:
    setup_logging()
    cfg = _coerce_config(config)
    if workers is not None:
        cfg = merge_config(cfg, {"runtime": {"workers": workers}})
    worker_count = _worker_count(cfg)
    output_dir = Path(outdir or cfg.get("project", {}).get("outdir", DEFAULT_OUTDIR))
    source_dir = Path(cfg.get("project", {}).get("source_dir", "data/alien_sources"))
    metadata_dir = output_dir / "metadata"
    qc_dir = output_dir / "qc"
    gmt_dir = output_dir / "gmt"

    if dry_run:
        LOGGER.info("Dry run: would build ALIEN GMT outputs under %s with %d worker process(es).", output_dir, worker_count)
        return BuildResult(output_dir, tuple(_configured_namespace_names(cfg)), 0, 0, tuple())

    ensure_dirs(metadata_dir, qc_dir, gmt_dir)
    warnings: list[str] = []
    LOGGER.info("Starting ALIEN GMT build")
    LOGGER.info("Output directory: %s", output_dir)
    LOGGER.info("Source cache directory: %s", source_dir)
    LOGGER.info("Local worker processes: %d", worker_count)

    source_terms = load_source_memberships(cfg, force_download=force_download)
    if source_terms.empty:
        raise RuntimeError("No source gene-set memberships were loaded.")
    LOGGER.info("Loaded %d source gene memberships across %d terms", len(source_terms), source_terms["term_id"].nunique())
    _check_term_id_collisions(source_terms, metadata_dir, cfg, warnings)

    LOGGER.info("Loading mapping resources")
    hgnc_maps = load_hgnc(source_dir, force=force_download)
    ncbi_maps = load_ncbi_gene_maps(cfg.get("ncbi_gene", {}), source_dir, force=force_download)
    targets = _prepare_targets(cfg, source_dir, metadata_dir, source_terms, force_download, warnings)
    write_gene_symbol_master(hgnc_maps, metadata_dir / "gene_symbol_master.tsv.gz")

    LOGGER.info("Projecting gene sets to configured namespaces")
    mapped_pre_filter, unmapped, ambiguous, non_gene, ensembl_archive_audit, ncbi_gene_audit, mapping_status_summary = map_namespaces(
        source_terms, hgnc_maps, ncbi_maps, targets, cfg, source_dir, workers=worker_count
    )
    _report_ensembl_archive_audit(ensembl_archive_audit, cfg)

    LOGGER.info("Applying size and redundancy filters")
    size_filtered, filtered, removed_size, removed_redundancy, redundancy_summary = _filter_namespaces(
        mapped_pre_filter, cfg, worker_count
    )

    LOGGER.info("Writing combined GMT outputs")
    for namespace, terms in filtered.items():
        write_gmt(terms, gmt_dir / f"{namespace}.gmt")
        if not terms:
            warnings.append(f"Empty output for namespace {namespace}.")

    LOGGER.info("Writing metadata and QC reports")
    write_tsv_gz(build_term_manifest(source_terms, filtered, removed_size, removed_redundancy), metadata_dir / "term_manifest.tsv.gz")
    write_tsv_gz(_target_gene_restrictions(targets), metadata_dir / "target_gene_restrictions.tsv.gz")
    write_tsv_gz(removed_size, metadata_dir / "removed_terms_size_filter.tsv.gz")
    write_tsv_gz(removed_redundancy, metadata_dir / "removed_terms_redundancy.tsv.gz")
    write_tsv_gz(_with_columns(unmapped, ["target_namespace", "family", "term_id", "source_tag", "input_gene", "attempted_current_symbol", "reason", "possible_matches"]), metadata_dir / "unmapped_genes.tsv.gz")
    write_tsv_gz(_with_columns(non_gene, ["target_namespace", "family", "term_id", "source_tag", "input_gene", "reason"]), metadata_dir / "dropped_non_gene_tokens.tsv.gz")
    write_tsv_gz(
        _with_columns(
            ambiguous,
            ["target_namespace", "family", "term_id", "source_tag", "input_gene", "mapping_stage", "possible_symbols", "possible_ensembl_ids", "action"],
        ),
        metadata_dir / "ambiguous_gene_mappings.tsv.gz",
    )
    write_tsv_gz(
        _with_columns(
            ensembl_archive_audit,
            [
                "target_namespace",
                "family",
                "term_id",
                "source_tag",
                "input_symbol",
                "source_ensembl_id",
                "source_ensembl_stable",
                "archive_found",
                "archive_is_current",
                "archive_release",
                "archive_assembly",
                "archive_latest",
                "archive_possible_replacements",
                "archive_error",
                "rescued_ensembl_id",
                "rescue_status",
            ],
        ),
        metadata_dir / "ensembl_archive_audit.tsv.gz",
    )
    write_tsv_gz(
        _with_columns(
            ncbi_gene_audit,
            [
                "target_namespace",
                "family",
                "term_id",
                "source_tag",
                "input_gene",
                "ncbi_stage",
                "ncbi_status",
                "ncbi_current_symbol",
                "ncbi_gene_ids",
                "ncbi_ensembl_ids",
                "possible_symbols",
                "possible_ensembl_ids",
                "resolved_ensembl_id",
                "action",
            ],
        ),
        metadata_dir / "ncbi_gene_rescue_audit.tsv.gz",
    )
    _ensembl_archive_summary(ensembl_archive_audit).to_csv(qc_dir / "ensembl_archive_summary.tsv", sep="\t", index=False)
    _non_gene_summary(non_gene).to_csv(qc_dir / "non_gene_token_summary.tsv", sep="\t", index=False)
    _ncbi_gene_summary(ncbi_gene_audit).to_csv(qc_dir / "ncbi_gene_rescue_summary.tsv", sep="\t", index=False)
    _target_universe_summary(targets, source_terms).to_csv(qc_dir / "target_universe_summary.tsv", sep="\t", index=False)
    mapping_status_summary.to_csv(qc_dir / "mapping_status_summary.tsv", sep="\t", index=False)
    provenance = build_source_provenance(cfg, source_dir, targets, warnings)
    write_provenance(provenance, metadata_dir / "source_provenance.json")
    write_qc(
        qc_dir,
        collection_summary(mapped_pre_filter, size_filtered, filtered),
        mapping_summary(source_terms, mapped_pre_filter, unmapped, ambiguous, non_gene, removed_size, mapping_status_summary),
        redundancy_summary,
        warnings,
    )

    namespaces = tuple(filtered)
    LOGGER.info("Done. Wrote ALIEN GMT build outputs under %s.", output_dir)
    return BuildResult(
        outdir=output_dir,
        namespaces=namespaces,
        n_source_memberships=len(source_terms),
        n_source_terms=int(source_terms["term_id"].nunique()),
        warnings=tuple(warnings),
    )


def _coerce_config(config: str | Path | dict[str, Any] | None) -> dict[str, Any]:
    if config is None or isinstance(config, (str, Path)):
        return load_config(config)
    return merge_config(load_config(None), config)


def _configured_namespace_names(cfg: dict[str, Any]) -> list[str]:
    names = ["symbols"] if cfg.get("outputs", {}).get("include_symbols", True) else []
    names.extend(str(target["name"]) for target in cfg.get("targets", []))
    return names


def _worker_count(cfg: dict[str, Any]) -> int:
    workers = int(cfg.get("runtime", {}).get("workers", 1) or 1)
    if workers < 1:
        raise ValueError("workers must be a positive integer")
    return workers


def _check_term_id_collisions(
    source_terms: pd.DataFrame,
    metadata_dir: Path,
    cfg: dict[str, Any],
    warnings: list[str],
) -> None:
    audit = _term_id_collision_audit(source_terms)
    if audit.empty:
        return
    options = cfg.get("term_id_collisions", {})
    if isinstance(options, str):
        action = options
    elif isinstance(options, dict):
        action = str(options.get("action", "error"))
    else:
        raise ValueError("term_id_collisions must be a mapping or one of: error, warn, merge")
    action = action.lower().strip()
    if action not in {"error", "warn", "merge"}:
        raise ValueError("term_id_collisions.action must be one of: error, warn, merge")
    audit = audit.assign(action=action)
    audit_path = metadata_dir / "term_id_collisions.tsv"
    ensure_dirs(audit_path.parent)
    audit.to_csv(audit_path, sep="\t", index=False)
    n_terms = int(audit["term_id"].nunique())
    n_signatures = len(audit)
    message = (
        f"Detected {n_signatures} conflicting term/source identities across {n_terms} term_id value(s); "
        f"see {audit_path}."
    )
    if action == "error":
        raise RuntimeError(f"{message} Rename or prefix the colliding term IDs, or set term_id_collisions.action explicitly.")
    warnings.append(f"{message} Continuing because term_id_collisions.action is {action!r}.")
    LOGGER.warning("%s Continuing because term_id_collisions.action is %r.", message, action)


def _term_id_collision_audit(source_terms: pd.DataFrame) -> pd.DataFrame:
    columns = _term_identity_columns(source_terms)
    audit_columns = [
        "term_id",
        "n_identities_for_term",
        "collision_identity_id",
        "n_rows",
        "n_gene_symbols",
        "sample_gene_symbols",
        *columns,
    ]
    if source_terms.empty or "term_id" not in source_terms:
        return pd.DataFrame(columns=audit_columns)

    working = source_terms.copy()
    for column in ["term_id", "gene_symbol", *columns]:
        if column not in working:
            working[column] = ""
        working[column] = working[column].fillna("").astype(str)

    identity = working[["term_id", *columns]].drop_duplicates()
    counts = identity.groupby("term_id", sort=False).size().rename("n_identities_for_term").reset_index()
    collided_ids = set(counts.loc[counts["n_identities_for_term"] > 1, "term_id"])
    if not collided_ids:
        return pd.DataFrame(columns=audit_columns)

    collided = working[working["term_id"].isin(collided_ids)]
    grouped = (
        collided.groupby(["term_id", *columns], dropna=False, sort=False)
        .agg(
            n_rows=("gene_symbol", "size"),
            n_gene_symbols=("gene_symbol", _n_unique_nonempty),
            sample_gene_symbols=("gene_symbol", _sample_values),
        )
        .reset_index()
    )
    grouped = grouped.merge(counts, on="term_id", how="left")
    grouped = grouped.sort_values(["term_id", *columns], kind="mergesort").reset_index(drop=True)
    grouped["collision_identity_id"] = grouped.groupby("term_id").cumcount() + 1
    return grouped[audit_columns]


def _term_identity_columns(source_terms: pd.DataFrame) -> list[str]:
    return [column for column in TERM_ID_IDENTITY_COLUMNS if column in source_terms.columns]


def _n_unique_nonempty(values: pd.Series) -> int:
    return len({str(value) for value in values if str(value)})


def _sample_values(values: pd.Series, limit: int = 12) -> str:
    unique = sorted({str(value) for value in values if str(value)})
    return ";".join(unique[:limit])


def _prepare_targets(
    cfg: dict[str, Any],
    source_dir: Path,
    metadata_dir: Path,
    source_terms: pd.DataFrame,
    force_download: bool,
    warnings: list[str],
) -> list[dict[str, Any]]:
    prepared: list[dict[str, Any]] = []
    for target in cfg.get("targets", []):
        name = str(target["name"])
        target_type = str(target.get("type", "ensembl_gtf"))
        if target_type != "ensembl_gtf":
            raise ValueError(f"Target {name} has unsupported type {target_type!r}.")
        restriction = target.get("restrict_to")
        if _empty_config_value(restriction):
            restriction = target.get("gene_universe")
        universe_path, universe_column, universe_ids = _target_restriction(restriction)
        universe, universe_table, id_type = read_gene_universe(universe_path, column=universe_column, ids=universe_ids)
        if universe is None:
            LOGGER.info("No dataset restriction supplied for %s; using all annotation genes.", name)
        LOGGER.info("%s universe ID type: %s", name, id_type)
        annotation, annotation_label, annotation_path = load_ensembl_gtf_target(source_dir, target, force_download)
        annotation = mark_universe(annotation, universe)
        write_mapping(annotation, metadata_dir / f"gene_mapping_{name}.tsv.gz")
        prepared.append(
            {
                **target,
                "name": name,
                "type": target_type,
                "source_path": universe_path,
                "source_column": universe_column,
                "universe": universe,
                "universe_table": universe_table,
                "annotation": annotation,
                "annotation_path": annotation_path,
                "annotation_label": annotation_label,
            }
        )
    return prepared


def _target_restriction(value: object) -> tuple[Path | None, str | None, list[str] | tuple[str, ...] | set[str] | None]:
    if _empty_config_value(value):
        return None, None, None
    if isinstance(value, dict):
        path = _optional_path(value.get("path") or value.get("file"))
        column = _optional_string(value.get("column") or value.get("gene_column") or value.get("field"))
        ids = value.get("ids", value.get("genes"))
        if ids is not None:
            if isinstance(ids, str):
                ids = [ids]
            return path, column, ids
        if path is None:
            raise ValueError("Target restrict_to dictionaries must define either path or ids.")
        return path, column, None
    if isinstance(value, (list, tuple, set)):
        return None, None, value
    return _optional_path(value), None, None


def _optional_string(value: object) -> str | None:
    if _empty_config_value(value):
        return None
    return str(value)


def _optional_path(value: object) -> Path | None:
    if _empty_config_value(value):
        return None
    return Path(str(value))


def _empty_config_value(value: object) -> bool:
    return value is None or (isinstance(value, str) and value.strip() in {"", "null"})


def _filter_namespaces(
    mapped_pre_filter: dict[str, dict[str, dict[str, Any]]],
    cfg: dict,
    workers: int,
) -> tuple[dict[str, dict], dict[str, dict], pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    if workers <= 1:
        return _filter_namespaces_serial(mapped_pre_filter, cfg)
    return _filter_namespaces_parallel(mapped_pre_filter, cfg, workers)


def _filter_namespaces_serial(
    mapped_pre_filter: dict[str, dict[str, dict[str, Any]]],
    cfg: dict,
) -> tuple[dict[str, dict], dict[str, dict], pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    filtered: dict[str, dict] = {}
    size_filtered: dict[str, dict] = {}
    removed_size_frames: list[pd.DataFrame] = []
    removed_redundancy_frames: list[pd.DataFrame] = []
    redundancy_summary_frames: list[pd.DataFrame] = []
    for namespace, terms in mapped_pre_filter.items():
        LOGGER.info("Filtering namespace %s with %d projected terms", namespace, len(terms))
        size_kept, removed_size = apply_size_and_broad_filters(terms, namespace, cfg)
        final_terms, removed_redundancy, redundancy_summary = remove_redundancy(size_kept, namespace, cfg)
        size_filtered[namespace] = size_kept
        filtered[namespace] = final_terms
        removed_size_frames.append(removed_size)
        removed_redundancy_frames.append(removed_redundancy)
        redundancy_summary_frames.append(redundancy_summary)
    return (
        size_filtered,
        filtered,
        _concat_any(removed_size_frames),
        _concat_any(removed_redundancy_frames),
        _concat_any(redundancy_summary_frames),
    )


def _filter_namespaces_parallel(
    mapped_pre_filter: dict[str, dict[str, dict[str, Any]]],
    cfg: dict,
    workers: int,
) -> tuple[dict[str, dict], dict[str, dict], pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    namespace_order = list(mapped_pre_filter)
    if not namespace_order:
        return {}, {}, pd.DataFrame(), pd.DataFrame(), pd.DataFrame()
    LOGGER.info("Applying size filters with %d workers", min(workers, len(namespace_order)))
    size_results: dict[str, tuple[dict[str, dict[str, Any]], pd.DataFrame]] = {}
    with ProcessPoolExecutor(max_workers=min(workers, len(namespace_order)), mp_context=process_pool_context()) as executor:
        futures = {
            executor.submit(_size_filter_namespace_task, namespace, mapped_pre_filter[namespace], cfg): namespace
            for namespace in namespace_order
        }
        for future in as_completed(futures):
            namespace, size_kept, removed_size = future.result()
            size_results[namespace] = (size_kept, removed_size)

    size_filtered = {namespace: size_results[namespace][0] for namespace in namespace_order}
    removed_size_frames = [size_results[namespace][1] for namespace in namespace_order]
    redundancy_tasks: list[tuple[str, str, dict[str, dict[str, Any]], dict]] = []
    for namespace in namespace_order:
        families = sorted({record["meta"].get("family", "") for record in size_filtered[namespace].values()})
        for family in families:
            subset = {
                term_id: record
                for term_id, record in size_filtered[namespace].items()
                if record["meta"].get("family", "") == family
            }
            redundancy_tasks.append((namespace, family, subset, cfg))

    filtered = {namespace: {} for namespace in namespace_order}
    if not redundancy_tasks:
        return size_filtered, filtered, _concat_any(removed_size_frames), pd.DataFrame(), pd.DataFrame()

    LOGGER.info("Removing redundancy with %d workers across %d namespace/family groups", min(workers, len(redundancy_tasks)), len(redundancy_tasks))
    redundancy_results: list[tuple[str, str, dict[str, dict[str, Any]], pd.DataFrame, pd.DataFrame]] = []
    with ProcessPoolExecutor(max_workers=min(workers, len(redundancy_tasks)), mp_context=process_pool_context()) as executor:
        futures = [executor.submit(_redundancy_family_task, *task) for task in redundancy_tasks]
        for future in as_completed(futures):
            redundancy_results.append(future.result())

    namespace_order_index = {namespace: i for i, namespace in enumerate(namespace_order)}
    removed_redundancy_frames: list[pd.DataFrame] = []
    redundancy_summary_frames: list[pd.DataFrame] = []
    for namespace, _, final_terms, removed_redundancy, redundancy_summary in sorted(
        redundancy_results,
        key=lambda item: (namespace_order_index[item[0]], item[1]),
    ):
        filtered[namespace].update(final_terms)
        removed_redundancy_frames.append(removed_redundancy)
        redundancy_summary_frames.append(redundancy_summary)

    return (
        size_filtered,
        filtered,
        _concat_any(removed_size_frames),
        _concat_any(removed_redundancy_frames),
        _concat_any(redundancy_summary_frames),
    )


def _size_filter_namespace_task(
    namespace: str,
    terms: dict[str, dict[str, Any]],
    cfg: dict,
) -> tuple[str, dict[str, dict[str, Any]], pd.DataFrame]:
    size_kept, removed_size = apply_size_and_broad_filters(terms, namespace, cfg)
    return namespace, size_kept, removed_size


def _redundancy_family_task(
    namespace: str,
    family: str,
    terms: dict[str, dict[str, Any]],
    cfg: dict,
) -> tuple[str, str, dict[str, dict[str, Any]], pd.DataFrame, pd.DataFrame]:
    final_terms, removed_redundancy, redundancy_summary = remove_redundancy(terms, namespace, cfg)
    return namespace, family, final_terms, removed_redundancy, redundancy_summary


def _target_gene_restrictions(targets: list[dict[str, Any]]) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for target in targets:
        namespace = str(target["name"])
        source_path = target.get("source_path")
        source_column = target.get("source_column")
        universe_table = target.get("universe_table", pd.DataFrame())
        universe = target.get("universe")
        annotation = target["annotation"]
        annotation_ids = {strip_ensembl_version(gene) for gene in annotation.get("ensembl_gene_id", pd.Series(dtype=str))}
        if universe_table.empty and universe:
            universe_table = pd.DataFrame({"input_gene_id": sorted(universe), "ensembl_gene_id": sorted(universe), "id_type": "ensembl_stable"})
        for _, row in universe_table.iterrows():
            stable = strip_ensembl_version(row.get("ensembl_gene_id", ""))
            if not stable:
                continue
            rows.append(
                {
                    "target_namespace": namespace,
                    "source_path": str(source_path or ""),
                    "source_column": str(source_column or ""),
                    "input_gene_id": row.get("input_gene_id", ""),
                    "ensembl_gene_id": stable,
                    "id_type": row.get("id_type", ""),
                    "has_annotation_metadata": stable in annotation_ids,
                }
            )
    return pd.DataFrame(rows).drop_duplicates()


def _target_universe_summary(targets: list[dict[str, Any]], source_terms: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    source_ids = {
        strip_ensembl_version(value)
        for value in source_terms.get("gene_ensembl_from_source", pd.Series(dtype=str))
        if strip_ensembl_version(value)
    }
    for target in targets:
        namespace = str(target["name"])
        universe = target.get("universe")
        annotation = target["annotation"]
        universe_ids = {strip_ensembl_version(gene) for gene in universe} if universe else {
            strip_ensembl_version(gene) for gene in annotation.get("ensembl_gene_id", pd.Series(dtype=str))
        }
        annotation_ids = {strip_ensembl_version(gene) for gene in annotation.get("ensembl_gene_id", pd.Series(dtype=str))}
        annotated_ids = universe_ids & annotation_ids
        missing_annotation = universe_ids - annotated_ids
        rows.append(
            {
                "target_namespace": namespace,
                "matrix_universe_size": len(universe_ids),
                "annotation_helper_size": len(annotation_ids),
                "ids_in_both": len(annotated_ids),
                "matrix_ids_missing_annotation_metadata": len(missing_annotation),
                "annotation_ids_absent_from_matrix": len(annotation_ids - universe_ids),
                "source_ensembl_ids_kept_without_annotation": len(source_ids & missing_annotation),
            }
        )
    return pd.DataFrame(rows)


def _concat_any(frames: list[pd.DataFrame]) -> pd.DataFrame:
    frames = [frame for frame in frames if frame is not None and not frame.empty]
    return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()


def _with_columns(df: pd.DataFrame, columns: list[str]) -> pd.DataFrame:
    if df.empty:
        return pd.DataFrame(columns=columns)
    for col in columns:
        if col not in df:
            df[col] = ""
    return df[columns]


def _report_ensembl_archive_audit(audit: pd.DataFrame, cfg: dict) -> None:
    if audit.empty:
        LOGGER.info("Ensembl archive audit: no source Ensembl IDs needed archive fallback.")
        return
    summary = audit.groupby(["target_namespace", "rescue_status"]).size().reset_index(name="n")
    LOGGER.info("Ensembl archive fallback/audit summary:\n%s", summary.to_string(index=False))
    rescued_statuses = {
        "archive_ensembl_rescue",
        "hgnc_current_symbol",
        "hgnc_previous_symbol",
        "hgnc_alias_symbol",
        "hgnc_ensembl_tiebreak",
        "manual_symbol_repair",
        "ncbi_gene_history",
        "ncbi_gene_info_synonym",
    }
    unresolved = audit[~audit["rescue_status"].isin(rescued_statuses)]
    if unresolved.empty:
        LOGGER.info("Ensembl archive audit: all %d problematic source Ensembl mappings were resolved by archive or symbol fallback.", len(audit))
        return
    limit = int(cfg.get("ensembl_archive", {}).get("unresolved_print_limit", 20))
    LOGGER.warning("Ensembl archive audit: %d problematic source Ensembl mappings were not resolved.", len(unresolved))
    for _, row in unresolved.head(limit).iterrows():
        LOGGER.warning(
            "Unrescued Ensembl ID example: namespace=%s term=%s source=%s symbol=%s ensembl=%s status=%s",
            row.get("target_namespace", ""),
            row.get("term_id", ""),
            row.get("source_tag", ""),
            row.get("input_symbol", ""),
            row.get("source_ensembl_stable", ""),
            row.get("rescue_status", ""),
        )


def _ensembl_archive_summary(audit: pd.DataFrame) -> pd.DataFrame:
    if audit.empty:
        return pd.DataFrame(columns=["target_namespace", "rescue_status", "n"])
    return audit.groupby(["target_namespace", "rescue_status"]).size().reset_index(name="n")


def _non_gene_summary(non_gene: pd.DataFrame) -> pd.DataFrame:
    columns = ["target_namespace", "family", "source_tag", "reason", "n"]
    if non_gene.empty:
        return pd.DataFrame(columns=columns)
    return non_gene.groupby(columns[:-1]).size().reset_index(name="n")


def _ncbi_gene_summary(audit: pd.DataFrame) -> pd.DataFrame:
    columns = ["target_namespace", "ncbi_stage", "ncbi_status", "action", "n"]
    if audit.empty:
        return pd.DataFrame(columns=columns)
    return audit.groupby(columns[:-1]).size().reset_index(name="n")
