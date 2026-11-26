"""Every table ALIEN reads or writes, declared once.

``docs/outputs.md`` is generated from :data:`OUTPUT_FILES` (``python scripts/generate_outputs_doc.py``),
so the documentation cannot drift from the code.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Final

#: Bumped whenever an output file or column changes incompatibly.
OUTPUT_SCHEMA_VERSION: Final = 3

#: Bumped whenever a derived cache (normalized MSigDB parquet, NCBI maps, parsed GTF) changes layout.
CACHE_SCHEMA_VERSION: Final = 1

CANONICAL_COLUMNS: Final = (
    "term_id",
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
    "gene_symbol",
    "gene_ensembl_from_source",
    "metadata_json",
)
#: Columns that identify a term; two sources producing the same term_id with different values collide.
TERM_ID_IDENTITY_COLUMNS: Final = tuple(
    column for column in CANONICAL_COLUMNS if column not in {"term_id", "gene_symbol", "gene_ensembl_from_source"}
)
TERM_META_COLUMNS: Final = tuple(column for column in CANONICAL_COLUMNS if column not in {"gene_symbol"})

ROW_IDENTITY: Final = ("target_namespace", "family", "term_id", "source_tag")
UNMAPPED_COLUMNS: Final = (*ROW_IDENTITY, "input_gene", "attempted_current_symbol", "reason", "possible_matches")
NON_GENE_COLUMNS: Final = (*ROW_IDENTITY, "input_gene", "reason")
AMBIGUOUS_COLUMNS: Final = (
    *ROW_IDENTITY,
    "input_gene",
    "mapping_stage",
    "possible_symbols",
    "possible_ensembl_ids",
    "action",
)
ARCHIVE_AUDIT_COLUMNS: Final = (
    *ROW_IDENTITY,
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
)
NCBI_AUDIT_COLUMNS: Final = (
    *ROW_IDENTITY,
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
)
MAPPING_STATUS_COLUMNS: Final = ("target_namespace", "family", "source_tag", "mapping_status", "n_gene_members")
SIZE_REMOVED_COLUMNS: Final = (
    "target_namespace",
    "family",
    "term_id",
    "original_name",
    "source_tag",
    "original_size_symbols",
    "mapped_size",
    "reason",
    "exclusion_rule",
)
SYMBOL_RESCUE_COLUMNS: Final = (
    "target_namespace",
    "source_tag",
    "input_gene",
    "source_ensembl_id",
    "mapping_status",
    "mapped_symbol",
    "ensembl_gene_id",
    "n_memberships",
    "n_terms",
    "example_term_id",
)
REDUNDANCY_REMOVED_COLUMNS: Final = (
    "target_namespace",
    "family",
    "term_id",
    "representative_term_id",
    "redundancy_cluster_id",
    "reason",
    "mapped_size",
)
REDUNDANCY_SUMMARY_COLUMNS: Final = (
    "target_namespace",
    "family",
    "n_input",
    "n_after_exact",
    "n_after_jaccard",
    "n_removed_exact",
    "n_removed_jaccard",
    "jaccard_cutoff",
)
TERM_MANIFEST_COLUMNS: Final = (
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
)
SOURCE_MANIFEST_COLUMNS: Final = (
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
)
TERM_ID_COLLISION_COLUMNS: Final = (
    "term_id",
    "n_identities_for_term",
    "collision_identity_id",
    "n_rows",
    "n_gene_symbols",
    "sample_gene_symbols",
    *TERM_ID_IDENTITY_COLUMNS,
    "action",
    "resolved_term_id",
)
ANNOTATION_COLUMNS: Final = (
    "gencode_version",
    "ensembl_gene_id_versioned",
    "ensembl_gene_id",
    "gene_symbol",
    "gene_biotype",
    "seqname",
    "start",
    "end",
    "strand",
    "source_gtf",
    "target_id_metadata_source",
    "is_in_output_genes",
)
GENE_ID_LIST_COLUMNS: Final = ("input_gene_id", "ensembl_gene_id", "gene_symbol", "id_type")
TARGET_GENE_ID_COLUMNS: Final = (
    "target_namespace",
    "source_path",
    "source_id_column",
    "source_symbol_column",
    "input_gene_id",
    "ensembl_gene_id",
    "gene_symbol",
    "id_type",
    "has_annotation_metadata",
)
METADATA_FALLBACK_COLUMNS: Final = (
    "target_namespace",
    "fallback_index",
    "mode",
    "annotation",
    "annotation_path",
    "n_fallback_genes",
    "n_missing_before",
    "n_rows_added",
    "n_missing_after",
)
TARGET_NAMESPACE_SUMMARY_COLUMNS: Final = (
    "target_namespace",
    "output_gene_source",
    "output_genes_size",
    "gene_filter_size",
    "effective_target_size",
    "annotation_helper_size",
    "primary_annotation_helper_size",
    "metadata_fallback_helper_size",
    "output_genes_with_annotation_metadata",
    "annotation_metadata_coverage",
    "target_ids_missing_annotation_metadata",
    "annotation_ids_absent_from_target",
    "source_ensembl_ids_kept_without_annotation",
)
COLLECTION_SUMMARY_COLUMNS: Final = (
    "target_namespace",
    "family",
    "source_tag",
    "n_terms_before_filter",
    "n_terms_after_size_filter",
    "n_terms_after_redundancy",
    "median_size",
    "min_size",
    "max_size",
    "n_unique_genes",
)
MAPPING_SUMMARY_COLUMNS: Final = (
    "target_namespace",
    "family",
    "source_tag",
    "n_gene_members_original",
    "n_non_gene_token_dropped",
    "n_gene_members_considered_for_mapping",
    "n_gene_members_mapped_raw",
    "n_gene_members_mapped_unique",
    "n_gene_members_collapsed_after_mapping",
    "n_gene_members_unmapped",
    "n_ambiguous",
    "mapping_rate_raw",
    "mapping_rate_unique",
    "n_terms_lost_because_too_small",
    "n_terms_lost_because_too_large",
)
ARCHIVE_SUMMARY_COLUMNS: Final = ("target_namespace", "rescue_status", "n")
NON_GENE_SUMMARY_COLUMNS: Final = ("target_namespace", "family", "source_tag", "reason", "n")
NCBI_SUMMARY_COLUMNS: Final = ("target_namespace", "ncbi_stage", "ncbi_status", "action", "n")


@dataclass(frozen=True)
class OutputFile:
    path: str
    description: str
    columns: tuple[str, ...] = ()
    modes: tuple[str, ...] = ("full",)


OUTPUT_FILES: Final = (
    OutputFile(
        "gmt/<namespace>.gmt",
        "One GMT per namespace (symbols plus each target). Terms sorted by ID; genes sorted; empty terms omitted.",
        modes=("full", "minimal", "gmt"),
    ),
    OutputFile(
        "MANIFEST.json",
        "Every file in the build with sha256 and size, the ALIEN version, output schema version and mode. "
        "Written last; `alien verify` re-hashes against it.",
        modes=("full", "minimal", "gmt"),
    ),
    OutputFile(
        "metadata/effective_config.yml",
        "The validated configuration used, without project.outdir.",
        modes=("full", "minimal"),
    ),
    OutputFile(
        "metadata/source_provenance.json",
        "ALIEN version, packages, platform, config sha256, every input file and remote resource with hashes, "
        "redundancy settings, worker count and warnings. The only file with a timestamp.",
        modes=("full", "minimal"),
    ),
    OutputFile(
        "metadata/source_manifest.tsv",
        "One row per loaded source collection.",
        SOURCE_MANIFEST_COLUMNS,
        ("full", "minimal"),
    ),
    OutputFile(
        "metadata/term_manifest.tsv.gz",
        "One row per (term, namespace) with the namespace-specific outcome.",
        TERM_MANIFEST_COLUMNS,
        ("full", "minimal"),
    ),
    OutputFile(
        "metadata/target_metadata_fallbacks.tsv",
        "What each annotation metadata fallback added.",
        METADATA_FALLBACK_COLUMNS,
        ("full", "minimal"),
    ),
    OutputFile(
        "metadata/target_output_genes.tsv.gz",
        "Configured output_genes per target.",
        TARGET_GENE_ID_COLUMNS,
        ("full", "minimal"),
    ),
    OutputFile(
        "metadata/target_gene_filter.tsv.gz",
        "Configured gene_filter IDs per target.",
        TARGET_GENE_ID_COLUMNS,
        ("full", "minimal"),
    ),
    OutputFile(
        "metadata/gene_mapping_<target>.tsv.gz",
        "The annotation helper used for one target (GTF, fallbacks and output-gene symbol rows).",
        ANNOTATION_COLUMNS,
    ),
    OutputFile("metadata/gene_symbol_master.tsv.gz", "The HGNC complete set used for symbol resolution."),
    OutputFile(
        "metadata/removed_terms_size_filter.tsv.gz",
        "Terms removed by size or broad-term filters.",
        SIZE_REMOVED_COLUMNS,
    ),
    OutputFile(
        "metadata/removed_terms_redundancy.tsv.gz",
        "Terms removed as exact or Jaccard duplicates.",
        REDUNDANCY_REMOVED_COLUMNS,
    ),
    OutputFile("metadata/unmapped_genes.tsv.gz", "Memberships that could not be mapped.", UNMAPPED_COLUMNS),
    OutputFile(
        "metadata/symbol_rescues.tsv.gz",
        "Each distinct source token mapped by anything other than a source Ensembl ID or an exact current symbol "
        "(previous symbol, alias, case, legacy name, NCBI, archive), with how many memberships and terms it reached.",
        SYMBOL_RESCUE_COLUMNS,
    ),
    OutputFile(
        "metadata/dropped_non_gene_tokens.tsv.gz",
        "Memberships dropped as configured non-gene tokens.",
        NON_GENE_COLUMNS,
    ),
    OutputFile(
        "metadata/ambiguous_gene_mappings.tsv.gz",
        "Memberships dropped because a symbol matched several genes.",
        AMBIGUOUS_COLUMNS,
    ),
    OutputFile(
        "metadata/ensembl_archive_audit.tsv.gz",
        "Every source Ensembl ID checked against the Ensembl archive and its outcome.",
        ARCHIVE_AUDIT_COLUMNS,
    ),
    OutputFile(
        "metadata/ncbi_gene_rescue_audit.tsv.gz", "Every NCBI Gene rescue attempt and its outcome.", NCBI_AUDIT_COLUMNS
    ),
    OutputFile(
        "metadata/term_id_collisions.tsv",
        "Only when two sources produce the same term_id with different identities.",
        TERM_ID_COLLISION_COLUMNS,
    ),
    OutputFile(
        "qc/collection_summary.tsv", "Term counts and sizes after each filter stage.", COLLECTION_SUMMARY_COLUMNS
    ),
    OutputFile("qc/collection_summary.json", "The same table as JSON records."),
    OutputFile(
        "qc/mapping_summary.tsv",
        "Mapping counts and rates per namespace, family and source.",
        MAPPING_SUMMARY_COLUMNS,
        ("full", "minimal"),
    ),
    OutputFile("qc/mapping_status_summary.tsv", "Membership counts per mapping status.", MAPPING_STATUS_COLUMNS),
    OutputFile(
        "qc/redundancy_summary.tsv",
        "Term counts through exact and Jaccard redundancy removal.",
        REDUNDANCY_SUMMARY_COLUMNS,
        ("full", "minimal"),
    ),
    OutputFile(
        "qc/target_namespace_summary.tsv",
        "Target ID-set sizes and annotation coverage.",
        TARGET_NAMESPACE_SUMMARY_COLUMNS,
        ("full", "minimal"),
    ),
    OutputFile("qc/ensembl_archive_summary.tsv", "Archive audit counts by status.", ARCHIVE_SUMMARY_COLUMNS),
    OutputFile("qc/non_gene_token_summary.tsv", "Dropped non-gene tokens by source.", NON_GENE_SUMMARY_COLUMNS),
    OutputFile("qc/ncbi_gene_rescue_summary.tsv", "NCBI rescue counts by stage and outcome.", NCBI_SUMMARY_COLUMNS),
    OutputFile("qc/warnings.txt", "Build warnings, one per line.", modes=("full", "minimal")),
)
