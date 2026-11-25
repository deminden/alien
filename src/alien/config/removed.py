"""Options that existed in 0.1.x but never had an effect."""

from __future__ import annotations

#: Options that existed in 0.1.x but never had an effect; using one is an error naming the replacement.
REMOVED_OPTIONS: dict[str, dict[str, str]] = {
    "redundancy": {
        "apply_after_output_gene_projection": "redundancy always runs after target projection",
        "nested_overlap_filter": "the nested-overlap filter was never implemented; remove the block",
    },
    "filtering": {
        "remove_empty_after_output_gene_projection": "empty terms are always dropped; use min_size_default",
    },
    "gene_mapping": {
        "strip_ensembl_version": "Ensembl versions are always stripped",
        "allow_ambiguous_alias_mapping": "ambiguous aliases are always dropped and audited",
        "prefer_gencode_exact_gene_name": "annotation gene names are always preferred",
        "prefer_current_hgnc_symbol": "current HGNC symbols are always preferred",
        "drop_unmapped_from_ensg_gmt": "unmapped genes are never written to Ensembl GMTs",
    },
    "ensembl_archive": {"unresolved_example_threshold": "use ensembl_archive.unresolved_print_limit"},
}
