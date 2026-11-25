from __future__ import annotations

from copy import deepcopy
from typing import Any

DEFAULT_CONFIG: dict[str, Any] = {
    "project": {
        "name": "alien",
        "source_dir": "data/alien_sources",
        "outdir": "data/alien_gmt",
    },
    "outputs": {
        "include_symbols": True,
    },
    "term_id_collisions": {
        "action": "error",
    },
    "sources": [
        {
            "type": "msigdb_remote",
            "version": "2026.1",
            "db_species": "HS",
            "collection": "C2",
        },
        {
            "type": "symbol_gmt",
            "name": "local_symbolic_libraries",
            "paths": [],
            "source_tag": "LOCAL",
            "family": "biology_process_pathway",
            "aspect": "pathway",
        },
    ],
    "targets": [
        {
            "name": "human_gencode49",
            "type": "ensembl_gtf",
            "annotation": {
                "source": "GENCODE",
                "version": "49",
            },
        }
    ],
    "filtering": {
        "min_size_default": 10,
        "max_size_default": 500,
        "max_size_disease": 500,
        "max_size_cancer": 500,
        "drop_broad_disease_terms": True,
        "broad_disease_regex": [
            "^disease$",
            "^diseases$",
            "^cancer$",
            "^neoplasm$",
            "^neoplasms$",
            "^malignant neoplasm$",
            "^malignant neoplasms$",
            "^carcinoma$",
            "^tumor$",
            "^tumour$",
            "^genetic disease$",
            "^genetic disorder$",
        ],
    },
    "redundancy": {
        "exact_duplicate_removal": True,
        "jaccard_cutoff": 0.85,
        "apply_per_family": True,
    },
    "source_priority": {
        "biology_process_pathway": ["REACTOME", "WIKIPATHWAYS", "KEGG_MEDICUS", "GOBP"],
        "biology_function_location": ["GOMF", "GOCC"],
        "disease_phenotype": [
            "HPO",
            "DISGENET",
            "CLINVAR",
            "GWAS_CATALOG",
            "JENSEN_DISEASES_CURATED",
            "JENSEN_DISEASES_EXPERIMENTAL",
            "JENSEN_DISEASES_TEXTMINING",
        ],
        "cancer_dependency_state": [
            "C9",
            "C6",
            "C4_3CA",
            "C4_CGN",
            "C4_CM",
            "DEPMAP",
            "CCLE",
            "NCI60",
        ],
    },
    "gene_mapping": {
        "use_hgnc_alias_mapping": True,
        "preserve_unmapped_symbols_in_symbol_gmt": True,
        "non_gene_tokens": [
            "disease",
            "diseases",
            "phenotype",
            "phenotypes",
            "group",
            "groups",
        ],
        "non_gene_token_regex": [],
        "manual_symbol_repairs": {},
    },
    "ncbi_gene": {
        "enabled": True,
        "allowed_source_tags": ["CCLE", "NCI60", "DEPMAP"],
        "gene_info_url": "https://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/Homo_sapiens.gene_info.gz",
        "gene_history_url": "https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene_history.gz",
        "use_gene_info_synonyms": True,
    },
    "ensembl_archive": {
        "enabled": True,
        "endpoint": "https://rest.ensembl.org/archive/id",
        "batch_size": 50,
        "connect_timeout_seconds": 30,
        "timeout_seconds": 120,
        "retry_attempts": 3,
        "retry_backoff_seconds": 5,
        "unresolved_print_limit": 20,
    },
}


def default_config() -> dict[str, Any]:
    return deepcopy(DEFAULT_CONFIG)
