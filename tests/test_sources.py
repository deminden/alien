import pandas as pd

from alien.config import load_config
from alien.sources import read_canonical_memberships, read_msigdb_like_table, read_symbol_gmt_source


def test_canonical_table_normalization_with_source_gene_id(tmp_path):
    path = tmp_path / "canonical.tsv"
    pd.DataFrame(
        {
            "term_id": ["Term 1", "Term 1"],
            "gene_symbol": ["TP53", "GENE2"],
            "gene_id": ["ENSG00000141510.18", "ENSG000002.1"],
            "gene_id_namespace": ["ensembl_gene", "ensembl_gene"],
        }
    ).to_csv(path, sep="\t", index=False)

    frame = read_canonical_memberships(path, {"source": "Example", "source_tag": "EX", "family": "biology_process_pathway"})

    assert frame["term_id"].tolist() == ["Term_1", "Term_1"]
    assert frame["source_tag"].unique().tolist() == ["EX"]
    assert frame["gene_ensembl_from_source"].tolist() == ["ENSG00000141510.18", "ENSG000002.1"]


def test_msigdb_like_table_normalization(tmp_path):
    path = tmp_path / "msigdb.tsv"
    pd.DataFrame(
        {
            "gs_name": ["PATHWAY_A", "PATHWAY_A"],
            "gs_description": ["desc", "desc"],
            "gene_symbol": ["TP53", "GENE2"],
            "ensembl_gene": ["ENSG00000141510.18", "ENSG000002.1"],
        }
    ).to_csv(path, sep="\t", index=False)

    frame = read_msigdb_like_table(path, {"source_tag": "REACTOME", "collection": "C2", "family": "biology_process_pathway", "aspect": "pathway"})

    assert frame["term_id"].unique().tolist() == ["MSIGDB_REACTOME__PATHWAY_A"]
    assert frame["gene_symbol"].tolist() == ["TP53", "GENE2"]
    assert frame["gene_ensembl_from_source"].tolist() == ["ENSG00000141510.18", "ENSG000002.1"]


def test_symbol_gmt_source_normalization(tmp_path):
    path = tmp_path / "library.gmt"
    path.write_text("TERM_A\tdesc\tTP53\tGENE2\n", encoding="utf-8")

    frame = read_symbol_gmt_source(path, {"source": "Example", "source_tag": "EX", "family": "disease_phenotype", "aspect": "disease"})

    assert frame["term_id"].unique().tolist() == ["EX__TERM_A"]
    assert frame["gene_symbol"].tolist() == ["TP53", "GENE2"]
    assert frame["gene_ensembl_from_source"].tolist() == ["", ""]


def test_default_filter_lists_can_be_extended(tmp_path):
    config = tmp_path / "config.yml"
    config.write_text(
        """
gene_mapping:
  extra_non_gene_tokens:
    - signature
filtering:
  extra_broad_disease_regex:
    - "^syndrome$"
sources:
  - type: symbol_gmt
    path: example.gmt
targets: []
""",
        encoding="utf-8",
    )

    cfg = load_config(config)

    assert "signature" in cfg["gene_mapping"]["non_gene_tokens"]
    assert "^syndrome$" in cfg["filtering"]["broad_disease_regex"]
