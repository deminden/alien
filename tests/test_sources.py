import hashlib
import gzip
import zipfile

import pandas as pd
import rdata

from alien.config import load_config
from alien import sources
from alien.sources import load_ensembl_gtf_target, read_canonical_memberships, read_msigdb_like_table, read_msigdb_remote, read_symbol_gmt_source


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


def test_gencode_target_version_downloads_without_explicit_path(tmp_path, monkeypatch):
    def fake_download(url, path, force=False):
        assert url == sources.GENCODE_URLS["47"]
        path.parent.mkdir(parents=True, exist_ok=True)
        with gzip.open(path, "wt", encoding="utf-8") as handle:
            handle.write(
                'chr1\tALIEN\tgene\t1\t10\t.\t+\t.\tgene_id "ENSG00000141510.18"; gene_name "TP53"; gene_type "protein_coding";\n'
            )

    monkeypatch.setattr(sources, "download_file", fake_download)

    annotation, label, path = load_ensembl_gtf_target(
        tmp_path,
        {
            "name": "human_gencode47",
            "type": "ensembl_gtf",
            "annotation": {"source": "GENCODE", "version": "47"},
        },
    )

    assert label == "GENCODE v47"
    assert path == tmp_path / "gencode" / "gencode.v47.annotation.gtf.gz"
    assert annotation["ensembl_gene_id"].tolist() == ["ENSG00000141510"]


def test_ensembl_gtf_target_allows_custom_attribute_names(tmp_path):
    gtf = tmp_path / "custom.gtf"
    gtf.write_text(
        'chr1\tALIEN\tgene\t1\t10\t.\t+\t.\tID "ENSG00000141510.18"; Name "TP53"; biotype "protein_coding";\n'
        'chr1\tALIEN\tgene\t20\t30\t.\t+\t.\tID=ENSG000002.1;Name=GENE2;biotype=lncRNA;\n',
        encoding="utf-8",
    )

    annotation, _, _ = load_ensembl_gtf_target(
        tmp_path,
        {
            "name": "custom_annotation",
            "type": "ensembl_gtf",
            "annotation": {
                "source": "Custom",
                "path": str(gtf),
                "attributes": {
                    "gene_id": "ID",
                    "gene_name": "Name",
                    "gene_biotype": "biotype",
                },
            },
        },
    )

    assert annotation["ensembl_gene_id"].tolist() == ["ENSG00000141510", "ENSG000002"]
    assert annotation["gene_symbol"].tolist() == ["TP53", "GENE2"]
    assert annotation["gene_biotype"].tolist() == ["protein_coding", "lncRNA"]


def test_msigdb_remote_reads_cached_rds_release(tmp_path):
    cache_dir = tmp_path / "remote"
    release_dir = tmp_path / "release"
    release_dir.mkdir()
    summary_path = release_dir / "msigdb.test.summary.rds"
    c2_path = release_dir / "c2.rds"
    rdata.write_rds(
        summary_path,
        pd.DataFrame(
            {
                "db_target_species": ["HS"],
                "gs_collection": ["C2"],
                "gs_subcollection": ["CP:REACTOME"],
                "df_rds": ["c2.rds"],
            }
        ).astype(object),
    )
    rdata.write_rds(
        c2_path,
        pd.DataFrame(
            {
                "db_version": ["test", "test"],
                "db_target_species": ["HS", "HS"],
                "db_gene_symbol": ["TP53", "GENE2"],
                "db_ensembl_gene": ["ENSG00000141510", "ENSG000002"],
                "gs_id": ["M1", "M1"],
                "gs_name": ["REACTOME_SIGNAL", "REACTOME_SIGNAL"],
                "gs_description": ["desc", "desc"],
                "gs_collection": ["C2", "C2"],
                "gs_subcollection": ["CP:REACTOME", "CP:REACTOME"],
                "gs_url": ["https://example.test", "https://example.test"],
            }
        ).astype(object),
    )
    zip_path = cache_dir / "msigdb.test.zip"
    cache_dir.mkdir()
    with zipfile.ZipFile(zip_path, "w") as archive:
        archive.write(summary_path, arcname=summary_path.name)
        archive.write(c2_path, arcname=c2_path.name)
    digest = hashlib.md5(zip_path.read_bytes()).hexdigest()

    frame = read_msigdb_remote(
        {
            "cache_dir": str(cache_dir),
            "db_species": "HS",
            "collection": "C2",
            "release": {
                "zip_url": "https://example.test/msigdb.test.zip",
                "zip_md5": digest,
                "zip_name": "msigdb.test.zip",
                "summary_rds": "msigdb.test.summary.rds",
            },
        },
        tmp_path,
    )

    assert frame["term_id"].unique().tolist() == ["MSIGDB_REACTOME__REACTOME_SIGNAL"]
    assert frame["source_tag"].unique().tolist() == ["REACTOME"]
    assert frame["gene_symbol"].tolist() == ["TP53", "GENE2"]
    assert frame["gene_ensembl_from_source"].tolist() == ["ENSG00000141510", "ENSG000002"]


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
