import pandas as pd
import pytest

from alien import build
from alien.build import _term_id_collision_audit


def test_small_build_is_deterministic_and_writes_combined_gmts(tmp_path):
    source = tmp_path / "sources"
    out1 = tmp_path / "out1"
    out2 = tmp_path / "out2"
    source_tsv = tmp_path / "source.tsv"
    symbol_gmt = tmp_path / "symbolic.gmt"
    gtf = tmp_path / "genes.gtf"
    output_genes = tmp_path / "output_genes.tsv"
    _write_hgnc(source)
    _write_source_table(source_tsv)
    symbol_gmt.write_text("SYMBOL_TERM\tdesc\tTP53\tGENE3\n", encoding="utf-8")
    _write_gtf(gtf)
    output_genes.write_text("feature\nENSG00000141510.18\nENSG000002.1\nENSG000003.1\n", encoding="utf-8")

    cfg = {
        "project": {"source_dir": str(source)},
        "sources": [
            {
                "type": "msigdb_tsv",
                "path": str(source_tsv),
                "source_tag": "REACTOME",
                "collection": "C2",
                "family": "biology_process_pathway",
                "aspect": "pathway",
            },
            {
                "type": "symbol_gmt",
                "path": str(symbol_gmt),
                "source": "Example symbolic",
                "source_tag": "LOCAL",
                "family": "biology_process_pathway",
                "aspect": "pathway",
            },
        ],
        "targets": [
            {
                "name": "human_test",
                "type": "ensembl_gtf",
                "annotation": {
                    "source": "TestAnnotation",
                    "version": "test",
                    "path": str(gtf),
                },
                "output_genes": {"path": str(output_genes), "id_column": "feature"},
            }
        ],
        "filtering": {
            "min_size_default": 1,
            "max_size_default": 10,
            "max_size_disease": 10,
            "max_size_cancer": 10,
        },
        "source_priority": {"biology_process_pathway": ["REACTOME", "LOCAL"]},
        "ncbi_gene": {"enabled": False},
        "ensembl_archive": {"enabled": False},
    }

    result1 = build(config=cfg, outdir=out1, workers=1)
    result2 = build(config=cfg, outdir=out2, workers=2)

    assert result1.namespaces == ("symbols", "human_test")
    assert (out1 / "gmt" / "symbols.gmt").exists()
    assert (out1 / "gmt" / "human_test.gmt").exists()
    assert (out1 / "metadata" / "source_manifest.tsv").exists()
    assert (out1 / "metadata" / "term_manifest.tsv.gz").exists()
    assert (out1 / "metadata" / "target_output_genes.tsv.gz").exists()
    assert (out1 / "qc" / "warnings.txt").exists()
    source_manifest = pd.read_csv(out1 / "metadata" / "source_manifest.tsv", sep="\t", dtype=str).fillna("")
    assert set(source_manifest["source_tag"]) == {"LOCAL", "REACTOME"}
    assert source_manifest.set_index("source_tag").loc["REACTOME", "n_terms"] == "1"
    mapping = pd.read_csv(out1 / "metadata" / "gene_mapping_human_test.tsv.gz", sep="\t", dtype=str)
    assert mapping["source_gtf"].unique().tolist() == [str(gtf)]
    assert (out1 / "gmt" / "symbols.gmt").read_text(encoding="utf-8") == (out2 / "gmt" / "symbols.gmt").read_text(encoding="utf-8")
    assert (out1 / "gmt" / "human_test.gmt").read_text(encoding="utf-8") == (out2 / "gmt" / "human_test.gmt").read_text(encoding="utf-8")


def test_term_id_collision_audit_ignores_membership_rows_for_same_term_identity():
    source_terms = pd.DataFrame(
        {
            "term_id": ["TERM", "TERM"],
            "source_tag": ["LOCAL", "LOCAL"],
            "collection": ["example", "example"],
            "original_name": ["Term", "Term"],
            "gene_symbol": ["TP53", "GENE2"],
        }
    )

    assert _term_id_collision_audit(source_terms).empty


def test_restrict_to_key_is_not_supported(tmp_path):
    source = tmp_path / "sources"
    source_tsv = tmp_path / "source.tsv"
    gtf = tmp_path / "genes.gtf"
    _write_hgnc(source)
    _write_source_table(source_tsv)
    _write_gtf(gtf)

    cfg = {
        "project": {"source_dir": str(source)},
        "outputs": {"include_symbols": False},
        "sources": [
            {
                "type": "msigdb_tsv",
                "path": str(source_tsv),
                "source_tag": "REACTOME",
                "collection": "C2",
                "family": "biology_process_pathway",
                "aspect": "pathway",
            }
        ],
        "targets": [
            {
                "name": "human_test",
                "type": "ensembl_gtf",
                "annotation": {"source": "TestAnnotation", "version": "test", "path": str(gtf)},
                "restrict_to": {"ids": ["ENSG00000141510"]},
            }
        ],
        "ncbi_gene": {"enabled": False},
        "ensembl_archive": {"enabled": False},
    }

    with pytest.raises(ValueError, match="use 'output_genes' or 'gene_filter' instead"):
        build(config=cfg, outdir=tmp_path / "out", workers=1)


def test_gene_universe_key_is_not_supported(tmp_path):
    source = tmp_path / "sources"
    source_tsv = tmp_path / "source.tsv"
    gtf = tmp_path / "genes.gtf"
    _write_hgnc(source)
    _write_source_table(source_tsv)
    _write_gtf(gtf)

    cfg = {
        "project": {"source_dir": str(source)},
        "outputs": {"include_symbols": False},
        "sources": [{"type": "msigdb_tsv", "path": str(source_tsv), "source_tag": "REACTOME"}],
        "targets": [
            {
                "name": "bad_target",
                "type": "ensembl_gtf",
                "annotation": {"source": "TestAnnotation", "version": "test", "path": str(gtf)},
                "gene_universe": {"ids": ["ENSG00000141510"]},
            }
        ],
        "ncbi_gene": {"enabled": False},
        "ensembl_archive": {"enabled": False},
    }

    with pytest.raises(ValueError, match="use 'output_genes' instead"):
        build(config=cfg, outdir=tmp_path / "out", workers=1)


def test_output_ids_key_is_not_supported(tmp_path):
    source = tmp_path / "sources"
    source_tsv = tmp_path / "source.tsv"
    gtf = tmp_path / "genes.gtf"
    _write_hgnc(source)
    _write_source_table(source_tsv)
    _write_gtf(gtf)

    cfg = {
        "project": {"source_dir": str(source)},
        "outputs": {"include_symbols": False},
        "sources": [{"type": "msigdb_tsv", "path": str(source_tsv), "source_tag": "REACTOME"}],
        "targets": [
            {
                "name": "bad_target",
                "type": "ensembl_gtf",
                "annotation": {"source": "TestAnnotation", "version": "test", "path": str(gtf)},
                "output_ids": {"ids": ["ENSG00000141510"]},
            }
        ],
        "ncbi_gene": {"enabled": False},
        "ensembl_archive": {"enabled": False},
    }

    with pytest.raises(ValueError, match="use 'output_genes' instead"):
        build(config=cfg, outdir=tmp_path / "out", workers=1)


def test_gene_filter_intersects_annotation_namespace(tmp_path):
    source = tmp_path / "sources"
    outdir = tmp_path / "out"
    source_tsv = tmp_path / "source.tsv"
    gtf = tmp_path / "genes.gtf"
    gene_filter = tmp_path / "filter.tsv"
    _write_hgnc(source)
    _write_source_table(source_tsv)
    _write_gtf(gtf)
    gene_filter.write_text("feature\nENSG000002.1\nENSG00000499999.1\n", encoding="utf-8")

    cfg = {
        "project": {"source_dir": str(source)},
        "outputs": {"include_symbols": False},
        "sources": [
            {
                "type": "msigdb_tsv",
                "path": str(source_tsv),
                "source_tag": "REACTOME",
                "collection": "C2",
                "family": "biology_process_pathway",
                "aspect": "pathway",
            }
        ],
        "targets": [
            {
                "name": "filtered_annotation",
                "type": "ensembl_gtf",
                "annotation": {"source": "TestAnnotation", "version": "test", "path": str(gtf)},
                "gene_filter": {"path": str(gene_filter), "id_column": "feature"},
            }
        ],
        "filtering": {"min_size_default": 1, "max_size_default": 10},
        "ncbi_gene": {"enabled": False},
        "ensembl_archive": {"enabled": False},
    }

    build(config=cfg, outdir=outdir, workers=1)

    gmt_line = (outdir / "gmt" / "filtered_annotation.gmt").read_text(encoding="utf-8").strip().split("\t")
    assert gmt_line[2:] == ["ENSG000002"]
    audit = pd.read_csv(outdir / "metadata" / "target_gene_filter.tsv.gz", sep="\t", dtype=str)
    assert audit["ensembl_gene_id"].tolist() == ["ENSG000002", "ENSG00000499999"]
    assert audit["has_annotation_metadata"].tolist() == ["True", "False"]
    summary = pd.read_csv(outdir / "qc" / "target_namespace_summary.tsv", sep="\t")
    assert summary.loc[0, "gene_filter_size"] == 2
    assert summary.loc[0, "effective_target_size"] == 1


def test_output_genes_symbol_column_augments_annotation_helper(tmp_path):
    source = tmp_path / "sources"
    outdir = tmp_path / "out"
    source_tsv = tmp_path / "source.tsv"
    gtf = tmp_path / "genes.gtf"
    output_genes = tmp_path / "output_genes.tsv"
    _write_hgnc(source)
    pd.DataFrame(
        {
            "gs_name": ["MATRIX_TERM"],
            "gs_description": ["desc"],
            "gene_symbol": ["MATRIXONLY"],
            "ensembl_gene": [""],
        }
    ).to_csv(source_tsv, sep="\t", index=False)
    _write_gtf(gtf)
    output_genes.write_text("feature_id\tgene_symbol\nENSG00000999999.1\tMATRIXONLY\n", encoding="utf-8")

    cfg = {
        "project": {"source_dir": str(source)},
        "outputs": {"include_symbols": False},
        "sources": [{"type": "msigdb_tsv", "path": str(source_tsv), "source_tag": "LOCAL"}],
        "targets": [
            {
                "name": "study_ids",
                "type": "ensembl_gtf",
                "annotation": {"source": "TestAnnotation", "version": "test", "path": str(gtf)},
                "output_genes": {"path": str(output_genes), "id_column": "feature_id", "symbol_column": "gene_symbol"},
            }
        ],
        "filtering": {"min_size_default": 1, "max_size_default": 10},
        "ncbi_gene": {"enabled": False},
        "ensembl_archive": {"enabled": False},
    }

    build(config=cfg, outdir=outdir, workers=1)

    assert "ENSG00000999999" in (outdir / "gmt" / "study_ids.gmt").read_text(encoding="utf-8")
    mapping = pd.read_csv(outdir / "metadata" / "gene_mapping_study_ids.tsv.gz", sep="\t", dtype=str).fillna("")
    assert mapping.set_index("ensembl_gene_id").loc["ENSG00000999999", "target_id_metadata_source"] == "output_genes"
    audit = pd.read_csv(outdir / "metadata" / "target_output_genes.tsv.gz", sep="\t", dtype=str).fillna("")
    assert audit.loc[0, "gene_symbol"] == "MATRIXONLY"
    assert audit.loc[0, "has_annotation_metadata"] == "False"
    warnings = (outdir / "qc" / "warnings.txt").read_text(encoding="utf-8")
    assert "1 output_genes lack annotation-helper metadata" in warnings


def test_metadata_fallback_fills_missing_output_metadata_without_overriding_primary(tmp_path):
    source = tmp_path / "sources"
    outdir = tmp_path / "out"
    source_tsv = tmp_path / "source.tsv"
    primary_gtf = tmp_path / "primary.gtf"
    fallback_gtf = tmp_path / "fallback.gtf"
    output_genes = tmp_path / "output_genes.tsv"
    _write_hgnc(source)
    pd.DataFrame(
        {
            "gs_name": ["SUPP_TERM", "SUPP_TERM"],
            "gs_description": ["desc", "desc"],
            "gene_symbol": ["DUP", "FILL"],
            "ensembl_gene": ["", ""],
        }
    ).to_csv(source_tsv, sep="\t", index=False)
    primary_gtf.write_text(
        'chr1\tALIEN\tgene\t1\t10\t.\t+\t.\tgene_id "ENSG000001.1"; gene_name "DUP"; gene_type "protein_coding";\n',
        encoding="utf-8",
    )
    fallback_gtf.write_text(
        'chr1\tALIEN\tgene\t20\t30\t.\t+\t.\tgene_id "ENSG000009.1"; gene_name "DUP"; gene_type "protein_coding";\n'
        'chr1\tALIEN\tgene\t40\t50\t.\t+\t.\tgene_id "ENSG000008.1"; gene_name "FILL"; gene_type "protein_coding";\n',
        encoding="utf-8",
    )
    output_genes.write_text("feature_id\nENSG000001.1\nENSG000009.1\nENSG000008.1\n", encoding="utf-8")

    cfg = {
        "project": {"source_dir": str(source)},
        "outputs": {"include_symbols": False},
        "sources": [{"type": "msigdb_tsv", "path": str(source_tsv), "source_tag": "LOCAL"}],
        "targets": [
            {
                "name": "fallback_target",
                "type": "ensembl_gtf",
                "annotation": {
                    "source": "Primary",
                    "version": "test",
                    "path": str(primary_gtf),
                    "metadata_fallbacks": [
                        {
                            "source": "Fallback",
                            "version": "test",
                            "path": str(fallback_gtf),
                        }
                    ],
                },
                "output_genes": {"path": str(output_genes), "id_column": "feature_id"},
            }
        ],
        "filtering": {"min_size_default": 1, "max_size_default": 10},
        "ncbi_gene": {"enabled": False},
        "ensembl_archive": {"enabled": False},
    }

    build(config=cfg, outdir=outdir, workers=1)

    gmt_line = (outdir / "gmt" / "fallback_target.gmt").read_text(encoding="utf-8").strip().split("\t")
    assert gmt_line[2:] == ["ENSG000001", "ENSG000008"]
    mapping = pd.read_csv(outdir / "metadata" / "gene_mapping_fallback_target.tsv.gz", sep="\t", dtype=str).fillna("")
    sources = mapping.set_index("ensembl_gene_id")["target_id_metadata_source"].to_dict()
    assert sources["ENSG000001"] == "annotation_gtf"
    assert sources["ENSG000008"] == "metadata_fallback"
    assert sources["ENSG000009"] == "metadata_fallback"
    fallback_audit = pd.read_csv(outdir / "metadata" / "target_metadata_fallbacks.tsv", sep="\t")
    assert fallback_audit.loc[0, "n_rows_added"] == 2
    assert fallback_audit.loc[0, "n_missing_after"] == 0
    summary = pd.read_csv(outdir / "qc" / "target_namespace_summary.tsv", sep="\t")
    assert summary.loc[0, "annotation_metadata_coverage"] == 1.0
    assert summary.loc[0, "metadata_fallback_helper_size"] == 2


def test_output_genes_and_gene_filter_are_mutually_exclusive(tmp_path):
    source = tmp_path / "sources"
    source_tsv = tmp_path / "source.tsv"
    gtf = tmp_path / "genes.gtf"
    _write_hgnc(source)
    _write_source_table(source_tsv)
    _write_gtf(gtf)

    cfg = {
        "project": {"source_dir": str(source)},
        "outputs": {"include_symbols": False},
        "sources": [{"type": "msigdb_tsv", "path": str(source_tsv), "source_tag": "REACTOME"}],
        "targets": [
            {
                "name": "bad_target",
                "type": "ensembl_gtf",
                "annotation": {"source": "TestAnnotation", "version": "test", "path": str(gtf)},
                "output_genes": {"ids": ["ENSG00000141510"]},
                "gene_filter": {"ids": ["ENSG00000141510"]},
            }
        ],
        "ncbi_gene": {"enabled": False},
        "ensembl_archive": {"enabled": False},
    }

    with pytest.raises(ValueError, match="only one of output_genes or gene_filter"):
        build(config=cfg, outdir=tmp_path / "out", workers=1)


def test_build_errors_and_writes_audit_for_conflicting_term_ids(tmp_path):
    source_a = tmp_path / "a.tsv"
    source_b = tmp_path / "b.tsv"
    outdir = tmp_path / "out"
    pd.DataFrame({"term_id": ["SHARED_TERM"], "gene_symbol": ["TP53"]}).to_csv(source_a, sep="\t", index=False)
    pd.DataFrame({"term_id": ["SHARED_TERM"], "gene_symbol": ["GENE2"]}).to_csv(source_b, sep="\t", index=False)

    cfg = {
        "project": {"source_dir": str(tmp_path / "sources")},
        "outputs": {"include_symbols": True},
        "sources": [
            {"type": "canonical_tsv", "path": str(source_a), "source": "Library A", "source_tag": "LIB_A", "collection": "a"},
            {"type": "canonical_tsv", "path": str(source_b), "source": "Library B", "source_tag": "LIB_B", "collection": "b"},
        ],
        "targets": [],
    }

    with pytest.raises(RuntimeError, match="conflicting term/source identities"):
        build(config=cfg, outdir=outdir, workers=1)

    audit_path = outdir / "metadata" / "term_id_collisions.tsv"
    assert audit_path.exists()
    audit = pd.read_csv(audit_path, sep="\t", dtype=str)
    assert audit["term_id"].unique().tolist() == ["SHARED_TERM"]
    assert sorted(audit["source_tag"].tolist()) == ["LIB_A", "LIB_B"]
    assert audit["action"].unique().tolist() == ["error"]


def _write_hgnc(source):
    hgnc = source / "hgnc"
    hgnc.mkdir(parents=True)
    pd.DataFrame(
        {
            "symbol": ["TP53", "GENE2", "GENE3", "MATRIXONLY", "DUP", "FILL"],
            "prev_symbol": ["", "", "", "", "", ""],
            "alias_symbol": ["", "", "", "", "", ""],
            "ensembl_gene_id": ["ENSG00000141510", "ENSG000002", "ENSG000003", "ENSG00000999999", "", ""],
        }
    ).to_csv(hgnc / "hgnc_complete_set.txt", sep="\t", index=False)


def _write_source_table(path):
    pd.DataFrame(
        {
            "gs_name": ["REACTOME_TERM", "REACTOME_TERM"],
            "gs_description": ["desc", "desc"],
            "gene_symbol": ["TP53", "GENE2"],
            "ensembl_gene": ["ENSG00000141510.18", "ENSG000002.1"],
        }
    ).to_csv(path, sep="\t", index=False)


def _write_gtf(path):
    path.write_text(
        'chr1\tALIEN\tgene\t1\t10\t.\t+\t.\tgene_id "ENSG00000141510.18"; gene_name "TP53"; gene_type "protein_coding";\n'
        'chr1\tALIEN\tgene\t20\t30\t.\t+\t.\tgene_id "ENSG000002.1"; gene_name "GENE2"; gene_type "protein_coding";\n'
        'chr1\tALIEN\tgene\t40\t50\t.\t+\t.\tgene_id "ENSG000003.1"; gene_name "GENE3"; gene_type "protein_coding";\n',
        encoding="utf-8",
    )
