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
    universe = tmp_path / "universe.tsv"
    _write_hgnc(source)
    _write_source_table(source_tsv)
    symbol_gmt.write_text("SYMBOL_TERM\tdesc\tTP53\tGENE3\n", encoding="utf-8")
    _write_gtf(gtf)
    universe.write_text("feature\nENSG00000141510.18\nENSG000002.1\nENSG000003.1\n", encoding="utf-8")

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
                "restrict_to": {"path": str(universe), "column": "feature"},
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
    assert (out1 / "metadata" / "term_manifest.tsv.gz").exists()
    assert (out1 / "qc" / "warnings.txt").exists()
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
            "symbol": ["TP53", "GENE2", "GENE3"],
            "prev_symbol": ["", "", ""],
            "alias_symbol": ["", "", ""],
            "ensembl_gene_id": ["ENSG00000141510", "ENSG000002", "ENSG000003"],
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
