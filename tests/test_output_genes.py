import pandas as pd
import pytest

from alien.sources import read_output_genes


def test_read_output_genes_from_parquet_index(tmp_path):
    path = tmp_path / "matrix.parquet"
    frame = pd.DataFrame(
        {"Description": ["TP53", "GENE2"], "sample_a": [1, 0]},
        index=pd.Index(["ENSG00000141510.18", "ENSG000002.1"], name="Name"),
    )
    frame.to_parquet(path)

    output_genes, rows, id_type = read_output_genes(path)

    assert output_genes == {"ENSG00000141510", "ENSG000002"}
    assert id_type == "ensembl_versioned"
    assert rows["input_gene_id"].tolist() == ["ENSG00000141510.18", "ENSG000002.1"]


def test_read_output_genes_from_tsv_gz_named_column(tmp_path):
    path = tmp_path / "expression.tsv.gz"
    pd.DataFrame(
        {
            "feature": ["ENSG00000141510.18", "ENSG000002.1"],
            "sample_a": [10, 5],
            "sample_b": [0, 3],
        }
    ).to_csv(path, sep="\t", index=False, compression="gzip")

    output_genes, rows, id_type = read_output_genes(path, id_column="feature")

    assert output_genes == {"ENSG00000141510", "ENSG000002"}
    assert id_type == "ensembl_versioned"
    assert rows["input_gene_id"].tolist() == ["ENSG00000141510.18", "ENSG000002.1"]


def test_read_output_genes_with_symbol_column(tmp_path):
    path = tmp_path / "expression.tsv.gz"
    pd.DataFrame(
        {
            "feature_id": ["ENSG00000141510.18", "ENSG000002.1"],
            "gene_symbol": ["TP53", "GENE2"],
            "sample_a": [10, 5],
        }
    ).to_csv(path, sep="\t", index=False, compression="gzip")

    output_genes, rows, id_type = read_output_genes(path, id_column="feature_id", symbol_column="gene_symbol")

    assert output_genes == {"ENSG00000141510", "ENSG000002"}
    assert id_type == "ensembl_versioned"
    assert rows["gene_symbol"].tolist() == ["TP53", "GENE2"]


def test_read_output_genes_unknown_column_errors(tmp_path):
    path = tmp_path / "expression.tsv.gz"
    pd.DataFrame({"feature": ["ENSG00000141510.18"]}).to_csv(path, sep="\t", index=False, compression="gzip")

    try:
        read_output_genes(path, id_column="missing")
    except ValueError as exc:
        assert "missing" in str(exc)
    else:
        raise AssertionError("Expected a missing-column error.")


def test_read_output_genes_rejects_symbol_only_id_column(tmp_path):
    path = tmp_path / "expression.tsv.gz"
    pd.DataFrame({"feature": ["TP53"], "gene_symbol": ["TP53"]}).to_csv(path, sep="\t", index=False, compression="gzip")

    with pytest.raises(ValueError, match="must provide Ensembl IDs"):
        read_output_genes(path, id_column="feature", symbol_column="gene_symbol")


def test_read_output_genes_from_inline_ids():
    output_genes, rows, id_type = read_output_genes(None, ids=["ENSG00000141510.18", "ENSG000002"])

    assert output_genes == {"ENSG00000141510", "ENSG000002"}
    assert id_type == "mixed"
    assert rows["input_gene_id"].tolist() == ["ENSG00000141510.18", "ENSG000002"]
