import pandas as pd

from alien.sources import read_gene_universe


def test_read_gene_universe_from_parquet_index(tmp_path):
    path = tmp_path / "matrix.parquet"
    frame = pd.DataFrame(
        {"Description": ["TP53", "GENE2"], "sample_a": [1, 0]},
        index=pd.Index(["ENSG00000141510.18", "ENSG000002.1"], name="Name"),
    )
    frame.to_parquet(path)

    universe, rows, id_type = read_gene_universe(path)

    assert universe == {"ENSG00000141510", "ENSG000002"}
    assert id_type == "ensembl_versioned"
    assert rows["input_gene_id"].tolist() == ["ENSG00000141510.18", "ENSG000002.1"]


def test_read_gene_universe_from_tsv_gz_named_column(tmp_path):
    path = tmp_path / "expression.tsv.gz"
    pd.DataFrame(
        {
            "feature": ["ENSG00000141510.18", "ENSG000002.1"],
            "sample_a": [10, 5],
            "sample_b": [0, 3],
        }
    ).to_csv(path, sep="\t", index=False, compression="gzip")

    universe, rows, id_type = read_gene_universe(path, column="feature")

    assert universe == {"ENSG00000141510", "ENSG000002"}
    assert id_type == "ensembl_versioned"
    assert rows["input_gene_id"].tolist() == ["ENSG00000141510.18", "ENSG000002.1"]


def test_read_gene_universe_unknown_column_errors(tmp_path):
    path = tmp_path / "expression.tsv.gz"
    pd.DataFrame({"feature": ["ENSG00000141510.18"]}).to_csv(path, sep="\t", index=False, compression="gzip")

    try:
        read_gene_universe(path, column="missing")
    except ValueError as exc:
        assert "missing" in str(exc)
    else:
        raise AssertionError("Expected a missing-column error.")


def test_read_gene_universe_from_inline_ids():
    universe, rows, id_type = read_gene_universe(None, ids=["ENSG00000141510.18", "ENSG000002"])

    assert universe == {"ENSG00000141510", "ENSG000002"}
    assert id_type == "mixed"
    assert rows["input_gene_id"].tolist() == ["ENSG00000141510.18", "ENSG000002"]
