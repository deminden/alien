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
