"""5.12, 5.13 and 5.14: honest columns, the right cache directory, and clean output-gene tables."""

from __future__ import annotations

import logging
from pathlib import Path

import pytest

from alien import build

from conftest import MiniWorld, read_tsv


def test_size_filter_table_has_no_duplicated_size_columns(mini_world: MiniWorld) -> None:
    gmt = mini_world.gmt("lib.gmt", {"T1": ["G1"], "T2": ["G2", "G3", "G4"]})
    result = build(config=mini_world.config(gmt))
    removed = read_tsv(result.outdir / "metadata" / "removed_terms_size_filter.tsv.gz")
    assert "mapped_size_before_output_genes" not in removed.columns
    assert "mapped_size" in removed.columns


def test_mapping_summary_counts_size_losses_per_source_tag(mini_world: MiniWorld) -> None:
    gmt_a = mini_world.gmt("a.gmt", {"SMALL": ["G1"], "OK": ["G2", "G3"]})
    gmt_b = mini_world.gmt("b.gmt", {"OK2": ["G4", "G5", "G6"]})
    cfg = mini_world.config(gmt_a)
    cfg["sources"].append({"type": "symbol_gmt", "path": str(gmt_b), "source_tag": "OTHER"})
    result = build(config=cfg)
    summary = read_tsv(result.outdir / "qc" / "mapping_summary.tsv")
    mini = summary[summary["target_namespace"].eq("mini")].set_index("source_tag")
    assert mini.loc["LOCAL", "n_terms_lost_because_too_small"] == "1"
    assert mini.loc["OTHER", "n_terms_lost_because_too_small"] == "0"


def test_archive_cache_defaults_to_alien_sources(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    from alien.mapping import _build_archive_resolver

    monkeypatch.chdir(tmp_path)
    resolver = _build_archive_resolver(None, {"ensembl_archive": {"enabled": True}})
    assert resolver is not None
    assert resolver.cache_path == Path("data/alien_sources/ensembl_archive/archive_id_cache.json")


def test_symbol_rows_in_output_genes_do_not_leak_into_ensembl_column(caplog: pytest.LogCaptureFixture) -> None:
    from alien.sources import read_output_genes

    with caplog.at_level(logging.WARNING, logger="alien"):
        ids, table, _ = read_output_genes(None, ids=["ENSG00000000101.2", "TP53", "ENSG00000000102"])
    assert ids == {"ENSG00000000101", "ENSG00000000102"}
    symbol_row = table[table["input_gene_id"].eq("TP53")].iloc[0]
    assert symbol_row["ensembl_gene_id"] == ""
    assert symbol_row["id_type"] == "symbol"
    assert "1 output gene row" in caplog.text
