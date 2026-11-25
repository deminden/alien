from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

from alien.cache import ResourceCache
from alien.config import load_config
from alien.config.schema import CanonicalTsvSource, MsigdbTsvSource, SymbolGmtSource
from alien.errors import ConfigError, SourceError
from alien.sources import load_source_memberships, registered_types, term_id_collision_audit
from alien.sources.local import read_canonical_memberships, read_symbol_gmt, source_paths
from alien.sources.msigdb import read_msigdb_like_table


def test_canonical_table_normalization_with_source_gene_id(tmp_path: Path) -> None:
    path = tmp_path / "canonical.tsv"
    pd.DataFrame(
        {
            "term_id": ["Term 1", "Term 1"],
            "gene_symbol": ["TP53", "GENE2"],
            "gene_id": ["ENSG00000141510.18", "ENSG000002.1"],
            "gene_id_namespace": ["ensembl_gene", "ensembl_gene"],
        }
    ).to_csv(path, sep="\t", index=False)
    spec = CanonicalTsvSource(type="canonical_tsv", source="Example", source_tag="EX", family="biology_process_pathway")
    frame = read_canonical_memberships(path, spec)
    assert frame["term_id"].tolist() == ["Term_1", "Term_1"]
    assert frame["source_tag"].unique().tolist() == ["EX"]
    assert frame["gene_ensembl_from_source"].tolist() == ["ENSG00000141510.18", "ENSG000002.1"]


def test_canonical_gene_id_without_namespace_keeps_only_ensembl(tmp_path: Path) -> None:
    path = tmp_path / "canonical.tsv"
    pd.DataFrame({"term_id": ["T", "T"], "gene_symbol": ["A", "B"], "gene_id": ["ENSG1.1", "7157"]}).to_csv(
        path, sep="\t", index=False
    )
    frame = read_canonical_memberships(path, CanonicalTsvSource(type="canonical_tsv"))
    assert frame["gene_ensembl_from_source"].tolist() == ["ENSG1.1", ""]
    assert frame["source"].unique().tolist() == ["canonical"]
    assert frame["source_tag"].unique().tolist() == ["CANONICAL"]


def test_canonical_table_requires_term_and_symbol(tmp_path: Path) -> None:
    path = tmp_path / "bad.tsv"
    pd.DataFrame({"term_id": ["T"]}).to_csv(path, sep="\t", index=False)
    with pytest.raises(SourceError, match="gene_symbol"):
        read_canonical_memberships(path, CanonicalTsvSource(type="canonical_tsv"))


def test_msigdb_like_table_normalization(tmp_path: Path) -> None:
    path = tmp_path / "msigdb.tsv"
    pd.DataFrame(
        {
            "gs_name": ["PATHWAY_A", "PATHWAY_A"],
            "gs_description": ["desc", "desc"],
            "gene_symbol": ["TP53", "GENE2"],
            "ensembl_gene": ["ENSG00000141510.18", "ENSG000002.1"],
        }
    ).to_csv(path, sep="\t", index=False)
    spec = MsigdbTsvSource(type="msigdb_tsv", source_tag="REACTOME", collection="C2", aspect="pathway")
    frame = read_msigdb_like_table(path, spec)
    assert frame["term_id"].unique().tolist() == ["MSIGDB_REACTOME__PATHWAY_A"]
    assert frame["gene_symbol"].tolist() == ["TP53", "GENE2"]
    assert frame["gene_ensembl_from_source"].tolist() == ["ENSG00000141510.18", "ENSG000002.1"]
    with pytest.raises(SourceError, match="must contain"):
        read_msigdb_like_table(path, MsigdbTsvSource(type="msigdb_tsv", term_name_column="nope"))


def test_symbol_gmt_source_normalization(tmp_path: Path) -> None:
    path = tmp_path / "library.gmt"
    path.write_text("TERM_A\tdesc\tTP53\tGENE2\n", encoding="utf-8")
    spec = SymbolGmtSource(
        type="symbol_gmt", source="Example", source_tag="EX", family="disease_phenotype", aspect="disease"
    )
    frame = read_symbol_gmt(path, spec)
    assert frame["term_id"].unique().tolist() == ["EX__TERM_A"]
    assert frame["gene_symbol"].tolist() == ["TP53", "GENE2"]
    assert frame["gene_ensembl_from_source"].tolist() == ["", ""]
    assert frame["collection"].unique().tolist() == ["library"]


def test_symbol_gmt_without_paths_reads_named_directory(tmp_path: Path) -> None:
    library = tmp_path / "mylib"
    library.mkdir()
    (library / "b.gmt").write_text("B\td\tX\n", encoding="utf-8")
    (library / "a.gmt").write_text("A\td\tY\n", encoding="utf-8")
    spec = SymbolGmtSource(type="symbol_gmt", name="mylib")
    assert source_paths(spec, tmp_path) == [library / "a.gmt", library / "b.gmt"]
    with pytest.raises(SourceError, match="does not exist"):
        source_paths(SymbolGmtSource(type="symbol_gmt", path=tmp_path / "missing.gmt"), tmp_path)


def test_disabled_sources_are_skipped_and_empty_builds_fail(tmp_path: Path) -> None:
    gmt = tmp_path / "x.gmt"
    gmt.write_text("T\td\tA\n", encoding="utf-8")
    config = load_config({"sources": [{"type": "symbol_gmt", "path": str(gmt), "enabled": False}], "targets": []})
    with pytest.raises(SourceError, match="No source gene-set memberships"):
        load_source_memberships(config, ResourceCache(tmp_path))
    assert {"symbol_gmt", "msigdb_remote", "enrichr_remote", "plugin"} <= set(registered_types())


def test_term_id_collision_audit_ignores_rows_of_the_same_identity() -> None:
    frame = pd.DataFrame(
        {
            "term_id": ["TERM", "TERM"],
            "source_tag": ["LOCAL", "LOCAL"],
            "collection": ["example", "example"],
            "original_name": ["Term", "Term"],
            "gene_symbol": ["TP53", "GENE2"],
        }
    )
    assert term_id_collision_audit(frame).empty


def test_default_filter_lists_can_be_extended(tmp_path: Path) -> None:
    config = tmp_path / "config.yml"
    config.write_text(
        "gene_mapping:\n  extra_non_gene_tokens: [signature]\nfiltering:\n  extra_broad_disease_regex: ['^syndrome$']\n"
        "sources:\n  - type: symbol_gmt\n    path: example.gmt\ntargets: []\n",
        encoding="utf-8",
    )
    cfg = load_config(config)
    assert "signature" in cfg.gene_mapping.non_gene_tokens and "disease" in cfg.gene_mapping.non_gene_tokens
    assert "^syndrome$" in cfg.filtering.broad_disease_regex


def test_config_without_sources_is_an_error() -> None:
    with pytest.raises(ConfigError, match="at least one source"):
        load_config()
