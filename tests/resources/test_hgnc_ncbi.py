from __future__ import annotations

import gzip
from pathlib import Path

import pytest

from alien.cache import ResourceCache
from alien.config.schema import HgncConfig, NcbiGeneConfig
from alien.errors import ResourceError, SourceError
from alien.resources.hgnc import HGNC_FALLBACK_URLS, HGNC_URL, build_hgnc_maps, hgnc_path, hgnc_urls, load_hgnc_table
from alien.resources.ncbi import load_ncbi_maps

from conftest import FakeServer

HGNC_TABLE = "symbol\tprev_symbol\talias_symbol\tensembl_gene_id\nTP53\tP53OLD\tLFS1\tENSG00000141510\n"


def test_hgnc_falls_back_to_mirror(tmp_path: Path, fake_http: FakeServer) -> None:
    fake_http.routes[HGNC_URL] = (503, b"down")
    fake_http.routes[HGNC_FALLBACK_URLS[0]] = HGNC_TABLE
    table = load_hgnc_table(ResourceCache(tmp_path), HgncConfig())
    maps = build_hgnc_maps(table)
    assert maps.previous == {"P53OLD": "TP53"} and maps.alias == {"LFS1": "TP53"}
    assert maps.current_ensembl_ids == {"ENSG00000141510"}
    assert fake_http.urls() == [HGNC_URL, HGNC_FALLBACK_URLS[0]]


def test_hgnc_all_urls_failing_is_a_resource_error(tmp_path: Path, fake_http: FakeServer) -> None:
    with pytest.raises(ResourceError, match="every URL failed"):
        load_hgnc_table(ResourceCache(tmp_path), HgncConfig())
    assert not list((tmp_path / "hgnc").glob("*.txt"))


def test_hgnc_pinned_archive_date(tmp_path: Path, fake_http: FakeServer) -> None:
    config = HgncConfig(archive_date="2024-01-01")
    assert all("2024-01-01" in url for url in hgnc_urls(config))
    assert hgnc_path(tmp_path, config).name == "hgnc_complete_set_2024-01-01.txt"
    fake_http.routes[hgnc_urls(config)[0]] = HGNC_TABLE
    cache = ResourceCache(tmp_path)
    load_hgnc_table(cache, config)
    assert cache.used["hgnc"]["upstream_version"] == "2024-01-01"
    assert hgnc_urls(HgncConfig(url="https://mirror.test/h.txt")) == ["https://mirror.test/h.txt"]


def test_hgnc_table_without_symbol_column(tmp_path: Path) -> None:
    path = hgnc_path(tmp_path, HgncConfig())
    path.parent.mkdir(parents=True)
    path.write_text("x\ty\n1\t2\n", encoding="utf-8")
    with pytest.raises(SourceError, match="not an HGNC complete set"):
        load_hgnc_table(ResourceCache(tmp_path), HgncConfig())


def _write_ncbi(root: Path) -> None:
    directory = root / "ncbi_gene"
    directory.mkdir(parents=True)
    info = "#tax_id\tGeneID\tSymbol\tSynonyms\tdbXrefs\n9606\t7\tGGG7\tSYN7|-\tEnsembl:ENSG7\n9606\t8\tHHH8\tSHARED\t-\n9606\t9\tIII9\tSHARED\t-\n10090\t1\tMouse\t-\t-\n"
    history = "#tax_id\tGeneID\tDiscontinued_GeneID\tDiscontinued_Symbol\n9606\t7\t70\tOLD7\n9606\t-\t71\tGONE\n9606\t99\t72\tORPHAN\n"
    (directory / "Homo_sapiens.gene_info.gz").write_bytes(gzip.compress(info.encode()))
    (directory / "gene_history.gz").write_bytes(gzip.compress(history.encode()))


def test_ncbi_maps_and_cache_key_includes_synonym_option(tmp_path: Path) -> None:
    _write_ncbi(tmp_path)
    cache = ResourceCache(tmp_path)
    with_synonyms = load_ncbi_maps(NcbiGeneConfig(), cache)
    history = with_synonyms.stages["ncbi_gene_history"]
    assert set(history) == {"OLD7"} and history["OLD7"].ensembl_gene_ids == ("ENSG7",)
    assert set(with_synonyms.stages["ncbi_gene_info_synonym"]) == {"SYN7"}
    assert with_synonyms.ambiguous["ncbi_gene_info_synonym"] == {"SHARED": ["HHH8", "III9"]}
    # 10.2: changing use_gene_info_synonyms must not reuse the other variant's cache.
    without = load_ncbi_maps(NcbiGeneConfig(use_gene_info_synonyms=False), cache)
    assert "ncbi_gene_info_synonym" not in without.stages
    assert len(list((tmp_path / "ncbi_gene").glob("ncbi_gene_rescue_maps__*.json"))) == 2
    assert load_ncbi_maps(NcbiGeneConfig(), cache) == with_synonyms


def test_ncbi_disabled_returns_empty_maps(tmp_path: Path) -> None:
    assert not load_ncbi_maps(NcbiGeneConfig(enabled=False), ResourceCache(tmp_path))
