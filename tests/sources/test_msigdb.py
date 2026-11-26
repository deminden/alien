from __future__ import annotations

import hashlib
import zipfile
from pathlib import Path
from typing import Any

import pandas as pd
import pytest
import rdata

from alien.cache import ResourceCache
from alien.config.schema import MsigdbRemoteSource
from alien.errors import ConfigError, ResourceError, SourceError
from alien.schema import CANONICAL_COLUMNS
from alien.sources.msigdb import ensure_release, load_msigdb_remote, release_info, remote_source_tag


def _load(tmp_path: Path, **overrides: Any) -> pd.DataFrame:
    return load_msigdb_remote(_tiny_release_spec(tmp_path, **overrides), ResourceCache(tmp_path))


def test_msigdb_remote_loads_human_database_and_expected_columns(tmp_path: Path) -> None:
    frame = _load(tmp_path)
    assert list(frame.columns) == list(CANONICAL_COLUMNS)
    assert set(frame["source_tag"]) == {"GOMF", "REACTOME"}
    assert set(frame["gene_symbol"]) == {"TP53", "GENE2", "RPLP0", "KINASE1"}
    assert set(frame["gene_ensembl_from_source"]) == {"ENSG00000141510", "ENSG000002", "ENSG000003", "ENSG000004"}
    assert frame["db_version"].unique().tolist() == ["test.1"]
    # The second load comes from the derived parquet cache and is identical.
    assert _load(tmp_path).equals(frame)


def test_msigdb_remote_parses_each_rds_once(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    from alien.sources import msigdb

    parsed: list[str] = []
    original = msigdb.read_rds_frame

    def counting(path: Path) -> pd.DataFrame:
        parsed.append(path.name)
        return original(path)

    monkeypatch.setattr(msigdb, "read_rds_frame", counting)
    first = _load(tmp_path, collection="C2")
    # A second source over the same collection file (different memberships cache key) reuses the parse.
    assert not _load(tmp_path, collection="C2", subcollection="CP:REACTOME").empty
    assert [name for name in parsed if "summary" not in name] == ["hs_c2.rds"]
    assert _load(tmp_path, collection="C2").equals(first)


def test_msigdb_remote_collection_filter_matches_summary(tmp_path: Path) -> None:
    frame = _load(tmp_path, collection="C2")
    assert set(frame["collection"]) == {"C2"}
    assert set(frame["term_id"]) == {"MSIGDB_REACTOME__REACTOME_SIGNAL", "MSIGDB_REACTOME__REACTOME_METABOLISM"}


def test_msigdb_remote_subcollection_accepts_full_or_short_name(tmp_path: Path) -> None:
    full = _load(tmp_path, collection="C5", subcollection="GO:MF")
    short = _load(tmp_path, collection="C5", subcollection="MF")
    assert full.equals(short)
    assert full["source_tag"].unique().tolist() == ["GOMF"]
    assert full["family"].unique().tolist() == ["biology_function_location"]
    assert full["aspect"].unique().tolist() == ["molecular_function"]


def test_msigdb_remote_db_species_is_case_insensitive_and_filters_database(tmp_path: Path) -> None:
    frame = _load(tmp_path, db_species="mm")
    assert set(frame["collection"]) == {"M8"}
    assert set(frame["source_tag"]) == {"M8"}
    assert frame["gene_symbol"].tolist() == ["Trp53", "Gene2"]


def test_msigdb_remote_raises_for_unknown_collection_or_species(tmp_path: Path) -> None:
    with pytest.raises(SourceError, match="No MSigDB remote collection files matched"):
        _load(tmp_path, collection="X")
    with pytest.raises(SourceError, match="No MSigDB remote collection files matched"):
        _load(tmp_path, db_species="RN")


def test_msigdb_remote_verifies_cached_archive_checksum(tmp_path: Path) -> None:
    spec = _tiny_release_spec(tmp_path)
    spec = spec.model_copy(update={"release": spec.release.model_copy(update={"zip_md5": "0" * 32})})
    with pytest.raises(ResourceError, match="MD5 verification"):
        load_msigdb_remote(spec, ResourceCache(tmp_path))


def test_msigdb_remote_requires_summary_after_extraction(tmp_path: Path) -> None:
    cache_dir = tmp_path / "remote"
    cache_dir.mkdir()
    zip_path = cache_dir / "msigdb.no_summary.zip"
    with zipfile.ZipFile(zip_path, "w") as archive:
        archive.writestr("payload.txt", "not an RDS file\n")
    spec = MsigdbRemoteSource(
        type="msigdb_remote",
        cache_dir=cache_dir,
        release={
            "zip_url": "https://example.test/msigdb.no_summary.zip",
            "zip_md5": hashlib.md5(zip_path.read_bytes()).hexdigest(),
            "zip_name": zip_path.name,
            "summary_rds": "msigdb.no_summary.summary.rds",
        },
    )
    with pytest.raises(ResourceError, match="summary"):
        ensure_release(spec, ResourceCache(tmp_path))


def test_msigdb_remote_refuses_zip_slip(tmp_path: Path) -> None:
    cache_dir = tmp_path / "remote"
    cache_dir.mkdir()
    zip_path = cache_dir / "evil.zip"
    with zipfile.ZipFile(zip_path, "w") as archive:
        archive.writestr("../escape.txt", "x")
    spec = MsigdbRemoteSource(
        type="msigdb_remote",
        cache_dir=cache_dir,
        release={"zip_url": "https://example.test/evil.zip", "zip_md5": hashlib.md5(zip_path.read_bytes()).hexdigest()},
    )
    with pytest.raises(ResourceError, match="outside"):
        ensure_release(spec, ResourceCache(tmp_path))
    assert not (tmp_path / "escape.txt").exists()


def test_msigdb_remote_force_reextracts_cached_release(tmp_path: Path) -> None:
    spec = _tiny_release_spec(tmp_path)
    release_dir = ensure_release(spec, ResourceCache(tmp_path))
    stale = release_dir / "stale.txt"
    stale.write_text("old extraction\n", encoding="utf-8")
    refreshed = ensure_release(spec.model_copy(update={"force": True}), ResourceCache(tmp_path))
    assert refreshed == release_dir
    assert not stale.exists()


def test_msigdb_remote_downloads_missing_archive(tmp_path: Path, fake_http: Any) -> None:
    spec = _tiny_release_spec(tmp_path)
    zip_path = tmp_path / "remote" / "msigdb.test.zip"
    fake_http.routes["https://example.test/msigdb.test.zip"] = zip_path.read_bytes()
    zip_path.unlink()
    frame = load_msigdb_remote(spec, ResourceCache(tmp_path / "fresh"))
    assert not frame.empty
    assert fake_http.urls() == ["https://example.test/msigdb.test.zip"]


def test_release_table_and_tags() -> None:
    builtin = release_info(MsigdbRemoteSource(type="msigdb_remote"))
    assert builtin["zip_name"] == "msigdb.2026.1.zip" and builtin["summary_rds"] == "msigdb.2026.1.summary.rds"
    with pytest.raises(ConfigError, match="Unknown MSigDB release"):
        release_info(MsigdbRemoteSource(type="msigdb_remote", version="1999.1"))
    assert (
        remote_source_tag("GO:BP") == "GOBP"
        and remote_source_tag("") == "MSIGDB"
        and remote_source_tag("CP:KEGG") == "KEGG"
    )


def _tiny_release_spec(tmp_path: Path, **overrides: Any) -> MsigdbRemoteSource:
    cache_dir = tmp_path / "remote"
    release_dir = tmp_path / "release_payload"
    release_dir.mkdir(exist_ok=True)
    files = _write_tiny_rds_release(release_dir)
    zip_path = cache_dir / "msigdb.test.zip"
    cache_dir.mkdir(exist_ok=True)
    if not zip_path.exists():
        with zipfile.ZipFile(zip_path, "w") as archive:
            for path in files:
                archive.write(path, arcname=path.name)
    data: dict[str, Any] = {
        "type": "msigdb_remote",
        "cache_dir": cache_dir,
        "db_species": "HS",
        "release": {
            "zip_url": "https://example.test/msigdb.test.zip",
            "zip_md5": hashlib.md5(zip_path.read_bytes()).hexdigest(),
            "zip_name": zip_path.name,
            "summary_rds": "msigdb.test.summary.rds",
        },
    }
    data.update(overrides)
    return MsigdbRemoteSource.model_validate(data)


def _write_tiny_rds_release(release_dir):
    summary = pd.DataFrame(
        {
            "db_version": ["test.1", "test.1", "test.1"],
            "db_target_species": ["HS", "HS", "MM"],
            "gs_collection": ["C2", "C5", "M8"],
            "gs_subcollection": ["CP:REACTOME", "GO:MF", "M8"],
            "gs_collection_name": ["Curated", "Ontology", "Mouse cell type"],
            "num_genesets": ["2", "1", "1"],
            "df_rds": ["hs_c2.rds", "hs_c5.rds", "mm_m8.rds"],
        }
    ).astype(object)
    hs_c2 = pd.DataFrame(
        {
            "db_version": ["test.1", "test.1", "test.1"],
            "db_target_species": ["HS", "HS", "HS"],
            "db_gene_symbol": ["TP53", "GENE2", "RPLP0"],
            "db_ncbi_gene": ["7157", "2", "6175"],
            "db_ensembl_gene": ["ENSG00000141510", "ENSG000002", "ENSG000003"],
            "source_gene": ["TP53", "GENE2", "RPLP0"],
            "gs_id": ["M1", "M1", "M2"],
            "gs_name": ["REACTOME_SIGNAL", "REACTOME_SIGNAL", "REACTOME_METABOLISM"],
            "gs_description": ["signal desc", "signal desc", "metabolism desc"],
            "gs_collection": ["C2", "C2", "C2"],
            "gs_subcollection": ["CP:REACTOME", "CP:REACTOME", "CP:REACTOME"],
            "gs_url": ["https://example.test/signal", "https://example.test/signal", "https://example.test/metabolism"],
        }
    ).astype(object)
    hs_c5 = pd.DataFrame(
        {
            "db_version": ["test.1"],
            "db_target_species": ["HS"],
            "db_gene_symbol": ["KINASE1"],
            "db_ncbi_gene": ["100"],
            "db_ensembl_gene": ["ENSG000004"],
            "source_gene": ["KINASE1"],
            "gs_id": ["M3"],
            "gs_name": ["GO_KINASE_ACTIVITY"],
            "gs_description": ["kinase desc"],
            "gs_collection": ["C5"],
            "gs_subcollection": ["GO:MF"],
            "gs_url": ["https://example.test/kinase"],
        }
    ).astype(object)
    mm_m8 = pd.DataFrame(
        {
            "db_version": ["test.1", "test.1"],
            "db_target_species": ["MM", "MM"],
            "db_gene_symbol": ["Trp53", "Gene2"],
            "db_ncbi_gene": ["22059", "2"],
            "db_ensembl_gene": ["ENSMUSG00000059552", "ENSMUSG00000000002"],
            "source_gene": ["Trp53", "Gene2"],
            "gs_id": ["MM1", "MM1"],
            "gs_name": ["MOUSE_CELL_TYPE", "MOUSE_CELL_TYPE"],
            "gs_description": ["mouse desc", "mouse desc"],
            "gs_collection": ["M8", "M8"],
            "gs_subcollection": ["M8", "M8"],
            "gs_url": ["https://example.test/mouse", "https://example.test/mouse"],
        }
    ).astype(object)
    paths = [
        release_dir / "msigdb.test.summary.rds",
        release_dir / "hs_c2.rds",
        release_dir / "hs_c5.rds",
        release_dir / "mm_m8.rds",
    ]
    for path, frame in zip(paths, [summary, hs_c2, hs_c5, mm_m8], strict=True):
        rdata.write_rds(path, frame)
    return paths
