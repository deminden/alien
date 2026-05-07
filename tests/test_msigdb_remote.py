import hashlib
import zipfile

import pandas as pd
import pytest
import rdata

from alien.sources import ensure_msigdb_remote_cache, read_msigdb_remote


def test_msigdb_remote_loads_human_database_and_expected_columns(tmp_path):
    spec = _tiny_release_spec(tmp_path)

    frame = read_msigdb_remote(spec, tmp_path)

    assert list(frame.columns) == [
        "term_id",
        "original_name",
        "display_name",
        "description",
        "source",
        "source_tag",
        "collection",
        "subcollection",
        "family",
        "aspect",
        "db_version",
        "source_url",
        "source_license_note",
        "gene_symbol",
        "gene_ensembl_from_source",
        "metadata_json",
    ]
    assert set(frame["source_tag"]) == {"GOMF", "REACTOME"}
    assert set(frame["gene_symbol"]) == {"TP53", "GENE2", "RPLP0", "KINASE1"}
    assert set(frame["gene_ensembl_from_source"]) == {"ENSG00000141510", "ENSG000002", "ENSG000003", "ENSG000004"}
    assert frame["db_version"].unique().tolist() == ["test.1"]


def test_msigdb_remote_collection_filter_matches_summary(tmp_path):
    spec = _tiny_release_spec(tmp_path, collection="C2")

    frame = read_msigdb_remote(spec, tmp_path)

    assert set(frame["collection"]) == {"C2"}
    assert set(frame["source_tag"]) == {"REACTOME"}
    assert set(frame["term_id"]) == {
        "MSIGDB_REACTOME__REACTOME_SIGNAL",
        "MSIGDB_REACTOME__REACTOME_METABOLISM",
    }


def test_msigdb_remote_subcollection_accepts_full_or_short_name(tmp_path):
    full = read_msigdb_remote(_tiny_release_spec(tmp_path, collection="C5", subcollection="GO:MF"), tmp_path)
    short = read_msigdb_remote(_tiny_release_spec(tmp_path, collection="C5", subcollection="MF"), tmp_path)

    assert full.equals(short)
    assert full["source_tag"].unique().tolist() == ["GOMF"]
    assert full["family"].unique().tolist() == ["biology_function_location"]
    assert full["aspect"].unique().tolist() == ["molecular_function"]


def test_msigdb_remote_db_species_is_case_insensitive_and_filters_database(tmp_path):
    frame = read_msigdb_remote(_tiny_release_spec(tmp_path, db_species="mm"), tmp_path)

    assert set(frame["collection"]) == {"M8"}
    assert set(frame["source_tag"]) == {"M8"}
    assert frame["gene_symbol"].tolist() == ["Trp53", "Gene2"]


def test_msigdb_remote_raises_for_unknown_collection_or_species(tmp_path):
    with pytest.raises(ValueError, match="No MSigDB remote collection files matched"):
        read_msigdb_remote(_tiny_release_spec(tmp_path, collection="X"), tmp_path)

    with pytest.raises(ValueError, match="No MSigDB remote collection files matched"):
        read_msigdb_remote(_tiny_release_spec(tmp_path, db_species="RN"), tmp_path)


def test_msigdb_remote_verifies_cached_archive_checksum(tmp_path):
    spec = _tiny_release_spec(tmp_path)
    spec["release"]["zip_md5"] = "0" * 32

    with pytest.raises(RuntimeError, match="MD5 verification"):
        read_msigdb_remote(spec, tmp_path)


def test_msigdb_remote_requires_summary_after_extraction(tmp_path):
    cache_dir = tmp_path / "remote"
    cache_dir.mkdir()
    zip_path = cache_dir / "msigdb.no_summary.zip"
    payload = tmp_path / "payload.txt"
    payload.write_text("not an RDS file\n", encoding="utf-8")
    with zipfile.ZipFile(zip_path, "w") as archive:
        archive.write(payload, arcname="payload.txt")
    release = {
        "zip_url": "https://example.test/msigdb.no_summary.zip",
        "zip_md5": hashlib.md5(zip_path.read_bytes()).hexdigest(),
        "zip_name": zip_path.name,
        "summary_rds": "msigdb.no_summary.summary.rds",
    }

    with pytest.raises(FileNotFoundError, match="summary RDS"):
        ensure_msigdb_remote_cache(cache_dir, release)


def test_msigdb_remote_force_reextracts_cached_release(tmp_path):
    spec = _tiny_release_spec(tmp_path)
    release = spec["release"]
    cache_dir = spec["cache_dir"]
    release_dir = ensure_msigdb_remote_cache(cache_dir, release)
    stale = release_dir / "stale.txt"
    stale.write_text("old extraction\n", encoding="utf-8")

    refreshed = ensure_msigdb_remote_cache(cache_dir, release, force=True)

    assert refreshed == release_dir
    assert not stale.exists()
    assert (release_dir / release["summary_rds"]).exists()


def _tiny_release_spec(tmp_path, **overrides):
    cache_dir = tmp_path / "remote"
    release_dir = tmp_path / "release_payload"
    release_dir.mkdir(exist_ok=True)
    files = _write_tiny_rds_release(release_dir)
    zip_path = cache_dir / "msigdb.test.zip"
    cache_dir.mkdir(exist_ok=True)
    with zipfile.ZipFile(zip_path, "w") as archive:
        for path in files:
            archive.write(path, arcname=path.name)
    spec = {
        "cache_dir": cache_dir,
        "db_species": "HS",
        "release": {
            "zip_url": "https://example.test/msigdb.test.zip",
            "zip_md5": hashlib.md5(zip_path.read_bytes()).hexdigest(),
            "zip_name": zip_path.name,
            "summary_rds": "msigdb.test.summary.rds",
        },
    }
    spec.update(overrides)
    return spec


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
    for path, frame in zip(paths, [summary, hs_c2, hs_c5, mm_m8]):
        rdata.write_rds(path, frame)
    return paths
