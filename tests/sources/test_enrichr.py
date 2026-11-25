from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pandas as pd
import pytest

from alien import build
from alien.cache import LockFile, ResourceCache
from alien.config import load_config
from alien.config.schema import EnrichrLibrary, EnrichrRemoteSource
from alien.errors import ConfigError, ResourceError, SourceError
from alien.reports.manifest import build_source_manifest
from alien.sources.enrichr import load_enrichr_remote, resolve_enrichr_library

from conftest import FakeServer

INDEX = "https://example.test/enrichr/datasetStatistics"
LIBRARY = "https://example.test/enrichr/geneSetLibrary?mode=text&libraryName={library}"


def _serve(server: FakeServer, term: str = "TERM_A") -> None:
    names = ["KEGG_2021_Human", "ClinVar_2019", "ClinVar_2022"]
    server.routes[INDEX] = json.dumps({"statistics": [{"libraryName": name, "numTerms": 1} for name in names]})
    for name in ("KEGG_2021_Human", "ClinVar_2022", "ClinVar_2019"):
        server.routes[LIBRARY.format(library=name)] = f"{term}\tdesc\tTP53\tGENE2\n"


def _spec(tmp_path: Path, **overrides: Any) -> EnrichrRemoteSource:
    data: dict[str, Any] = {
        "type": "enrichr_remote",
        "cache_dir": tmp_path / "enrichr",
        "dataset_statistics_endpoint": INDEX,
        "library_download_endpoint": LIBRARY,
        "libraries": [
            {"name": "KEGG_2021_Human", "source_tag": "KEGG", "family": "biology_process_pathway", "aspect": "pathway"}
        ],
    }
    data.update(overrides)
    return EnrichrRemoteSource.model_validate(data)


def test_exact_name_downloads_and_normalizes(tmp_path: Path, fake_http: FakeServer) -> None:
    _serve(fake_http)
    frame = load_enrichr_remote(_spec(tmp_path), ResourceCache(tmp_path))
    assert fake_http.urls() == [INDEX, LIBRARY.format(library="KEGG_2021_Human")]
    assert frame["term_id"].unique().tolist() == ["ENRICHR_KEGG_2021_Human__TERM_A"]
    assert frame["source_tag"].unique().tolist() == ["KEGG"]
    assert frame["collection"].unique().tolist() == ["KEGG_2021_Human"]
    assert frame["gene_symbol"].tolist() == ["TP53", "GENE2"]
    metadata = json.loads(frame["metadata_json"].iloc[0])
    assert (metadata["selected_library"], metadata["match_method"]) == ("KEGG_2021_Human", "exact")


def test_reuses_cache_and_force_refreshes(tmp_path: Path, fake_http: FakeServer) -> None:
    _serve(fake_http)
    first = load_enrichr_remote(_spec(tmp_path), ResourceCache(tmp_path))
    fake_http.calls.clear()
    assert load_enrichr_remote(_spec(tmp_path), ResourceCache(tmp_path)).equals(first)
    assert fake_http.calls == []
    _serve(fake_http, term="TERM_B")
    refreshed = load_enrichr_remote(_spec(tmp_path), ResourceCache(tmp_path, force=True))
    assert len(fake_http.calls) == 2
    assert refreshed["term_id"].unique().tolist() == ["ENRICHR_KEGG_2021_Human__TERM_B"]


def test_regex_chooses_latest_year_then_shortest_name() -> None:
    names = ["ClinVar_2019", "ClinVar_2022", "ClinVar_Long_2022", "ClinVar_2021"]
    selected, candidates, method = resolve_enrichr_library(
        names, EnrichrLibrary(name="ClinVar", match=r"^ClinVar.*[0-9]{4}$")
    )
    assert selected == "ClinVar_2022"
    assert candidates == ["ClinVar_2019", "ClinVar_2021", "ClinVar_2022", "ClinVar_Long_2022"]
    assert method == "regex"
    with pytest.raises(SourceError, match="No Enrichr libraries matched"):
        resolve_enrichr_library(names, EnrichrLibrary(match="^Nope"))


def test_source_manifest_records_regex_resolution(tmp_path: Path, fake_http: FakeServer) -> None:
    _serve(fake_http)
    library = {
        "name": "ClinVar",
        "match": r"^ClinVar_[0-9]{4}$",
        "source_tag": "CLINVAR",
        "family": "disease_phenotype",
    }
    frame = load_enrichr_remote(_spec(tmp_path, libraries=[library]), ResourceCache(tmp_path))
    manifest = build_source_manifest(frame)
    assert manifest.loc[0, "selected_library"] == "ClinVar_2022"
    assert manifest.loc[0, "configured_name"] == "ClinVar"
    assert manifest.loc[0, "match_method"] == "regex"
    assert manifest.loc[0, "candidate_libraries"] == '["ClinVar_2019", "ClinVar_2022"]'


def test_frozen_build_keeps_the_pinned_library(tmp_path: Path, fake_http: FakeServer) -> None:
    _serve(fake_http)
    library = {"match": r"^ClinVar_[0-9]{4}$"}
    cache = ResourceCache(tmp_path, lock=LockFile(path=tmp_path / "alien.lock.json"))
    load_enrichr_remote(_spec(tmp_path, libraries=[library]), cache)
    cache.lock.resolutions["enrichr:library::^ClinVar_[0-9]{4}$"] = "ClinVar_2019"
    cache.lock.save()
    frozen = ResourceCache(tmp_path, frozen=True, lock=LockFile.load(tmp_path / "alien.lock.json"))
    with pytest.raises(ResourceError, match=r"ClinVar_2019\.gmt is not cached"):  # pinned choice wins over the regex
        load_enrichr_remote(_spec(tmp_path, libraries=[library]), frozen)


def test_unresolved_library_and_optional_key(tmp_path: Path, fake_http: FakeServer) -> None:
    _serve(fake_http)
    with pytest.raises(SourceError, match="was not found"):
        load_enrichr_remote(_spec(tmp_path, libraries=[{"name": "Missing"}]), ResourceCache(tmp_path))
    with pytest.raises(ConfigError, match=r"optional.*not supported"):
        load_config({"sources": [{"type": "enrichr_remote", "libraries": ["KEGG"], "optional": True}]})
    with pytest.raises(ConfigError, match=r"optional.*not supported"):
        load_config({"sources": [{"type": "enrichr_remote", "libraries": [{"name": "KEGG", "optional": True}]}]})


def test_empty_library_download_is_rejected(tmp_path: Path, fake_http: FakeServer) -> None:
    _serve(fake_http)
    fake_http.routes[LIBRARY.format(library="KEGG_2021_Human")] = "   \n"
    with pytest.raises(ResourceError, match="was empty"):
        load_enrichr_remote(_spec(tmp_path), ResourceCache(tmp_path))
    assert not (tmp_path / "enrichr" / "KEGG_2021_Human.gmt").exists()


def test_single_library_shorthands() -> None:
    assert EnrichrRemoteSource.model_validate({"type": "enrichr_remote", "library": "X"}).libraries[0].name == "X"
    assert EnrichrRemoteSource.model_validate({"type": "enrichr_remote", "match": "^X"}).libraries[0].match == "^X"


def test_enrichr_build_smoke(tmp_path: Path, fake_http: FakeServer, mini_world: Any) -> None:
    _serve(fake_http)
    cfg = mini_world.config(mini_world.gmt("unused.gmt", {"T": ["G1"]}))
    cfg["sources"] = [_spec(tmp_path).model_dump(mode="json", exclude_none=True)]
    cfg["filtering"] = {"min_size_default": 1}
    result = build(config=cfg)
    manifest = pd.read_csv(result.outdir / "metadata" / "source_manifest.tsv", sep="\t", dtype=str).fillna("")
    assert manifest.loc[0, "source"] == "Enrichr"
    assert manifest.loc[0, "candidate_libraries"] == '["KEGG_2021_Human"]'
    assert "ENRICHR_KEGG_2021_Human__TERM_A" in result.gmt_paths["symbols"].read_text(encoding="utf-8")
