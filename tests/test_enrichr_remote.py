import json

import pandas as pd
import pytest

from alien import build
from alien import sources
from alien.reports import build_source_manifest
from alien.sources import read_enrichr_remote, resolve_enrichr_library


class FakeResponse:
    def __init__(self, text, status_code=200):
        self.text = text
        self.status_code = status_code

    def raise_for_status(self):
        if self.status_code >= 400:
            raise RuntimeError(f"HTTP {self.status_code}")


def test_enrichr_remote_exact_name_downloads_and_normalizes(tmp_path, monkeypatch):
    calls = []
    monkeypatch.setattr(sources.requests, "get", _fake_enrichr_get(calls))

    frame = read_enrichr_remote(_enrichr_spec(tmp_path), tmp_path)

    assert calls == [
        "https://example.test/enrichr/datasetStatistics",
        "https://example.test/enrichr/geneSetLibrary?mode=text&libraryName=KEGG_2021_Human",
    ]
    assert frame["term_id"].unique().tolist() == ["ENRICHR_KEGG_2021_Human__TERM_A"]
    assert frame["source"].unique().tolist() == ["Enrichr"]
    assert frame["source_tag"].unique().tolist() == ["KEGG"]
    assert frame["collection"].unique().tolist() == ["KEGG_2021_Human"]
    assert frame["gene_symbol"].tolist() == ["TP53", "GENE2"]
    metadata = json.loads(frame["metadata_json"].iloc[0])
    assert metadata["selected_library"] == "KEGG_2021_Human"
    assert metadata["match_method"] == "exact"


def test_enrichr_remote_reuses_cached_metadata_and_gmt(tmp_path, monkeypatch):
    monkeypatch.setattr(sources.requests, "get", _fake_enrichr_get([]))
    spec = _enrichr_spec(tmp_path)
    first = read_enrichr_remote(spec, tmp_path)

    def fail_get(*args, **kwargs):
        raise AssertionError("cache should avoid network")

    monkeypatch.setattr(sources.requests, "get", fail_get)
    second = read_enrichr_remote(spec, tmp_path)

    assert first.equals(second)


def test_enrichr_remote_force_refreshes_cached_files(tmp_path, monkeypatch):
    monkeypatch.setattr(sources.requests, "get", _fake_enrichr_get([], term_name="TERM_A"))
    spec = _enrichr_spec(tmp_path)
    read_enrichr_remote(spec, tmp_path)

    calls = []
    monkeypatch.setattr(sources.requests, "get", _fake_enrichr_get(calls, term_name="TERM_B"))
    refreshed = read_enrichr_remote(spec, tmp_path, force=True)

    assert calls == [
        "https://example.test/enrichr/datasetStatistics",
        "https://example.test/enrichr/geneSetLibrary?mode=text&libraryName=KEGG_2021_Human",
    ]
    assert refreshed["term_id"].unique().tolist() == ["ENRICHR_KEGG_2021_Human__TERM_B"]


def test_enrichr_remote_regex_match_chooses_latest_year_then_shortest_name():
    names = ["ClinVar_2019", "ClinVar_2022", "ClinVar_Long_2022", "ClinVar_2021"]

    selected, candidates, method = resolve_enrichr_library(names, {"name": "ClinVar", "match": r"^ClinVar.*[0-9]{4}$"})

    assert selected == "ClinVar_2022"
    assert candidates == ["ClinVar_2019", "ClinVar_2021", "ClinVar_2022", "ClinVar_Long_2022"]
    assert method == "regex"


def test_source_manifest_records_enrichr_regex_resolution(tmp_path, monkeypatch):
    monkeypatch.setattr(sources.requests, "get", _fake_enrichr_get([]))
    frame = read_enrichr_remote(
        _enrichr_spec(
            tmp_path,
            libraries=[
                {
                    "name": "ClinVar",
                    "match": r"^ClinVar_[0-9]{4}$",
                    "source_tag": "CLINVAR",
                    "family": "disease_phenotype",
                    "aspect": "disease",
                }
            ],
        ),
        tmp_path,
    )

    manifest = build_source_manifest(frame)

    assert manifest.loc[0, "collection"] == "ClinVar_2022"
    assert manifest.loc[0, "selected_library"] == "ClinVar_2022"
    assert manifest.loc[0, "configured_name"] == "ClinVar"
    assert manifest.loc[0, "match"] == r"^ClinVar_[0-9]{4}$"
    assert manifest.loc[0, "match_method"] == "regex"
    assert manifest.loc[0, "candidate_libraries"] == '["ClinVar_2019", "ClinVar_2022"]'


def test_enrichr_remote_unresolved_library_errors(tmp_path, monkeypatch):
    monkeypatch.setattr(sources.requests, "get", _fake_enrichr_get([]))

    with pytest.raises(ValueError, match="was not found"):
        read_enrichr_remote(_enrichr_spec(tmp_path, libraries=[{"name": "Missing"}]), tmp_path)


def test_optional_source_key_is_not_supported(tmp_path):
    with pytest.raises(ValueError, match="optional.*not supported"):
        read_enrichr_remote(_enrichr_spec(tmp_path, optional=True), tmp_path)

    with pytest.raises(ValueError, match="optional.*not supported"):
        read_enrichr_remote(_enrichr_spec(tmp_path, libraries=[{"name": "KEGG_2021_Human", "optional": True}]), tmp_path)


def test_enrichr_remote_build_smoke(tmp_path, monkeypatch):
    source_dir = tmp_path / "sources"
    outdir = tmp_path / "out"
    hgnc_dir = source_dir / "hgnc"
    hgnc_dir.mkdir(parents=True)
    pd.DataFrame(
        {
            "symbol": ["TP53", "GENE2"],
            "prev_symbol": ["", ""],
            "alias_symbol": ["", ""],
            "ensembl_gene_id": ["ENSG00000141510", "ENSG000002"],
        }
    ).to_csv(hgnc_dir / "hgnc_complete_set.txt", sep="\t", index=False)
    gtf = tmp_path / "genes.gtf"
    gtf.write_text(
        'chr1\tALIEN\tgene\t1\t10\t.\t+\t.\tgene_id "ENSG00000141510.18"; gene_name "TP53"; gene_type "protein_coding";\n'
        'chr1\tALIEN\tgene\t20\t30\t.\t+\t.\tgene_id "ENSG000002.1"; gene_name "GENE2"; gene_type "protein_coding";\n',
        encoding="utf-8",
    )
    monkeypatch.setattr(sources.requests, "get", _fake_enrichr_get([]))

    result = build(
        {
            "project": {"source_dir": str(source_dir)},
            "sources": [
                {
                    "type": "enrichr_remote",
                    "cache_dir": str(source_dir / "enrichr"),
                    "dataset_statistics_endpoint": "https://example.test/enrichr/datasetStatistics",
                    "library_download_endpoint": "https://example.test/enrichr/geneSetLibrary?mode=text&libraryName={library}",
                    "libraries": [
                        {
                            "name": "KEGG_2021_Human",
                            "source_tag": "KEGG",
                            "family": "biology_process_pathway",
                            "aspect": "pathway",
                        }
                    ],
                }
            ],
            "targets": [
                {
                    "name": "human_test",
                    "type": "ensembl_gtf",
                    "annotation": {"source": "Test", "version": "test", "path": str(gtf)},
                }
            ],
            "filtering": {
                "min_size_default": 1,
                "max_size_default": 10,
                "max_size_disease": 10,
                "max_size_cancer": 10,
            },
            "ncbi_gene": {"enabled": False},
            "ensembl_archive": {"enabled": False},
        },
        outdir=outdir,
        workers=1,
    )

    assert result.namespaces == ("symbols", "human_test")
    assert "ENRICHR_KEGG_2021_Human__TERM_A" in (outdir / "gmt" / "symbols.gmt").read_text(encoding="utf-8")
    assert "ENSG00000141510" in (outdir / "gmt" / "human_test.gmt").read_text(encoding="utf-8")
    source_manifest = pd.read_csv(outdir / "metadata" / "source_manifest.tsv", sep="\t", dtype=str).fillna("")
    assert source_manifest.loc[0, "source"] == "Enrichr"
    assert source_manifest.loc[0, "collection"] == "KEGG_2021_Human"
    assert source_manifest.loc[0, "selected_library"] == "KEGG_2021_Human"
    assert source_manifest.loc[0, "match_method"] == "exact"
    assert source_manifest.loc[0, "candidate_libraries"] == '["KEGG_2021_Human"]'


def _enrichr_spec(tmp_path, **overrides):
    spec = {
        "cache_dir": str(tmp_path / "enrichr"),
        "dataset_statistics_endpoint": "https://example.test/enrichr/datasetStatistics",
        "library_download_endpoint": "https://example.test/enrichr/geneSetLibrary?mode=text&libraryName={library}",
        "libraries": [
            {
                "name": "KEGG_2021_Human",
                "source_tag": "KEGG",
                "family": "biology_process_pathway",
                "aspect": "pathway",
            }
        ],
    }
    spec.update(overrides)
    return spec


def _fake_enrichr_get(calls, term_name="TERM_A"):
    def fake_get(url, timeout=None):
        calls.append(url)
        if url == "https://example.test/enrichr/datasetStatistics":
            return FakeResponse(
                json.dumps(
                    {
                        "statistics": [
                            {"libraryName": "KEGG_2021_Human", "numTerms": 1},
                            {"libraryName": "ClinVar_2019", "numTerms": 1},
                            {"libraryName": "ClinVar_2022", "numTerms": 1},
                        ]
                    }
                )
            )
        if url == "https://example.test/enrichr/geneSetLibrary?mode=text&libraryName=KEGG_2021_Human":
            return FakeResponse(f"{term_name}\tdesc\tTP53\tGENE2\n")
        if url == "https://example.test/enrichr/geneSetLibrary?mode=text&libraryName=ClinVar_2022":
            return FakeResponse(f"{term_name}\tdesc\tTP53\tGENE2\n")
        raise AssertionError(f"Unexpected URL: {url}")

    return fake_get
