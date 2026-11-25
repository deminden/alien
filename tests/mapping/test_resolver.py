from __future__ import annotations

from typing import Any

import pandas as pd

from alien.config import load_config
from alien.mapping import map_namespaces, resolve_symbol
from alien.mapping.resolver import choose_candidate
from alien.model import AnnotationCandidate
from alien.resources.ncbi import NcbiMaps, NcbiRecord

from conftest import hgnc_maps, make_target, membership

CONFIG = {"sources": [{"type": "symbol_gmt", "path": "unused.gmt"}], "outputs": {"include_symbols": False}}


class FakeArchive:
    def __init__(self) -> None:
        self.prefetched: list[str] = []
        self.lookups: list[str] = []

    def prefetch(self, stable_ids: list[str], namespace: str = "") -> None:
        self.prefetched += stable_ids

    def lookup(self, stable_id: str) -> dict[str, Any]:
        self.lookups.append(stable_id)
        replacements = ["ENSG00000111111"] if stable_id == "ENSG00000199999" else []
        return {
            "query_id": stable_id,
            "found": True,
            "archive_id": stable_id,
            "is_current": "0",
            "possible_replacement": replacements,
        }


def _map(
    rows: list[dict[str, str]], target: Any, hgnc: Any, ncbi: NcbiMaps | None = None, resolver: Any = None, **cfg: Any
):
    config = load_config({**CONFIG, **cfg})
    return map_namespaces(pd.DataFrame(rows), config, [target], hgnc=hgnc, ncbi=ncbi or NcbiMaps(), resolver=resolver)


def _statuses(result: Any) -> dict[str, int]:
    return dict(zip(result.status_summary["mapping_status"], result.status_summary["n_gene_members"], strict=True))


def test_symbol_resolution_and_ambiguity() -> None:
    hgnc = hgnc_maps([("TP53", "OLD1", "", "ENSG1"), ("A", "", "AMB", "ENSG2"), ("B", "", "AMB", "ENSG3")])
    assert resolve_symbol("TP53", hgnc).status == "exact_current_symbol"
    assert (resolve_symbol("OLD1", hgnc).symbol, resolve_symbol("OLD1", hgnc).status) == (
        "TP53",
        "hgnc_previous_symbol",
    )
    assert resolve_symbol("AMB", hgnc).status == "ambiguous_alias"
    assert resolve_symbol(" ", hgnc).status == "unmapped"


def test_empty_candidate_list_chooses_nothing() -> None:
    assert choose_candidate([]) is None
    coding = AnnotationCandidate("E1", "protein_coding")
    assert choose_candidate([coding, AnnotationCandidate("E2", "lncRNA")]) == coding


def test_ensembl_mapping_priority_and_output_gene_projection() -> None:
    hgnc = hgnc_maps(
        [
            ("TP53", "OLD1", "", "ENSG00000141510"),
            ("SRC", "", "", "ENSG00000999999"),
            ("A", "", "AMB", "ENSGA"),
            ("B", "", "AMB", "ENSGB"),
        ]
    )
    target = make_target(
        "human_test",
        [
            {"gene_symbol": "TP53", "ensembl_gene_id": "ENSG00000141510"},
            {"gene_symbol": "SRC", "ensembl_gene_id": "ENSG00000999999"},
            {"gene_symbol": "RESCUED", "ensembl_gene_id": "ENSG00000111111"},
        ],
        output_ids={"ENSG00000141510", "ENSG00000999999", "ENSG00000111111", "ENSG00000299999"},
    )
    rows = [
        membership("T1", "TP53"),
        membership("T2", "OLD1"),
        membership("T3", "AMB"),
        membership("T4", "SRC", "ENSG00000999999.7"),
        membership("T5", "MISSING"),
        membership("T6", "SRC", "ENSG00000STALE.1"),
        membership("T7", "UNKNOWN", "ENSG00000199999.1"),
        membership("T8", "MATRIXONLY", "ENSG00000299999.1"),
        membership("T9", "1-MAR"),
    ]
    result = _map(rows, target, hgnc, resolver=FakeArchive())
    terms = result.terms["human_test"]
    assert terms["T1"].genes == {"ENSG00000141510"}
    assert terms["T2"].genes == {"ENSG00000141510"}
    assert terms["T4"].genes == {"ENSG00000999999"}
    assert terms["T6"].genes == {"ENSG00000999999"}
    assert terms["T7"].genes == {"ENSG00000111111"}
    assert terms["T8"].genes == {"ENSG00000299999"}
    assert not {"T3", "T5", "T9"} & set(terms)
    assert set(result.ambiguous["input_gene"]) == {"AMB"}
    assert set(result.unmapped["input_gene"]) == {"MISSING", "1-MAR"}
    assert result.archive_audit["rescue_status"].tolist() == ["hgnc_current_symbol", "archive_ensembl_rescue"]
    assert _statuses(result) == {
        "direct_source_ensembl": 2,
        "hgnc_current_symbol": 2,
        "hgnc_previous_symbol": 1,
        "archive_ensembl_rescue": 1,
        "ambiguous_drop": 1,
        "unmapped_drop": 2,
    }


def test_archive_ignores_current_annotation_ids_outside_restricted_output_genes() -> None:
    hgnc = hgnc_maps([("OUTSIDE", "", "", "ENSG00000000001")])
    target = make_target(
        "human_test",
        [
            {"gene_symbol": "OUTSIDE", "ensembl_gene_id": "ENSG00000000001", "is_in_output_genes": False},
            {"gene_symbol": "INSIDE", "ensembl_gene_id": "ENSG00000000002"},
        ],
        output_ids={"ENSG00000000002"},
    )
    archive = FakeArchive()
    result = _map([membership("T", "OUTSIDE", "ENSG00000000001.1")], target, hgnc, resolver=archive)
    assert result.terms["human_test"] == {}
    assert result.unmapped["reason"].tolist() == ["not_in_output_genes"]
    assert archive.prefetched == [] and archive.lookups == []
    assert result.archive_audit.empty


def test_archive_ignores_current_hgnc_ids_absent_from_old_annotation() -> None:
    hgnc = hgnc_maps([("CURRENTONLY", "", "", "ENSG00000000003")])
    target = make_target("old", [{"gene_symbol": "INSIDE", "ensembl_gene_id": "ENSG00000000002"}], {"ENSG00000000002"})
    archive = FakeArchive()
    result = _map([membership("T", "CURRENTONLY", "ENSG00000000003.1")], target, hgnc, resolver=archive)
    assert result.unmapped["reason"].tolist() == ["not_in_output_genes"]
    assert archive.lookups == []


def test_multiple_targets() -> None:
    hgnc = hgnc_maps([("TP53", "", "", "ENSG00000141510"), ("GENE2", "", "", "ENSG000002")])
    rows = [membership("T", "TP53"), membership("T", "GENE2")]
    target_a = make_target(
        "a",
        [
            {"gene_symbol": "TP53", "ensembl_gene_id": "ENSG00000141510"},
            {"gene_symbol": "GENE2", "ensembl_gene_id": "ENSG000002"},
        ],
    )
    target_b = make_target("b", [{"gene_symbol": "TP53", "ensembl_gene_id": "ENSG00000141510"}], {"ENSG00000141510"})
    config = load_config({**CONFIG, "ensembl_archive": {"enabled": False}})
    result = map_namespaces(pd.DataFrame(rows), config, [target_a, target_b], hgnc=hgnc, ncbi=NcbiMaps(), resolver=None)
    assert result.terms["a"]["T"].genes == {"ENSG00000141510", "ENSG000002"}
    assert result.terms["b"]["T"].genes == {"ENSG00000141510"}
    assert result.unmapped["target_namespace"].tolist() == ["b"]
    assert set(result.status_summary["target_namespace"]) == {"a", "b"}


def test_non_gene_tokens_are_dropped_before_mapping() -> None:
    hgnc = hgnc_maps([("TP53", "", "", "ENSG1")])
    target = make_target("t", [{"gene_symbol": "TP53", "ensembl_gene_id": "ENSG1"}])
    result = _map([membership("T", "disease")], target, hgnc)
    assert result.terms["t"] == {}
    assert result.unmapped.empty and result.ambiguous.empty
    assert result.non_gene["reason"].tolist() == ["configured_non_gene_token"]


def test_hgnc_tiebreak_resolves_duplicated_annotation_symbol() -> None:
    hgnc = hgnc_maps([("COX2", "", "", "ENSG_COX2_GOOD")])
    target = make_target(
        "t",
        [
            {"gene_symbol": "COX2", "ensembl_gene_id": "ENSG_COX2_BAD"},
            {"gene_symbol": "COX2", "ensembl_gene_id": "ENSG_COX2_GOOD"},
        ],
    )
    result = _map([membership("T", "COX2")], target, hgnc)
    assert result.terms["t"]["T"].genes == {"ENSG_COX2_GOOD"}
    assert _statuses(result) == {"hgnc_ensembl_tiebreak": 1}


NCBI = NcbiMaps(stages={"ncbi_gene_history": {"OLDNCBI": NcbiRecord("NEWGENE", ("123",), ("ENSG_NEW",))}}, ambiguous={})


def test_ncbi_history_rescues_unmapped_legacy_symbol() -> None:
    hgnc = hgnc_maps([("NEWGENE", "", "", "ENSG_NEW")])
    target = make_target("t", [{"gene_symbol": "NEWGENE", "ensembl_gene_id": "ENSG_NEW"}])
    result = _map([membership("T", "OLDNCBI")], target, hgnc, NCBI, ncbi_gene={"allowed_source_tags": ["REACTOME"]})
    assert result.terms["t"]["T"].genes == {"ENSG_NEW"}
    assert result.ncbi_audit["action"].tolist() == ["rescued"]
    assert _statuses(result) == {"ncbi_gene_history": 1}


def test_ncbi_rescue_can_be_restricted_by_source_tag() -> None:
    hgnc = hgnc_maps([("NEWGENE", "", "", "ENSG_NEW")])
    target = make_target("t", [{"gene_symbol": "NEWGENE", "ensembl_gene_id": "ENSG_NEW"}])
    result = _map([membership("T", "OLDNCBI")], target, hgnc, NCBI, ncbi_gene={"allowed_source_tags": ["CCLE"]})
    assert result.terms["t"] == {}
    assert result.unmapped["reason"].tolist() == ["unmapped_symbol"]
    assert result.ncbi_audit.empty
