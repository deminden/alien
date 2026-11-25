from __future__ import annotations

from alien.config.schema import RedundancyConfig
from alien.filtering import remove_redundancy
from alien.model import TermRecord

PRIORITY = {"biology_process_pathway": ["REACTOME", "GOBP"]}


def _record(term_id: str, source: str, genes: list[str], name: str = "name") -> TermRecord:
    meta = {"family": "biology_process_pathway", "source_tag": source, "display_name": name, "original_name": name}
    return TermRecord(term_id, "", meta, set(genes))


def test_exact_duplicate_removed_and_source_priority_kept() -> None:
    terms = {"A": _record("A", "GOBP", ["1", "2"]), "B": _record("B", "REACTOME", ["1", "2"])}
    result = remove_redundancy(terms, "symbols", "biology_process_pathway", RedundancyConfig(), PRIORITY)
    assert list(result.kept) == ["B"]
    assert result.removed.iloc[0]["reason"] == "exact_duplicate"
    assert result.summary.iloc[0]["n_removed_exact"] == 1


def test_jaccard_cluster_and_cutoff() -> None:
    terms = {
        "A": _record("A", "GOBP", [f"G{i}" for i in range(20)]),
        "B": _record("B", "REACTOME", [f"G{i}" for i in range(19)] + ["GX"]),
        "C": _record("C", "REACTOME", list("ABCXYZ")),
    }
    result = remove_redundancy(
        terms, "symbols", "biology_process_pathway", RedundancyConfig(jaccard_cutoff=0.85), PRIORITY
    )
    assert set(result.kept) == {"B", "C"}
    assert result.removed[result.removed["reason"].eq("jaccard_duplicate")].iloc[0]["term_id"] == "A"


def test_generic_names_lose_ties() -> None:
    terms = {
        "A": _record("A", "GOBP", ["1", "2"], "regulation of thing"),
        "B": _record("B", "GOBP", ["1", "2"], "thing pathway long name"),
    }
    result = remove_redundancy(terms, "ns", "biology_process_pathway", RedundancyConfig(), PRIORITY)
    assert list(result.kept) == ["B"]
