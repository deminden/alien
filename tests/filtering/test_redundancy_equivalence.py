"""Property tests: the prefix-filter join finds exactly the brute-force pairs, and redundancy removal matches 0.1.x."""

from __future__ import annotations

import itertools
import sys
from pathlib import Path
from typing import Any

from hypothesis import given, settings
from hypothesis import strategies as st

from alien.config.schema import RedundancyConfig
from alien.filtering.redundancy import remove_redundancy, similar_pairs
from alien.model import TermRecord

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "legacy"))
import legacy_filtering

CUTOFFS = st.sampled_from([0.0, 0.25, 0.5, 0.6, 2 / 3, 0.7, 0.75, 0.8, 0.85, 0.9, 1.0])
GENE_SETS = st.lists(
    st.frozensets(st.integers(0, 14).map(lambda i: f"g{i}"), min_size=1, max_size=10), min_size=0, max_size=30
)


def brute_force(gene_sets: list[frozenset[str]], cutoff: float) -> list[tuple[int, int]]:
    pairs = []
    for i, j in itertools.combinations(range(len(gene_sets)), 2):
        a, b = gene_sets[i], gene_sets[j]
        if min(len(a), len(b)) / max(len(a), len(b)) <= cutoff:
            continue
        inter = len(a & b)
        union = len(a) + len(b) - inter
        if union and inter / union > cutoff:
            pairs.append((i, j))
    return pairs


@settings(max_examples=400, deadline=None)
@given(GENE_SETS, CUTOFFS)
def test_prefix_join_equals_brute_force(gene_sets: list[frozenset[str]], cutoff: float) -> None:
    assert similar_pairs(gene_sets, cutoff) == brute_force(gene_sets, cutoff)


@st.composite
def term_sets(draw: st.DrawFn) -> dict[str, Any]:
    sets = draw(GENE_SETS)
    terms = {}
    for index, genes in enumerate(sets):
        tag = draw(st.sampled_from(["A", "B", "C"]))
        name = draw(st.sampled_from(["regulation of x", "x pathway", "process y", "zz"]))
        family = draw(st.sampled_from(["f1", "f2"]))
        terms[f"T{index:02d}"] = {
            "genes": set(genes),
            "meta": {"family": family, "source_tag": tag, "display_name": name},
        }
    return terms


@settings(max_examples=300, deadline=None)
@given(term_sets(), CUTOFFS, st.booleans(), st.booleans())
def test_redundancy_matches_legacy(terms: dict[str, Any], cutoff: float, exact: bool, per_family: bool) -> None:
    priority = {"f1": ["B", "A"], "f2": ["C"]}
    cfg = {
        "redundancy": {"jaccard_cutoff": cutoff, "exact_duplicate_removal": exact, "apply_per_family": per_family},
        "source_priority": priority,
    }
    old_kept, old_removed, _ = legacy_filtering.remove_redundancy(dict(sorted(terms.items())), "ns", cfg)
    records = {
        t: TermRecord(t, "", {k: str(v) for k, v in r["meta"].items()}, set(r["genes"]))
        for t, r in sorted(terms.items())
    }
    config = RedundancyConfig(jaccard_cutoff=cutoff, exact_duplicate_removal=exact, apply_per_family=per_family)
    families = sorted({r.family for r in records.values()}) if per_family else (["*"] if records else [])
    kept: dict[str, Any] = {}
    removed: list[dict[str, Any]] = []
    for family in families:
        subset = {t: r for t, r in records.items() if family in ("*", r.family)}
        result = remove_redundancy(subset, "ns", family, config, priority)
        kept.update(result.kept)
        removed += result.removed.to_dict("records")
    assert list(kept) == list(old_kept)
    old_rows = [
        {**row, "mapped_size": row.pop("mapped_size_after_output_genes", row.get("mapped_size"))}
        for row in old_removed.to_dict("records")
    ]
    assert removed == old_rows


def test_representative_linkage_never_chains() -> None:
    # A~B (8/9) and B~C (8/10) but A~C is 7/10: single linkage chains all three; representative linkage does not.
    genes = {"A": set("abcdefgh"), "B": set("abcdefghi"), "C": set("bcdefghij")}
    tags = {"A": "P", "B": "X", "C": "X"}
    terms = {
        name: TermRecord(name, "", {"family": "f", "source_tag": tags[name], "display_name": name}, members)
        for name, members in genes.items()
    }
    priority = {"f": ["P"]}
    single = remove_redundancy(terms, "ns", "f", RedundancyConfig(jaccard_cutoff=0.75), priority)
    representative = remove_redundancy(
        terms, "ns", "f", RedundancyConfig(jaccard_cutoff=0.75, linkage="representative"), priority
    )
    assert set(single.kept) == {"A"}
    assert set(representative.kept) == {"A", "C"}
    for row in representative.removed.to_dict("records"):
        removed, kept = terms[row["term_id"]].genes, terms[row["representative_term_id"]].genes
        assert len(removed & kept) / len(removed | kept) > 0.75
