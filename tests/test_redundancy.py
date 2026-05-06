from alien.filtering import remove_redundancy


def _record(source, genes, name="name"):
    return {
        "genes": set(genes),
        "meta": {
            "family": "biology_process_pathway",
            "source_tag": source,
            "display_name": name,
            "original_name": name,
        },
    }


def _cfg(cutoff=0.85):
    return {
        "redundancy": {"jaccard_cutoff": cutoff},
        "source_priority": {"biology_process_pathway": ["REACTOME", "GOBP"]},
    }


def test_exact_duplicate_removed_and_source_priority_kept():
    terms = {
        "A": _record("GOBP", ["1", "2"]),
        "B": _record("REACTOME", ["1", "2"]),
    }
    kept, removed, summary = remove_redundancy(terms, "symbols", _cfg())
    assert list(kept) == ["B"]
    assert removed.iloc[0]["reason"] == "exact_duplicate"
    assert summary.iloc[0]["n_removed_exact"] == 1


def test_jaccard_cluster_and_cutoff():
    terms = {
        "A": _record("GOBP", [f"G{i}" for i in range(20)]),
        "B": _record("REACTOME", [f"G{i}" for i in range(19)] + ["GX"]),
        "C": _record("REACTOME", list("ABCXYZ")),
    }
    kept, removed, _ = remove_redundancy(terms, "symbols", _cfg(0.85))
    assert "B" in kept
    assert "A" not in kept
    assert "C" in kept
    assert removed[removed["reason"].eq("jaccard_duplicate")].iloc[0]["term_id"] == "A"
