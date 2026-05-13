import pandas as pd

from alien.mapping import EnsemblArchiveResolver, map_ensembl, map_namespaces, resolve_current_symbol
from alien.sources import build_hgnc_maps


def _term_df():
    return pd.DataFrame(
        [
            _row("T1", "TP53", ""),
            _row("T2", "OLD1", ""),
            _row("T3", "AMB", ""),
            _row("T4", "SRC", "ENSG00000999999.7"),
            _row("T5", "MISSING", ""),
            _row("T6", "SRC", "ENSG00000STALE.1"),
            _row("T7", "UNKNOWN", "ENSG00000199999.1"),
            _row("T8", "MATRIXONLY", "ENSG00000299999.1"),
            _row("T9", "1-MAR", ""),
        ]
    )


def _row(term_id, gene, ensembl):
    return {
        "term_id": term_id,
        "original_name": term_id,
        "display_name": term_id,
        "description": "",
        "source": "test",
        "source_tag": "REACTOME",
        "collection": "C2",
        "subcollection": "",
        "family": "biology_process_pathway",
        "aspect": "pathway",
        "db_version": "",
        "source_url": "",
        "source_license_note": "",
        "gene_symbol": gene,
        "gene_ensembl_from_source": ensembl,
        "metadata_json": "{}",
    }


def test_symbol_resolution_and_ambiguity():
    hgnc = build_hgnc_maps(
        pd.DataFrame(
            [
                {"symbol": "TP53", "prev_symbol": "OLD1", "alias_symbol": "", "ensembl_gene_id": "ENSG1"},
                {"symbol": "A", "prev_symbol": "", "alias_symbol": "AMB", "ensembl_gene_id": "ENSG2"},
                {"symbol": "B", "prev_symbol": "", "alias_symbol": "AMB", "ensembl_gene_id": "ENSG3"},
            ]
        )
    )
    assert resolve_current_symbol("TP53", hgnc)[1] == "exact_current_symbol"
    assert resolve_current_symbol("OLD1", hgnc)[:2] == ("TP53", "hgnc_previous_symbol")
    assert resolve_current_symbol("AMB", hgnc)[1] == "ambiguous_alias"


def test_ensembl_mapping_priority_and_universe_projection():
    hgnc = build_hgnc_maps(
        pd.DataFrame(
            [
                {"symbol": "TP53", "prev_symbol": "OLD1", "alias_symbol": "", "ensembl_gene_id": "ENSG00000141510"},
                {"symbol": "SRC", "prev_symbol": "", "alias_symbol": "", "ensembl_gene_id": "ENSG00000999999"},
                {"symbol": "A", "prev_symbol": "", "alias_symbol": "AMB", "ensembl_gene_id": "ENSGA"},
                {"symbol": "B", "prev_symbol": "", "alias_symbol": "AMB", "ensembl_gene_id": "ENSGB"},
            ]
        )
    )
    gencode = pd.DataFrame(
        [
            {"gene_symbol": "TP53", "ensembl_gene_id": "ENSG00000141510", "gene_biotype": "protein_coding", "is_in_expression_universe": True},
            {"gene_symbol": "SRC", "ensembl_gene_id": "ENSG00000999999", "gene_biotype": "protein_coding", "is_in_expression_universe": True},
            {"gene_symbol": "RESCUED", "ensembl_gene_id": "ENSG00000111111", "gene_biotype": "protein_coding", "is_in_expression_universe": True},
        ]
    )
    unmapped = []
    ambiguous = []
    non_gene = []
    archive_audit = []
    ncbi_audit = []
    mapping_status = {}
    mapped = map_ensembl(
        _term_df(),
        hgnc,
        {},
        gencode,
        {"ENSG00000141510", "ENSG00000999999", "ENSG00000111111", "ENSG00000299999"},
        "human_test",
        {"gene_mapping": {"manual_symbol_repairs": {}}},
        unmapped,
        ambiguous,
        non_gene,
        archive_audit,
        ncbi_audit,
        mapping_status,
        _FakeArchiveResolver(),
    )
    assert mapped["T1"]["genes"] == {"ENSG00000141510"}
    assert mapped["T2"]["genes"] == {"ENSG00000141510"}
    assert mapped["T4"]["genes"] == {"ENSG00000999999"}
    assert mapped["T6"]["genes"] == {"ENSG00000999999"}
    assert mapped["T7"]["genes"] == {"ENSG00000111111"}
    assert mapped["T8"]["genes"] == {"ENSG00000299999"}
    assert "T3" not in mapped
    assert "T5" not in mapped
    assert "T9" not in mapped
    assert {row["input_gene"] for row in ambiguous} == {"AMB"}
    assert {row["input_gene"] for row in unmapped} == {"MISSING", "1-MAR"}
    assert [row["rescue_status"] for row in archive_audit] == ["hgnc_current_symbol", "archive_ensembl_rescue"]
    assert mapping_status[("human_test", "biology_process_pathway", "REACTOME", "direct_source_ensembl")] == 2
    assert mapping_status[("human_test", "biology_process_pathway", "REACTOME", "hgnc_current_symbol")] == 2
    assert mapping_status[("human_test", "biology_process_pathway", "REACTOME", "hgnc_previous_symbol")] == 1
    assert mapping_status[("human_test", "biology_process_pathway", "REACTOME", "archive_ensembl_rescue")] == 1
    assert mapping_status[("human_test", "biology_process_pathway", "REACTOME", "ambiguous_drop")] == 1
    assert mapping_status[("human_test", "biology_process_pathway", "REACTOME", "unmapped_drop")] == 2


def test_archive_lookup_ignores_current_annotation_ids_outside_restricted_universe():
    hgnc = build_hgnc_maps(
        pd.DataFrame([{"symbol": "OUTSIDE", "prev_symbol": "", "alias_symbol": "", "ensembl_gene_id": "ENSG00000000001"}])
    )
    gencode = pd.DataFrame(
        [
            {"gene_symbol": "OUTSIDE", "ensembl_gene_id": "ENSG00000000001", "gene_biotype": "protein_coding", "is_in_expression_universe": False},
            {"gene_symbol": "INSIDE", "ensembl_gene_id": "ENSG00000000002", "gene_biotype": "protein_coding", "is_in_expression_universe": True},
        ]
    )
    unmapped = []
    archive_audit = []
    resolver = _FakeArchiveResolver()

    mapped = map_ensembl(
        pd.DataFrame([_row("T_OUTSIDE", "OUTSIDE", "ENSG00000000001.1")]),
        hgnc,
        {},
        gencode,
        {"ENSG00000000002"},
        "human_test",
        {"gene_mapping": {"manual_symbol_repairs": {}}},
        unmapped,
        [],
        [],
        archive_audit=archive_audit,
        archive_resolver=resolver,
    )

    assert mapped == {}
    assert unmapped[0]["reason"] == "not_in_target_universe"
    assert resolver.stable_ids == []
    assert resolver.lookups == []
    assert archive_audit == []


def test_archive_lookup_ignores_current_hgnc_ids_absent_from_old_annotation():
    hgnc = build_hgnc_maps(
        pd.DataFrame([{"symbol": "CURRENTONLY", "prev_symbol": "", "alias_symbol": "", "ensembl_gene_id": "ENSG00000000003"}])
    )
    gencode = pd.DataFrame(
        [{"gene_symbol": "INSIDE", "ensembl_gene_id": "ENSG00000000002", "gene_biotype": "protein_coding", "is_in_expression_universe": True}]
    )
    unmapped = []
    archive_audit = []
    resolver = _FakeArchiveResolver()

    mapped = map_ensembl(
        pd.DataFrame([_row("T_CURRENTONLY", "CURRENTONLY", "ENSG00000000003.1")]),
        hgnc,
        {},
        gencode,
        {"ENSG00000000002"},
        "old_annotation",
        {"gene_mapping": {"manual_symbol_repairs": {}}},
        unmapped,
        [],
        [],
        archive_audit=archive_audit,
        archive_resolver=resolver,
    )

    assert mapped == {}
    assert unmapped[0]["reason"] == "not_in_target_universe"
    assert resolver.stable_ids == []
    assert resolver.lookups == []
    assert archive_audit == []


def test_map_namespaces_parallel_maps_multiple_targets():
    hgnc = build_hgnc_maps(
        pd.DataFrame(
            [
                {"symbol": "TP53", "prev_symbol": "", "alias_symbol": "", "ensembl_gene_id": "ENSG00000141510"},
                {"symbol": "GENE2", "prev_symbol": "", "alias_symbol": "", "ensembl_gene_id": "ENSG000002"},
            ]
        )
    )
    term_df = pd.DataFrame([_row("T_PARALLEL", "TP53", ""), _row("T_PARALLEL", "GENE2", "")])
    targets = [
        {
            "name": "target_a",
            "annotation": pd.DataFrame(
                [
                    {"gene_symbol": "TP53", "ensembl_gene_id": "ENSG00000141510", "gene_biotype": "protein_coding", "is_in_expression_universe": True},
                    {"gene_symbol": "GENE2", "ensembl_gene_id": "ENSG000002", "gene_biotype": "protein_coding", "is_in_expression_universe": True},
                ]
            ),
            "universe": {"ENSG00000141510", "ENSG000002"},
        },
        {
            "name": "target_b",
            "annotation": pd.DataFrame(
                [{"gene_symbol": "TP53", "ensembl_gene_id": "ENSG00000141510", "gene_biotype": "protein_coding", "is_in_expression_universe": True}]
            ),
            "universe": {"ENSG00000141510"},
        },
    ]

    mapped, unmapped, ambiguous, non_gene, archive_audit, ncbi_audit, mapping_status = map_namespaces(
        term_df,
        hgnc,
        {},
        targets,
        {"outputs": {"include_symbols": False}, "gene_mapping": {"manual_symbol_repairs": {}}, "ensembl_archive": {"enabled": False}},
        workers=2,
    )

    assert mapped["target_a"]["T_PARALLEL"]["genes"] == {"ENSG00000141510", "ENSG000002"}
    assert mapped["target_b"]["T_PARALLEL"]["genes"] == {"ENSG00000141510"}
    assert unmapped["target_namespace"].tolist() == ["target_b"]
    assert ambiguous.empty
    assert non_gene.empty
    assert archive_audit.empty
    assert ncbi_audit.empty
    assert set(mapping_status["target_namespace"]) == {"target_a", "target_b"}


def test_configured_non_gene_tokens_are_dropped_before_mapping():
    hgnc = build_hgnc_maps(pd.DataFrame([{"symbol": "TP53", "prev_symbol": "", "alias_symbol": "", "ensembl_gene_id": "ENSG1"}]))
    unmapped = []
    ambiguous = []
    non_gene = []
    mapped = map_ensembl(
        pd.DataFrame([_row("T_DIRTY", "disease", "")]),
        hgnc,
        {},
        pd.DataFrame([{"gene_symbol": "TP53", "ensembl_gene_id": "ENSG1", "gene_biotype": "protein_coding", "is_in_expression_universe": True}]),
        {"ENSG1"},
        "human_test",
        {"gene_mapping": {"non_gene_tokens": ["disease"], "manual_symbol_repairs": {}}},
        unmapped,
        ambiguous,
        non_gene,
    )
    assert mapped == {}
    assert unmapped == []
    assert ambiguous == []
    assert non_gene[0]["reason"] == "configured_non_gene_token"


def test_hgnc_ensembl_tiebreak_resolves_ambiguous_gencode_symbol():
    hgnc = build_hgnc_maps(pd.DataFrame([{"symbol": "COX2", "prev_symbol": "", "alias_symbol": "", "ensembl_gene_id": "ENSG_COX2_GOOD"}]))
    gencode = pd.DataFrame(
        [
            {"gene_symbol": "COX2", "ensembl_gene_id": "ENSG_COX2_BAD", "gene_biotype": "protein_coding", "is_in_expression_universe": True},
            {"gene_symbol": "COX2", "ensembl_gene_id": "ENSG_COX2_GOOD", "gene_biotype": "protein_coding", "is_in_expression_universe": True},
        ]
    )
    mapping_status = {}
    mapped = map_ensembl(
        pd.DataFrame([_row("T_COX2", "COX2", "")]),
        hgnc,
        {},
        gencode,
        {"ENSG_COX2_BAD", "ENSG_COX2_GOOD"},
        "human_test",
        {"gene_mapping": {"manual_symbol_repairs": {}}},
        [],
        [],
        [],
        mapping_status=mapping_status,
    )
    assert mapped["T_COX2"]["genes"] == {"ENSG_COX2_GOOD"}
    assert mapping_status[("human_test", "biology_process_pathway", "REACTOME", "hgnc_ensembl_tiebreak")] == 1


def test_ncbi_history_rescues_unmapped_legacy_symbol():
    hgnc = build_hgnc_maps(pd.DataFrame([{"symbol": "NEWGENE", "prev_symbol": "", "alias_symbol": "", "ensembl_gene_id": "ENSG_NEW"}]))
    ncbi = {
        "stages": {
            "ncbi_gene_history": {
                "OLDNCBI": {
                    "current_symbol": "NEWGENE",
                    "current_gene_ids": ["123"],
                    "ensembl_gene_ids": ["ENSG_NEW"],
                }
            }
        },
        "ambiguous": {},
    }
    ncbi_audit = []
    mapping_status = {}
    mapped = map_ensembl(
        pd.DataFrame([_row("T_NCBI", "OLDNCBI", "")]),
        hgnc,
        ncbi,
        pd.DataFrame([{"gene_symbol": "NEWGENE", "ensembl_gene_id": "ENSG_NEW", "gene_biotype": "protein_coding", "is_in_expression_universe": True}]),
        {"ENSG_NEW"},
        "human_test",
        {"gene_mapping": {"manual_symbol_repairs": {}}, "ncbi_gene": {"allowed_source_tags": ["REACTOME"]}},
        [],
        [],
        [],
        ncbi_audit=ncbi_audit,
        mapping_status=mapping_status,
    )
    assert mapped["T_NCBI"]["genes"] == {"ENSG_NEW"}
    assert ncbi_audit[0]["action"] == "rescued"
    assert mapping_status[("human_test", "biology_process_pathway", "REACTOME", "ncbi_gene_history")] == 1


def test_ncbi_rescue_can_be_restricted_by_source_tag():
    hgnc = build_hgnc_maps(pd.DataFrame([{"symbol": "NEWGENE", "prev_symbol": "", "alias_symbol": "", "ensembl_gene_id": "ENSG_NEW"}]))
    ncbi = {
        "stages": {
            "ncbi_gene_history": {
                "OLDNCBI": {
                    "current_symbol": "NEWGENE",
                    "current_gene_ids": ["123"],
                    "ensembl_gene_ids": ["ENSG_NEW"],
                }
            }
        },
        "ambiguous": {},
    }
    unmapped = []
    mapped = map_ensembl(
        pd.DataFrame([_row("T_NCBI_BLOCKED", "OLDNCBI", "")]),
        hgnc,
        ncbi,
        pd.DataFrame([{"gene_symbol": "NEWGENE", "ensembl_gene_id": "ENSG_NEW", "gene_biotype": "protein_coding", "is_in_expression_universe": True}]),
        {"ENSG_NEW"},
        "human_test",
        {"gene_mapping": {"manual_symbol_repairs": {}}, "ncbi_gene": {"allowed_source_tags": ["CCLE"]}},
        unmapped,
        [],
        [],
    )
    assert mapped == {}
    assert unmapped[0]["reason"] == "unmapped_symbol"


class _FakeArchiveResolver:
    def __init__(self):
        self.stable_ids = []
        self.namespace = ""
        self.lookups = []

    def prefetch(self, stable_ids, namespace=""):
        self.stable_ids = stable_ids
        self.namespace = namespace

    def lookup(self, stable_id):
        self.lookups.append(stable_id)
        if stable_id == "ENSG00000199999":
            return {
                "found": True,
                "archive_id": stable_id,
                "is_current": "0",
                "release": "75",
                "assembly": "GRCh38",
                "latest": "",
                "possible_replacement": ["ENSG00000111111"],
            }
        return {"found": True, "archive_id": stable_id, "is_current": "0", "possible_replacement": []}


def test_archive_resolver_retries_cached_errors(tmp_path):
    resolver = _RecordingArchiveResolver(tmp_path / "cache.json")
    resolver.cache = {
        "ENSG_SUCCESS": {"query_id": "ENSG_SUCCESS", "found": True},
        "ENSG_ERROR": {"query_id": "ENSG_ERROR", "found": False, "error": "timeout"},
    }

    resolver.prefetch(["ENSG_SUCCESS", "ENSG_ERROR", "ENSG_NEW"])

    assert resolver.calls == [["ENSG_ERROR", "ENSG_NEW"]]
    assert resolver.cache["ENSG_SUCCESS"] == {"query_id": "ENSG_SUCCESS", "found": True}
    assert resolver.cache["ENSG_ERROR"]["found"] is True
    assert resolver.cache["ENSG_NEW"]["found"] is True


def test_archive_resolver_splits_failed_batches(tmp_path):
    resolver = _RecordingArchiveResolver(tmp_path / "cache.json")
    resolver.fail_multi_id_batches = True

    resolver.prefetch(["ENSG_A", "ENSG_B"])

    assert resolver.calls == [["ENSG_A", "ENSG_B"], ["ENSG_A"], ["ENSG_B"]]
    assert resolver.cache["ENSG_A"]["found"] is True
    assert resolver.cache["ENSG_B"]["found"] is True


class _RecordingArchiveResolver(EnsemblArchiveResolver):
    def __init__(self, cache_path):
        super().__init__(
            cache_path=cache_path,
            endpoint="https://example.test/archive/id",
            batch_size=50,
            timeout_seconds=120,
            connect_timeout_seconds=30,
            retry_attempts=1,
            retry_backoff_seconds=0,
        )
        self.calls = []
        self.fail_multi_id_batches = False

    def _post_archive_batch(self, stable_ids):
        self.calls.append(list(stable_ids))
        if self.fail_multi_id_batches and len(stable_ids) > 1:
            raise TimeoutError("simulated timeout")
        return [{"id": stable_id, "is_current": "1"} for stable_id in stable_ids]
