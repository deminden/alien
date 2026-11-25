"""Property test: the new resolver/projection equals the frozen 0.1.x map_ensembl on random inputs."""

from __future__ import annotations

import sys
from pathlib import Path
from typing import Any

import pandas as pd
from hypothesis import HealthCheck, given, settings
from hypothesis import strategies as st

from alien.config import load_config
from alien.mapping import map_namespaces
from alien.resources.ncbi import NcbiMaps, NcbiRecord

from conftest import hgnc_maps, make_target, membership

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "legacy"))
import legacy_mapping
from legacy_support import legacy_hgnc_maps

GENES = [f"G{i}" for i in range(6)]
EXTRA = ["P0", "P1", "A0", "A1", "JUNK", "disease", "", "N0", "N1", "N2"]
ENSEMBL = [f"ENSG0000000000{i}" for i in range(6)]
RETIRED = ["ENSG00000000090", "ENSG00000000091"]
TAGS = ["DEPMAP", "REACTOME"]


@st.composite
def worlds(draw: st.DrawFn) -> dict[str, Any]:
    hgnc = []
    for index, symbol in enumerate(GENES):
        prev = "|".join(draw(st.lists(st.sampled_from(["P0", "P1"]), max_size=2, unique=True)))
        alias = "|".join(draw(st.lists(st.sampled_from(["A0", "A1"]), max_size=2, unique=True)))
        ensembl = draw(st.sampled_from([ENSEMBL[index], ENSEMBL[index], ""]))
        hgnc.append((symbol, prev, alias, ensembl))
    annotation: list[dict[str, Any]] = []
    for index, symbol in enumerate(GENES):
        for ensembl in draw(st.lists(st.sampled_from([ENSEMBL[index], *ENSEMBL]), min_size=0, max_size=2, unique=True)):
            annotation.append(  # noqa: PERF401 - each element draws fresh values
                {
                    "gene_symbol": symbol,
                    "ensembl_gene_id": ensembl,
                    "gene_biotype": draw(st.sampled_from(["protein_coding", "lncRNA"])),
                    "is_in_output_genes": draw(st.booleans()),
                    "target_id_metadata_source": draw(
                        st.sampled_from(["annotation_gtf", "metadata_fallback", "output_genes"])
                    ),
                }
            )
    output_ids = draw(st.one_of(st.none(), st.sets(st.sampled_from(ENSEMBL + RETIRED[:1]), min_size=1)))
    rows = [
        membership(
            f"T{draw(st.integers(0, 4))}",
            draw(st.sampled_from(GENES + EXTRA)),
            draw(st.sampled_from(["", "", *ENSEMBL, *RETIRED, f"{ENSEMBL[0]}.3"])),
            source_tag=draw(st.sampled_from(TAGS)),
        )
        for _ in range(draw(st.integers(1, 25)))
    ]
    ncbi_stages: dict[str, dict[str, Any]] = {"ncbi_gene_history": {}, "ncbi_gene_info_synonym": {}}
    ncbi_ambiguous: dict[str, dict[str, list[str]]] = {"ncbi_gene_history": {}, "ncbi_gene_info_synonym": {}}
    for symbol in ["N0", "N1", "N2", "JUNK"]:
        stage = draw(st.sampled_from(list(ncbi_stages)))
        if draw(st.booleans()):
            ncbi_ambiguous[stage][symbol] = ["G0", "G1"]
        else:
            target = draw(st.sampled_from(GENES))
            ensembl_ids = draw(st.lists(st.sampled_from(ENSEMBL), max_size=2, unique=True))
            ncbi_stages[stage][symbol] = (target, sorted(ensembl_ids))
    archive = {
        RETIRED[0]: {
            "found": True,
            "archive_id": RETIRED[0],
            "is_current": "0",
            "latest": "",
            "possible_replacement": draw(st.lists(st.sampled_from(ENSEMBL), max_size=2, unique=True)),
        },
        RETIRED[1]: {"found": False},
    }
    repairs = draw(
        st.dictionaries(st.sampled_from(["JUNK", "P0", "A1"]), st.sampled_from([*GENES, "NOPE"]), max_size=2)
    )
    return {
        "hgnc": hgnc,
        "annotation": annotation,
        "output_ids": output_ids,
        "rows": rows,
        "ncbi": (ncbi_stages, ncbi_ambiguous),
        "archive": archive,
        "repairs": repairs,
        "use_alias": draw(st.booleans()),
        "allowed": draw(st.sampled_from([["DEPMAP"], []])),
        "non_gene": ["disease"],
    }


class _Archive:
    def __init__(self, records: dict[str, dict[str, Any]]) -> None:
        self.records = records

    def prefetch(self, stable_ids: list[str], namespace: str = "") -> None:
        return None

    def lookup(self, stable_id: str) -> dict[str, Any]:
        return {"query_id": stable_id, **self.records.get(stable_id, {"found": False})}


def _legacy(world: dict[str, Any]) -> dict[str, Any]:
    stages, ambiguous = world["ncbi"]
    ncbi = {
        "stages": {
            stage: {
                s: {"current_symbol": t, "current_gene_ids": ["1"], "ensembl_gene_ids": e}
                for s, (t, e) in records.items()
            }
            for stage, records in stages.items()
        },
        "ambiguous": ambiguous,
    }
    cfg = {
        "gene_mapping": {
            "manual_symbol_repairs": world["repairs"],
            "use_hgnc_alias_mapping": world["use_alias"],
            "non_gene_tokens": world["non_gene"],
        },
        "ncbi_gene": {"allowed_source_tags": world["allowed"]},
    }
    target = make_target("t", world["annotation"], world["output_ids"])
    lists: dict[str, list[dict[str, Any]]] = {k: [] for k in ("unmapped", "ambiguous", "non_gene", "archive", "ncbi")}
    status: dict[tuple[str, str, str, str], int] = {}
    terms = legacy_mapping.map_ensembl(
        pd.DataFrame(world["rows"]),
        legacy_hgnc_maps(world["hgnc"]),
        ncbi,
        target.annotation,
        set(world["output_ids"]) if world["output_ids"] else None,
        "t",
        cfg,
        lists["unmapped"],
        lists["ambiguous"],
        lists["non_gene"],
        lists["archive"],
        lists["ncbi"],
        status,
        _Archive(world["archive"]),
        prefetch_archive=False,
    )
    return {
        "terms": {term_id: sorted(record["genes"]) for term_id, record in terms.items()},
        "status": {(key[2], key[3]): value for key, value in status.items() if value},
        **{name: [{k: str(v) for k, v in row.items()} for row in rows] for name, rows in lists.items()},
    }


def _new(world: dict[str, Any]) -> dict[str, Any]:
    stages, ambiguous = world["ncbi"]
    ncbi = NcbiMaps(
        stages={
            stage: {s: NcbiRecord(t, ("1",), tuple(e)) for s, (t, e) in records.items()}
            for stage, records in stages.items()
        },
        ambiguous=ambiguous,
    )
    config = load_config(
        {
            "sources": [{"type": "symbol_gmt", "path": "x.gmt"}],
            "outputs": {"include_symbols": False},
            "gene_mapping": {
                "manual_symbol_repairs": world["repairs"],
                "use_hgnc_alias_mapping": world["use_alias"],
                "non_gene_tokens": world["non_gene"],
            },
            "ncbi_gene": {"allowed_source_tags": world["allowed"]},
        }
    )
    target = make_target("t", world["annotation"], world["output_ids"])
    result = map_namespaces(
        pd.DataFrame(world["rows"]),
        config,
        [target],
        hgnc=hgnc_maps(world["hgnc"]),
        ncbi=ncbi,
        resolver=_Archive(world["archive"]),  # type: ignore[arg-type]
    )
    tables = {
        "unmapped": result.unmapped,
        "ambiguous": result.ambiguous,
        "non_gene": result.non_gene,
        "archive": result.archive_audit,
        "ncbi": result.ncbi_audit,
    }
    return {
        "terms": {term_id: sorted(record.genes) for term_id, record in result.terms["t"].items()},
        "status": {
            (row["source_tag"], row["mapping_status"]): int(row["n_gene_members"])
            for row in result.status_summary.to_dict("records")
        },
        **{
            name: [{k: str(v) for k, v in row.items()} for row in frame.to_dict("records")]
            for name, frame in tables.items()
        },
    }


@settings(max_examples=200, deadline=None, suppress_health_check=[HealthCheck.too_slow])
@given(worlds())
def test_new_projection_equals_legacy_map_ensembl(world: dict[str, Any]) -> None:
    legacy, new = _legacy(world), _new(world)
    for key in legacy:
        assert new[key] == legacy[key], key
