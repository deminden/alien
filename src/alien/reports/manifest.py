"""Term manifest (one row per term and namespace) and source manifest (one row per source collection)."""

from __future__ import annotations

import json
from typing import Any

import pandas as pd

from ..model import RemovalReason, TermSet
from ..schema import SOURCE_MANIFEST_COLUMNS, TERM_MANIFEST_COLUMNS

_META_COLUMNS = TERM_MANIFEST_COLUMNS[2:13]
_SOURCE_GROUP = ("source", "source_tag", "collection", "subcollection", "family", "aspect", "db_version")


def _removals(*frames: pd.DataFrame) -> dict[tuple[str, str], dict[str, Any]]:
    out: dict[tuple[str, str], dict[str, Any]] = {}
    for frame in frames:
        for raw in frame.to_dict("records"):
            row = {str(key): value for key, value in raw.items()}
            out.setdefault((str(row["target_namespace"]), str(row["term_id"])), row)
    return out


def build_term_manifest(
    memberships: pd.DataFrame,
    mapped: dict[str, TermSet],
    final: dict[str, TermSet],
    removed_size: pd.DataFrame,
    removed_redundancy: pd.DataFrame,
) -> pd.DataFrame:
    """One row per (term, namespace) with the namespace-specific outcome."""
    if memberships.empty:
        return pd.DataFrame(columns=list(TERM_MANIFEST_COLUMNS))
    meta = memberships.drop_duplicates("term_id").set_index("term_id")
    original_sizes = memberships.groupby("term_id")["gene_symbol"].nunique().to_dict()
    removed = _removals(removed_size, removed_redundancy)
    rows: list[dict[str, Any]] = []
    for term_id in sorted(meta.index):
        term_meta = meta.loc[term_id].to_dict()
        for namespace, terms in mapped.items():
            removal = removed.get((namespace, term_id), {})
            kept = term_id in final.get(namespace, {})
            reason = str(removal.get("reason", "")) or ("" if kept else str(RemovalReason.NO_MAPPED_GENES))
            rows.append(
                {
                    "term_id": term_id,
                    "target_namespace": namespace,
                    **{column: term_meta.get(column, "") for column in _META_COLUMNS},
                    "original_size_symbols": original_sizes.get(term_id, 0),
                    "size_after_mapping": len(terms[term_id].genes) if term_id in terms else 0,
                    "kept": kept,
                    "removed_reason": reason,
                    "representative_term_id": removal.get("representative_term_id", "") or term_id,
                    "redundancy_cluster_id": removal.get("redundancy_cluster_id", ""),
                }
            )
    return pd.DataFrame(rows, columns=list(TERM_MANIFEST_COLUMNS))


def _parse(value: object) -> dict[str, Any]:
    try:
        parsed = json.loads(str(value or "").strip() or "{}")
    except json.JSONDecodeError:
        return {}
    return parsed if isinstance(parsed, dict) else {}


def _joined(records: list[dict[str, Any]], key: str) -> str:
    return ";".join(sorted({str(r.get(key, "")).strip() for r in records if str(r.get(key, "")).strip()}))


def _list_values(records: list[dict[str, Any]], key: str) -> str:
    values: set[str] = set()
    for record in records:
        raw = record.get(key, [])
        if isinstance(raw, list):
            values.update(str(value).strip() for value in raw if str(value).strip())
        elif str(raw or "").strip():
            values.add(str(raw).strip())
    return json.dumps(sorted(values), ensure_ascii=True) if values else ""


def _single(series: pd.Series) -> str:
    values = sorted({str(value).strip() for value in series if str(value).strip()})
    return values[0] if len(values) == 1 else ";".join(values)


def build_source_manifest(memberships: pd.DataFrame) -> pd.DataFrame:
    """Resolved source collections: one row per (source, tag, collection, ...) with term and gene counts."""
    if memberships.empty:
        return pd.DataFrame(columns=list(SOURCE_MANIFEST_COLUMNS))
    frame = memberships.copy()
    for column in (*_SOURCE_GROUP, "source_url", "source_license_note", "metadata_json", "term_id", "gene_symbol"):
        if column not in frame:
            frame[column] = ""
    rows = []
    for key, group in frame.groupby(list(_SOURCE_GROUP), dropna=False, sort=True):
        records = [_parse(value) for value in group["metadata_json"].drop_duplicates()]
        rows.append(
            {
                **dict(zip(_SOURCE_GROUP, key, strict=True)),
                "selected_library": _joined(records, "selected_library"),
                "configured_name": _joined(records, "configured_name"),
                "match": _joined(records, "match"),
                "match_method": _joined(records, "match_method"),
                "candidate_libraries": _list_values(records, "candidate_libraries"),
                "n_terms": int(group["term_id"].nunique()),
                "n_memberships": len(group),
                "n_unique_symbols": int(group["gene_symbol"].nunique()),
                "source_url": _single(group["source_url"]),
                "source_license_note": _single(group["source_license_note"]),
            }
        )
    return pd.DataFrame(rows, columns=list(SOURCE_MANIFEST_COLUMNS))
