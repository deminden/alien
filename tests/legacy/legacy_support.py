# ruff: noqa
# mypy: ignore-errors
"""Minimal stand-ins for the 0.1.x utils used by the frozen legacy modules."""

from __future__ import annotations

import logging
import multiprocessing as mp
import re

import pandas as pd

LOGGER = logging.getLogger("alien.legacy")


def strip_ensembl_version(value):
    if value is None or pd.isna(value):
        return ""
    text = str(value).strip()
    if not text or text.upper() in {"NA", "NAN", "NULL"}:
        return ""
    match = re.match(r"^(ENSG[0-9]+)(?:\.[0-9]+)?(?:_PAR_Y)?$", text)
    return match.group(1) if match else text


def clean_field(value):
    text = "" if value is None or pd.isna(value) else str(value)
    return re.sub(r"[\t\r\n]+", " ", text).strip()


def description_from_record(record):
    fields = [
        ("source", record.get("source_tag") or record.get("source")),
        ("collection", record.get("collection")),
        ("aspect", record.get("aspect")),
        ("original_name", record.get("original_name")),
        ("db_version", record.get("db_version")),
    ]
    return ";".join(f"{key}={clean_field(value)}" for key, value in fields if clean_field(value))


def ensure_dirs(*paths):
    for path in paths:
        path.mkdir(parents=True, exist_ok=True)


def write_json(obj, path):
    raise NotImplementedError


def process_pool_context():
    return mp.get_context("spawn")


def legacy_hgnc_maps(rows):
    current, previous, alias, ensembl = {}, {}, {}, {}
    for symbol, prev, ali, ens in rows:
        current[symbol] = symbol
        if ens:
            ensembl[symbol] = ens
        for item in [x for x in prev.split("|") if x]:
            previous.setdefault(item, set()).add(symbol)
        for item in [x for x in ali.split("|") if x]:
            alias.setdefault(item, set()).add(symbol)
    return {
        "current": current,
        "previous": {k: next(iter(v)) for k, v in previous.items() if len(v) == 1},
        "alias": {k: next(iter(v)) for k, v in alias.items() if len(v) == 1},
        "ambiguous": {
            "previous": {k: sorted(v) for k, v in previous.items() if len(v) > 1},
            "alias": {k: sorted(v) for k, v in alias.items() if len(v) > 1},
        },
        "current_to_ensembl": ensembl,
    }
