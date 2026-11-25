"""Enrichr libraries resolved by exact name or regex against Enrichr's library index."""

from __future__ import annotations

import json
import logging
import re
from functools import partial
from pathlib import Path
from typing import Final
from urllib.parse import quote

import pandas as pd

from ..cache import ResourceCache
from ..config.schema import EnrichrLibrary, EnrichrRemoteSource
from ..errors import ResourceError, SourceError
from ..model import sanitize_term_id
from .frames import concat_memberships, rows_from_gmt
from .registry import register

LOGGER = logging.getLogger("alien")
ENRICHR_BASE_URL: Final = "https://maayanlab.cloud/Enrichr"
DATASET_STATISTICS_ENDPOINT: Final = f"{ENRICHR_BASE_URL}/datasetStatistics"
LIBRARY_DOWNLOAD_ENDPOINT: Final = f"{ENRICHR_BASE_URL}/geneSetLibrary?mode=text&libraryName={{library}}"
ENRICHR_LICENSE: Final = "Enrichr library provenance and upstream licenses apply."


def library_names(stats: object) -> list[str]:
    if isinstance(stats, dict) and "statistics" in stats:
        stats = stats["statistics"]
    items = stats if isinstance(stats, list) else []
    names = [
        item.get("libraryName") or item.get("library_name") or item.get("name")
        for item in items
        if isinstance(item, dict)
    ]
    return sorted({str(name) for name in names if name})


def _latest_year(name: str) -> int:
    years = [int(year) for year in re.findall(r"(20[0-9]{2})", name)]
    return max(years) if years else 0


def resolve_enrichr_library(names: list[str], library: EnrichrLibrary) -> tuple[str, list[str], str]:
    """(selected library, candidates, match method). A regex picks the newest year, then the shortest name."""
    name = (library.name or "").strip()
    if name and name in names:
        return name, [name], "exact"
    match = (library.match or "").strip()
    if match:
        pattern = re.compile(match, re.IGNORECASE)
        candidates = sorted({candidate for candidate in names if pattern.search(candidate)})
        if candidates:
            return min(candidates, key=lambda value: (-_latest_year(value), len(value), value)), candidates, "regex"
        raise SourceError(f"No Enrichr libraries matched {name or match!r}.")
    raise SourceError(f"Enrichr library {name!r} was not found in datasetStatistics metadata.")


def _non_empty(path: Path) -> None:
    if not path.read_text(encoding="utf-8").strip():
        path.unlink()
        raise ResourceError(f"Downloaded Enrichr library {path.name} was empty.")


def _library_frame(
    spec: EnrichrRemoteSource, library: EnrichrLibrary, names: list[str], cache: ResourceCache, cache_dir: Path
) -> pd.DataFrame:
    resolved, candidates, method = resolve_enrichr_library(names, library)
    # Frozen builds take the library chosen when the lock was written, even if a newer one now matches.
    selected = cache.resolve(f"enrichr:library:{library.name or ''}:{library.match or ''}", partial(str, resolved))
    endpoint = spec.library_download_endpoint or LIBRARY_DOWNLOAD_ENDPOINT
    gmt_path = cache.fetch(
        f"enrichr:{selected}",
        endpoint.format(library=quote(selected)),
        cache_dir / f"{sanitize_term_id(selected)}.gmt",
        timeout=spec.timeout_seconds,
        upstream_version=selected,
        validate=_non_empty,
    )
    metadata = {
        "selected_library": selected,
        "configured_name": library.name or "",
        "match": library.match or "",
        "match_method": method,
        "candidate_libraries": candidates,
    }
    constants = {
        "source": "Enrichr",
        "source_tag": library.source_tag or spec.source_tag or sanitize_term_id(selected).upper(),
        "collection": selected,
        "subcollection": "",
        "family": library.family or spec.family or "biology_process_pathway",
        "aspect": library.aspect or spec.aspect or "library",
        "db_version": "",
        "source_url": ENRICHR_BASE_URL,
        "source_license_note": ENRICHR_LICENSE,
        "metadata_json": json.dumps(metadata, sort_keys=True),
    }
    return rows_from_gmt(gmt_path, f"ENRICHR_{sanitize_term_id(selected)}", constants)


@register("enrichr_remote")
def load_enrichr_remote(spec: EnrichrRemoteSource, cache: ResourceCache) -> pd.DataFrame:
    cache_dir = spec.cache_dir or cache.root / "enrichr"
    index_path = cache.fetch(
        "enrichr:datasetStatistics",
        spec.dataset_statistics_endpoint or DATASET_STATISTICS_ENDPOINT,
        cache_dir / "datasetStatistics.json",
        timeout=spec.timeout_seconds,
    )
    names = library_names(json.loads(index_path.read_text(encoding="utf-8")))
    if not names:
        raise SourceError("Enrichr metadata did not contain any library names.")
    return concat_memberships(_library_frame(spec, library, names, cache, cache_dir) for library in spec.libraries)
