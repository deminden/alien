"""Ensembl archive stable-ID lookups for retired source Ensembl IDs, cached on disk."""

from __future__ import annotations

import hashlib
import json
import logging
from collections.abc import Iterable
from pathlib import Path
from typing import Any, Final

from ..cache import ResourceCache
from ..config.schema import EnsemblArchiveConfig
from ..errors import ResourceError
from ..io import http
from ..io.atomic import atomic_write_json
from ..model import ensembl_stable_id

LOGGER = logging.getLogger("alien")
CACHE_KEY: Final = "ensembl_archive:lookups"
TRUE_VALUES: Final = frozenset({"1", "true", "True"})
#: Bookkeeping fields that do not change what a lookup means for mapping.
UNPINNED_FIELDS: Final = frozenset({"checked_release"})


def normalize_archive_record(stable_id: str, item: dict[str, Any]) -> dict[str, Any]:
    """Keep only the stable fields needed for rescue and audit."""
    replacements = []
    for replacement in item.get("possible_replacement") or []:
        replacement_id = (
            replacement.get("id") or replacement.get("stable_id") if isinstance(replacement, dict) else replacement
        )
        if stable := ensembl_stable_id(replacement_id):
            replacements.append(stable)
    return {
        "query_id": stable_id,
        "found": True,
        "archive_id": ensembl_stable_id(item.get("id", stable_id)),
        "is_current": str(item.get("is_current", "")),
        "release": str(item.get("release", "")),
        "assembly": str(item.get("assembly", "")),
        "latest": ensembl_stable_id(item.get("latest", "")),
        "possible_replacement": sorted(set(replacements)),
    }


def archive_target_matches(info: dict[str, Any], target_ids: frozenset[str] | set[str]) -> list[str]:
    """Archive replacement candidates that are in the target ID set."""
    if not info or not info.get("found"):
        return []
    candidates = set(info.get("possible_replacement") or [])
    if latest := ensembl_stable_id(info.get("latest", "")):
        candidates.add(latest)
    archive_id = ensembl_stable_id(info.get("archive_id", ""))
    if str(info.get("is_current", "")) in TRUE_VALUES and archive_id:
        candidates.add(archive_id)
    return sorted(candidate for candidate in candidates if candidate in target_ids)


def unresolved_status(info: dict[str, Any], matches: list[str], fallback: str) -> str:
    if matches:
        return f"archive_ambiguous_target_replacement__{fallback}"
    if info.get("error"):
        return f"archive_query_failed__{fallback}"
    if info and not info.get("found"):
        return f"archive_not_found__{fallback}"
    if info:
        return f"archive_no_target_replacement__{fallback}"
    return f"archive_not_checked__{fallback}"


class EnsemblArchiveResolver:
    """Look up retired Ensembl IDs; results are cached and saved atomically after every batch.

    "Not found" answers are stamped with the Ensembl release that gave them and are re-checked once a
    later release is available.
    """

    def __init__(self, config: EnsemblArchiveConfig, cache: ResourceCache) -> None:
        self.config = config
        self.resources = cache
        self.cache_path = cache.root / "ensembl_archive" / "archive_id_cache.json"
        self.cache: dict[str, dict[str, Any]] = self._load()
        self._current_release: int | None = None

    def _load(self) -> dict[str, dict[str, Any]]:
        if not self.cache_path.exists():
            return {}
        try:
            data = json.loads(self.cache_path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError) as exc:
            LOGGER.warning("Ignoring unreadable Ensembl archive cache %s: %s", self.cache_path, exc)
            return {}
        return data if isinstance(data, dict) else {}

    def lookup(self, stable_id: str) -> dict[str, Any]:
        return self.cache.get(stable_id, {"query_id": stable_id, "found": False})

    def current_release(self) -> int | None:
        if self._current_release is None:
            base = self.config.endpoint.rsplit("/archive", 1)[0]
            try:
                releases = http.get_json(f"{base}/info/data/?content-type=application/json").get("releases", [])
                self._current_release = max(int(release) for release in releases) if releases else None
            except (ResourceError, AttributeError, TypeError, ValueError) as exc:
                LOGGER.warning("Could not read the current Ensembl release: %s", exc)
        return self._current_release

    def _stale(self, stable_ids: Iterable[str]) -> list[str]:
        not_found = [
            i
            for i in stable_ids
            if i in self.cache and not self.cache[i].get("found") and not self.cache[i].get("error")
        ]
        if not not_found or not self.resources.network_allowed:
            return []
        release = self.current_release()
        if release is None:
            return []
        return [i for i in not_found if int(self.cache[i].get("checked_release") or 0) < release]

    def prefetch(self, stable_ids: list[str], namespace: str = "") -> None:
        """Fetch uncached IDs, retry cached failures, and re-check "not found" answers from older releases."""
        label = f" for {namespace}" if namespace else ""
        missing = [i for i in stable_ids if i not in self.cache or self.cache[i].get("error")]
        missing += self._stale(stable_ids)
        if not missing:
            LOGGER.info("Ensembl archive audit%s: %d problematic IDs already cached", label, len(stable_ids))
            self._pin(stable_ids, namespace)
            return
        if self.resources.frozen:
            raise ResourceError(
                f"{len(missing)} Ensembl IDs{label} are not in the archive lookup cache and --frozen forbids network "
                "access; run `alien fetch` to update it."
            )
        if not self.resources.network_allowed:
            LOGGER.warning(
                "Offline: %d Ensembl IDs%s are not in the archive cache and are left unchecked.", len(missing), label
            )
            for stable_id in missing:
                self.cache.setdefault(stable_id, {"query_id": stable_id, "found": False, "error": "offline"})
            return
        LOGGER.info("Auditing %d Ensembl IDs%s against the Ensembl archive endpoint", len(missing), label)
        self.current_release()  # stamps "not found" answers so a later release re-checks them
        for start in range(0, len(missing), self.config.batch_size):
            batch = missing[start : start + self.config.batch_size]
            LOGGER.info(
                "Fetching Ensembl archive batch%s %d-%d of %d", label, start + 1, start + len(batch), len(missing)
            )
            self._fetch_batch(batch)
            self._save()
        self._pin(stable_ids, namespace)

    def _pin(self, stable_ids: list[str], namespace: str) -> None:
        """Pin the lookups this namespace uses. The cache file is shared and grows, so its own hash is not pinned."""
        records = {
            stable_id: {key: value for key, value in self.lookup(stable_id).items() if key not in UNPINNED_FIELDS}
            for stable_id in sorted(set(stable_ids))
        }
        digest = hashlib.sha256(json.dumps(records, sort_keys=True, default=str).encode("utf-8")).hexdigest()
        key = f"{CACHE_KEY}:{namespace}" if namespace else CACHE_KEY
        self.resources.record_content(key, self.config.endpoint, self.cache_path, digest, len(records))

    def _fetch_batch(self, stable_ids: list[str]) -> None:
        try:
            payload = http.post_json(
                self.config.endpoint,
                {"id": stable_ids},
                timeout=self.config.timeout_seconds,
                connect_timeout=self.config.connect_timeout_seconds,
                retries=self.config.retry_attempts,
                backoff=self.config.retry_backoff_seconds,
            )
        except ResourceError as exc:
            if len(stable_ids) > 1:
                middle = len(stable_ids) // 2
                LOGGER.warning("Splitting failed Ensembl archive batch of %d IDs: %s", len(stable_ids), exc)
                self._fetch_batch(stable_ids[:middle])
                self._fetch_batch(stable_ids[middle:])
                return
            LOGGER.warning("Ensembl archive lookup failed for %s: %s", stable_ids[0], exc)
            self.cache[stable_ids[0]] = {"query_id": stable_ids[0], "found": False, "error": str(exc)}
            return
        returned: set[str] = set()
        for item in payload if isinstance(payload, list) else []:
            if isinstance(item, dict) and "id" in item:
                stable_id = ensembl_stable_id(item["id"])
                returned.add(stable_id)
                self.cache[stable_id] = normalize_archive_record(stable_id, item)
        release = self._current_release
        for stable_id in set(stable_ids) - returned:
            record: dict[str, Any] = {"query_id": stable_id, "found": False}
            if release is not None:
                record["checked_release"] = release
            self.cache[stable_id] = record

    def _save(self) -> None:
        with self.resources.file_lock(self.cache_path):
            on_disk = self._load()
            on_disk.update(self.cache)  # keep entries written concurrently by other builds
            persisted = {key: value for key, value in on_disk.items() if value.get("error") != "offline"}
            atomic_write_json(self.cache_path, dict(sorted(persisted.items())))
            self.cache = on_disk


def build_archive_resolver(config: EnsemblArchiveConfig, cache: ResourceCache) -> EnsemblArchiveResolver | None:
    return EnsemblArchiveResolver(config, cache) if config.enabled else None


def cache_path(source_dir: Path) -> Path:
    return source_dir / "ensembl_archive" / "archive_id_cache.json"
