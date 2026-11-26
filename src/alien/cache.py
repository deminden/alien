"""Resource cache: file layout, cross-process locking, the lock file, and offline/frozen modes.

* ``offline``: never touch the network; everything needed must already be cached.
* ``frozen``: offline, and every cached resource must match ``alien.lock.json`` byte for byte.
* ``alien fetch`` downloads what a config needs and writes or updates the lock file.
"""

from __future__ import annotations

import hashlib
import json
import logging
from collections.abc import Callable, Iterator, Mapping, Sequence
from contextlib import contextmanager, suppress
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Final, TypedDict, Unpack

from filelock import FileLock

from .errors import ResourceError
from .io import http
from .io.atomic import atomic_write_json
from .schema import CACHE_SCHEMA_VERSION

LOGGER = logging.getLogger("alien")
LOCK_FILE_NAME: Final = "alien.lock.json"
LOCK_FORMAT_VERSION: Final = 1
LOCK_TIMEOUT_SECONDS: Final = 3600


def file_sha256(path: Path) -> str:
    """sha256 of ``path``, memoised in a ``.sha256`` sidecar keyed by size and mtime."""
    stat = path.stat()
    sidecar = path.with_name(f".{path.name}.sha256")
    stamp = f"{stat.st_size}:{stat.st_mtime_ns}"
    if sidecar.exists():
        cached_stamp, _, cached_digest = sidecar.read_text(encoding="utf-8").strip().partition(" ")
        if cached_stamp == stamp and len(cached_digest) == 64:  # noqa: PLR2004 - hex sha256 length
            return cached_digest
    digest, _ = http.file_digests(path)
    with suppress(OSError):  # read-only caches still work; the digest is just recomputed next time
        sidecar.write_text(f"{stamp} {digest}\n", encoding="utf-8")
    return digest


def derived_cache_path(directory: Path, label: str, inputs: Mapping[str, Any], suffix: str = ".parquet") -> Path:
    """A cache file name keyed by everything that affects its content, including the cache schema version."""
    payload = json.dumps({"cache_schema_version": CACHE_SCHEMA_VERSION, **inputs}, sort_keys=True, default=str)
    digest = hashlib.sha256(payload.encode("utf-8")).hexdigest()[:16]
    return directory / f"{label}__{digest}{suffix}"


class FetchOptions(TypedDict, total=False):
    """Optional arguments of :meth:`ResourceCache.fetch`."""

    sha256: str | None
    md5: str | None
    timeout: float
    upstream_version: str
    validate: Callable[[Path], None] | None


@dataclass
class LockFile:
    """``alien.lock.json``: pinned remote resources and resolved choices (such as Enrichr library names)."""

    path: Path | None = None
    resources: dict[str, dict[str, Any]] = field(default_factory=dict)
    resolutions: dict[str, str] = field(default_factory=dict)

    @classmethod
    def load(cls, path: Path) -> LockFile:
        if not path.exists():
            return cls(path=path)
        try:
            data = json.loads(path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError) as exc:
            raise ResourceError(f"Lock file {path} is unreadable: {exc}") from exc
        if not isinstance(data, dict) or data.get("lock_version") != LOCK_FORMAT_VERSION:
            raise ResourceError(f"Lock file {path} has an unsupported format.")
        return cls(path=path, resources=dict(data.get("resources", {})), resolutions=dict(data.get("resolutions", {})))

    def to_dict(self) -> dict[str, Any]:
        return {
            "lock_version": LOCK_FORMAT_VERSION,
            "resources": dict(sorted(self.resources.items())),
            "resolutions": dict(sorted(self.resolutions.items())),
        }

    def save(self) -> None:
        if self.path is None:
            raise ResourceError("No lock file path configured.")
        atomic_write_json(self.path, self.to_dict())


@dataclass
class ResourceCache:
    """Every remote resource goes through :meth:`fetch`; every pinned choice through :meth:`resolve`."""

    root: Path
    force: bool = False
    offline: bool = False
    frozen: bool = False
    lock: LockFile = field(default_factory=LockFile)
    used: dict[str, dict[str, Any]] = field(default_factory=dict)

    @property
    def network_allowed(self) -> bool:
        return not (self.offline or self.frozen)

    @contextmanager
    def file_lock(self, path: Path) -> Iterator[None]:
        """Serialize writers of ``path`` across processes (Snakemake, Nextflow, SLURM arrays)."""
        path.parent.mkdir(parents=True, exist_ok=True)
        with FileLock(str(path.with_name(f".{path.name}.lock")), timeout=LOCK_TIMEOUT_SECONDS):
            yield

    def fetch(self, key: str, urls: str | Sequence[str], dest: Path, **options: Unpack[FetchOptions]) -> Path:
        """Return ``dest``, downloading it from the first working URL unless it is cached."""
        sha256, md5 = options.get("sha256"), options.get("md5")
        timeout, upstream_version = options.get("timeout", 3600), options.get("upstream_version", "")
        validate = options.get("validate")
        candidates = [urls] if isinstance(urls, str) else list(urls)
        with self.file_lock(dest):
            if dest.exists() and not self.force:
                self._check_cached(key, candidates[0], dest, sha256=sha256, md5=md5, upstream_version=upstream_version)
                return dest
            if not self.network_allowed:
                mode = "--frozen" if self.frozen else "--offline"
                raise ResourceError(
                    f"{key}: {dest} is not cached and network access is disabled ({mode}). "
                    "Run `alien fetch CONFIG` first."
                )
            record = self._download(key, candidates, dest, sha256=sha256, md5=md5, timeout=timeout)
            if validate is not None:
                validate(dest)
            self._record(key, record.to_dict(), dest, upstream_version)
        return dest

    def _download(
        self, key: str, urls: Sequence[str], dest: Path, *, sha256: str | None, md5: str | None, timeout: float
    ) -> http.FetchRecord:
        errors: list[str] = []
        for url in urls:
            try:
                return http.fetch(url, dest, sha256=sha256, md5=md5, timeout=timeout)
            except ResourceError as exc:
                errors.append(str(exc))
                LOGGER.warning("%s", exc)
        raise ResourceError(f"{key}: every URL failed:\n  " + "\n  ".join(errors))

    def _check_cached(
        self, key: str, url: str, dest: Path, *, sha256: str | None, md5: str | None, upstream_version: str
    ) -> None:
        pinned = self.lock.resources.get(key)
        digest = file_sha256(dest)
        if sha256 and digest != sha256.lower():
            raise ResourceError(f"{key}: cached {dest} failed sha256 verification; delete it to re-download.")
        if md5 and http.file_digests(dest)[1] != md5.lower():
            raise ResourceError(f"{key}: cached {dest} failed MD5 verification; delete it to re-download.")
        if self.frozen:
            if pinned is None:
                raise ResourceError(f"{key} is not pinned in {self.lock.path}; run `alien fetch` to update the lock.")
            if pinned.get("sha256") != digest:
                raise ResourceError(
                    f"{key}: cached {dest} differs from {self.lock.path} "
                    f"(sha256 expected {pinned.get('sha256')}, found {digest}). "
                    "Restore the file or re-run `alien fetch` to update the lock."
                )
        record: dict[str, Any] = (
            dict(pinned) if pinned and pinned.get("sha256") == digest else {"url": url, "sha256": digest}
        )
        record.setdefault("size", dest.stat().st_size)
        self._record(key, record, dest, upstream_version)

    def _record(self, key: str, record: dict[str, Any], dest: Path, upstream_version: str) -> None:
        entry = {**record, "path": self._relative(dest)}
        if upstream_version:
            entry["upstream_version"] = upstream_version
        self.used[key] = entry
        if not self.frozen:
            self.lock.resources[key] = entry

    def record_content(self, key: str, url: str, dest: Path, digest: str, n_records: int) -> None:
        """Pin the records a build read from a cache file that other builds also add to.

        ``digest`` covers only those records, so the pin survives unrelated additions to the file.
        """
        if self.frozen:
            pinned = self.lock.resources.get(key)
            if pinned is None:
                raise ResourceError(f"{key} is not pinned in {self.lock.path}; run `alien fetch` to update the lock.")
            if pinned.get("sha256") != digest:
                raise ResourceError(
                    f"{key}: the records this build reads from {dest} differ from {self.lock.path}; "
                    "re-run `alien fetch`."
                )
        self._record(key, {"url": url, "sha256": digest, "n_records": n_records, "scope": "records_used"}, dest, "")

    def resolve(self, key: str, compute: Callable[[], str]) -> str:
        """A pinned choice: frozen builds reuse the lock; other builds compute and record it."""
        if self.frozen:
            if key not in self.lock.resolutions:
                raise ResourceError(f"{key} is not pinned in {self.lock.path}; run `alien fetch` to update the lock.")
            return self.lock.resolutions[key]
        value = compute()
        self.lock.resolutions[key] = value
        return value

    def _relative(self, path: Path) -> str:
        try:
            return path.resolve().relative_to(self.root.resolve()).as_posix()
        except ValueError:
            return str(path)
