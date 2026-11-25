"""The one HTTP client: retries with backoff (honouring Retry-After), streaming, hashing, atomic replace."""

from __future__ import annotations

import hashlib
import logging
import os
import tempfile
from collections.abc import Callable
from dataclasses import asdict, dataclass
from datetime import UTC, datetime
from pathlib import Path
from typing import Any, Final

import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

from ..errors import ResourceError
from ..version import __version__

LOGGER = logging.getLogger("alien")
USER_AGENT: Final = f"bioalien/{__version__} (+https://github.com/deminden/alien)"
RETRY_STATUSES: Final = (429, 500, 502, 503, 504)
CHUNK: Final = 1024 * 1024


@dataclass(frozen=True)
class FetchRecord:
    """What was downloaded: enough to pin and verify it later."""

    url: str
    final_url: str
    sha256: str
    size: int
    etag: str = ""
    last_modified: str = ""
    retrieved_at: str = ""
    md5: str = ""

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


def build_session(retries: int = 3, backoff: float = 2.0) -> requests.Session:
    session = requests.Session()
    retry = Retry(
        total=retries,
        connect=retries,
        read=retries,
        backoff_factor=backoff,
        status_forcelist=RETRY_STATUSES,
        allowed_methods=frozenset({"GET", "HEAD", "POST"}),
        respect_retry_after_header=True,
    )
    adapter = HTTPAdapter(max_retries=retry)
    session.mount("http://", adapter)
    session.mount("https://", adapter)
    session.headers["User-Agent"] = USER_AGENT
    return session


#: Replaced in tests to serve fixtures without a network.
SESSION_FACTORY: Callable[..., requests.Session] = build_session


def fetch(
    url: str,
    dest: Path,
    *,
    sha256: str | None = None,
    md5: str | None = None,
    timeout: float = 3600,
    connect_timeout: float = 15,
) -> FetchRecord:
    """Stream ``url`` to a unique temporary file beside ``dest``, verify it, then ``os.replace`` it."""
    dest.parent.mkdir(parents=True, exist_ok=True)
    LOGGER.info("Downloading %s", url)
    sha = hashlib.sha256()
    md5_digest = hashlib.md5(usedforsecurity=False)
    size = 0
    fd, tmp_name = tempfile.mkstemp(prefix=f".{dest.name}.", suffix=".part", dir=dest.parent)
    tmp = Path(tmp_name)
    try:
        with (
            SESSION_FACTORY() as session,
            session.get(url, stream=True, timeout=(connect_timeout, timeout)) as response,
        ):
            response.raise_for_status()
            with os.fdopen(fd, "wb") as handle:
                for chunk in response.iter_content(chunk_size=CHUNK):
                    if chunk:
                        sha.update(chunk)
                        md5_digest.update(chunk)
                        size += len(chunk)
                        handle.write(chunk)
            headers = response.headers
            final_url = str(response.url or url)
        _verify(url, "sha256", sha.hexdigest(), sha256)
        _verify(url, "MD5", md5_digest.hexdigest(), md5)
        tmp.replace(dest)
    except requests.RequestException as exc:
        tmp.unlink(missing_ok=True)
        raise ResourceError(f"Could not download {url}: {exc}") from exc
    except BaseException:
        tmp.unlink(missing_ok=True)
        raise
    return FetchRecord(
        url=url,
        final_url=final_url,
        sha256=sha.hexdigest(),
        size=size,
        etag=str(headers.get("ETag", "")),
        last_modified=str(headers.get("Last-Modified", "")),
        retrieved_at=datetime.now(UTC).isoformat(timespec="seconds"),
        md5=md5_digest.hexdigest(),
    )


def _verify(url: str, kind: str, actual: str, expected: str | None) -> None:
    if expected and actual.lower() != expected.lower():
        raise ResourceError(f"Download from {url} failed {kind} verification (expected {expected}, got {actual}).")


def get_json(url: str, *, timeout: float = 120, connect_timeout: float = 15) -> Any:
    try:
        with SESSION_FACTORY() as session:
            response = session.get(url, timeout=(connect_timeout, timeout))
            response.raise_for_status()
            return response.json()
    except (requests.RequestException, ValueError) as exc:
        raise ResourceError(f"Could not read {url}: {exc}") from exc


def post_json(
    url: str,
    payload: Any,
    *,
    timeout: float = 120,
    connect_timeout: float = 30,
    retries: int = 3,
    backoff: float = 2.0,
) -> Any:
    try:
        with SESSION_FACTORY(retries=retries, backoff=backoff) as session:
            response = session.post(
                url,
                json=payload,
                headers={"Content-Type": "application/json", "Accept": "application/json"},
                timeout=(connect_timeout, timeout),
            )
            response.raise_for_status()
            return response.json()
    except (requests.RequestException, ValueError) as exc:
        raise ResourceError(f"POST {url} failed: {exc}") from exc


def file_digests(path: Path) -> tuple[str, str]:
    """(sha256, md5) of a file on disk."""
    sha = hashlib.sha256()
    md5_digest = hashlib.md5(usedforsecurity=False)
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(CHUNK), b""):
            sha.update(chunk)
            md5_digest.update(chunk)
    return sha.hexdigest(), md5_digest.hexdigest()
