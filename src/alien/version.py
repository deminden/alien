"""The installed ALIEN version, read from package metadata (single source: pyproject.toml)."""

from __future__ import annotations

from importlib.metadata import PackageNotFoundError, version

try:
    __version__ = version("bioalien")
except PackageNotFoundError:  # pragma: no cover - running from a source tree without installation
    __version__ = "0+unknown"
