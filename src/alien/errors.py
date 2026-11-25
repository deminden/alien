"""ALIEN exception hierarchy. The CLI maps each class to an exit code."""

from __future__ import annotations


class AlienError(Exception):
    """Base class for every error ALIEN raises on purpose."""

    exit_code = 1


class ConfigError(AlienError, ValueError):
    """The configuration is invalid or refers to something that does not exist."""

    exit_code = 2


class SourceError(AlienError, ValueError):
    """An input gene-set source or local data file is malformed or missing."""

    exit_code = 3


class TermIdCollisionError(SourceError):
    """Two sources produced the same term_id with different identities."""


class ResourceError(AlienError, RuntimeError):
    """A remote resource or cache could not be fetched, verified or read."""

    exit_code = 4
