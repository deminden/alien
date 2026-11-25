"""ALIEN: Audited Library Integration for External Namespaces."""

from __future__ import annotations

import logging

from .config import AlienConfig, list_presets, load_config, load_preset_config, preset_path
from .pipeline import BuildResult, build, fetch
from .version import __version__

logging.getLogger("alien").addHandler(logging.NullHandler())

__all__ = [
    "AlienConfig",
    "BuildResult",
    "__version__",
    "build",
    "fetch",
    "list_presets",
    "load_config",
    "load_preset_config",
    "preset_path",
]
