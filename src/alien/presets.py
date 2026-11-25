"""Deprecated location; import from :mod:`alien.config` instead. Removed in 0.3."""

from __future__ import annotations

import warnings

from .config.presets import (
    PRESETS,
    Preset,
    is_preset_name,
    list_presets,
    load_preset_config,
    preset_dir,
    preset_path,
    resolve_preset,
)

warnings.warn("alien.presets is deprecated; import from alien.config instead.", DeprecationWarning, stacklevel=2)

__all__ = [
    "PRESETS",
    "Preset",
    "is_preset_name",
    "list_presets",
    "load_preset_config",
    "preset_dir",
    "preset_path",
    "resolve_preset",
]
