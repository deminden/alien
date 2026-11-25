"""Configuration: typed schema, layered loading and bundled presets."""

from __future__ import annotations

from .loader import (
    ConfigInput,
    finalize_config,
    load_config,
    merge_config,
    read_config_file,
    read_layer,
    validate_config,
    validate_layers,
)
from .presets import Preset, is_preset_name, list_presets, load_preset_config, preset_path, resolve_preset
from .schema import AlienConfig

__all__ = [
    "AlienConfig",
    "ConfigInput",
    "Preset",
    "finalize_config",
    "is_preset_name",
    "list_presets",
    "load_config",
    "load_preset_config",
    "merge_config",
    "preset_path",
    "read_config_file",
    "read_layer",
    "resolve_preset",
    "validate_config",
    "validate_layers",
]
