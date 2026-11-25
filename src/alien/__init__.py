from __future__ import annotations

from .build import BuildResult, build
from .config import load_config
from .presets import list_presets, load_preset_config, preset_path

__all__ = ["BuildResult", "build", "list_presets", "load_config", "load_preset_config", "preset_path"]
__version__ = "0.1.6"
