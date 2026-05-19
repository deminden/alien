from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from .config import load_config


@dataclass(frozen=True)
class Preset:
    name: str
    filename: str
    description: str
    aliases: tuple[str, ...] = ()

    @property
    def path(self) -> Path:
        return preset_dir() / self.filename


PRESETS: tuple[Preset, ...] = (
    Preset(
        name="pathways",
        filename="pathways.yml",
        description="Reactome, WikiPathways, KEGG MEDICUS, and GO biological process terms.",
    ),
    Preset(
        name="function_location",
        filename="function_location.yml",
        description="GO molecular function and cellular component terms.",
        aliases=("function",),
    ),
    Preset(
        name="disease_phenotype",
        filename="disease_phenotype.yml",
        description="HPO, DisGeNET, ClinVar, GWAS Catalog, and Jensen disease libraries.",
        aliases=("disease",),
    ),
    Preset(
        name="cancer_dependency",
        filename="cancer_dependency.yml",
        description="Cancer and dependency signatures for TCGA recount3.",
        aliases=("cancer",),
    ),
)

_PRESET_BY_NAME = {preset.name: preset for preset in PRESETS}
_PRESET_ALIASES = {alias: preset.name for preset in PRESETS for alias in preset.aliases}


def preset_dir() -> Path:
    return Path(__file__).with_name("presets")


def list_presets() -> tuple[Preset, ...]:
    return PRESETS


def resolve_preset(name: str) -> Preset:
    key = name.strip()
    canonical = _PRESET_ALIASES.get(key, key)
    try:
        return _PRESET_BY_NAME[canonical]
    except KeyError as exc:
        choices = sorted([preset.name for preset in PRESETS] + list(_PRESET_ALIASES))
        raise ValueError(f"Unknown ALIEN preset {name!r}. Available presets and aliases: {', '.join(choices)}") from exc


def is_preset_name(name: str) -> bool:
    return name.strip() in _PRESET_BY_NAME or name.strip() in _PRESET_ALIASES


def preset_path(name: str) -> Path:
    return resolve_preset(name).path


def load_preset_config(name: str, overrides: Sequence[str | Path | dict[str, Any]] | None = None) -> dict[str, Any]:
    return load_config(preset_path(name), overrides=overrides)
