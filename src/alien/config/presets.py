"""Bundled build presets."""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path

from ..errors import ConfigError
from . import loader
from .loader import ConfigInput, load_config
from .schema import AlienConfig


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
        description="MSigDB Hallmark, Reactome, WikiPathways, KEGG MEDICUS, and GO biological process terms.",
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
        name="cancer_programs",
        filename="cancer_programs.yml",
        description="MSigDB C6 oncogenic, C4 3CA tumour and C9 signatures, for tumour expression data.",
        aliases=("cancer", "cancer_dependency"),
    ),
    Preset(
        name="cell_line_signatures",
        filename="cell_line_signatures.yml",
        description="Opt-in per-cell-line sets: DepMap CRISPR dependencies, CCLE and NCI-60 expression.",
        aliases=("cell_lines",),
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
    except KeyError:
        choices = sorted([preset.name for preset in PRESETS] + list(_PRESET_ALIASES))
        raise ConfigError(
            f"Unknown ALIEN preset {name!r}. Available presets and aliases: {', '.join(choices)}"
        ) from None


def is_preset_name(name: str) -> bool:
    return name.strip() in _PRESET_BY_NAME or name.strip() in _PRESET_ALIASES


def preset_path(name: str) -> Path:
    return resolve_preset(name).path


def load_preset_config(name: str, overrides: Sequence[ConfigInput] | None = None) -> AlienConfig:
    return load_config(preset_path(name), overrides=overrides, resolve_relative_paths=False)


def _preset_or_none(name: str) -> Path | None:
    return preset_path(name) if is_preset_name(name) else None


loader.PRESET_RESOLVER = _preset_or_none
