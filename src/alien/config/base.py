"""Shared pydantic base model and validation helpers for the configuration schema."""

from __future__ import annotations

import warnings
from typing import Annotated, Any

from pydantic import BaseModel, ConfigDict, Field

from .removed import REMOVED_OPTIONS

PositiveInt = Annotated[int, Field(ge=1)]
NonNegativeInt = Annotated[int, Field(ge=0)]


def _deprecated(old: str, new: str) -> None:
    warnings.warn(f"Config key {old!r} is deprecated; use {new!r}.", DeprecationWarning, stacklevel=4)


def _rename(data: Any, aliases: dict[str, str], where: str) -> Any:
    """Map deprecated key spellings to their canonical names, warning once per key."""
    if not isinstance(data, dict):
        return data
    data = dict(data)
    for old, new in aliases.items():
        if old in data:
            _deprecated(f"{where}.{old}", f"{where}.{new}")
            value = data.pop(old)
            data.setdefault(new, value)
    return data


def _reject_removed(section: str, data: Any) -> Any:
    if isinstance(data, dict):
        for key, replacement in REMOVED_OPTIONS.get(section, {}).items():
            if key in data:
                raise ValueError(f"option {section}.{key} was removed: {replacement}")
    return data


def _blank_to_none(value: Any) -> Any:
    if isinstance(value, str) and value.strip() in {"", "null"}:
        return None
    return value


class Model(BaseModel):
    model_config = ConfigDict(extra="forbid", frozen=True, populate_by_name=True, use_attribute_docstrings=True)
