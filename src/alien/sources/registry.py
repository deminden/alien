"""The source-type registry: ``type:`` value -> loader returning canonical memberships."""

from __future__ import annotations

from collections.abc import Callable
from typing import Any

import pandas as pd

from ..cache import ResourceCache

Loader = Callable[[Any, ResourceCache], pd.DataFrame]
_REGISTRY: dict[str, Loader] = {}


def register[SpecT](
    *type_names: str,
) -> Callable[[Callable[[SpecT, ResourceCache], pd.DataFrame]], Callable[[SpecT, ResourceCache], pd.DataFrame]]:
    """Register a loader for one or more ``type:`` values."""

    def decorator(
        loader: Callable[[SpecT, ResourceCache], pd.DataFrame],
    ) -> Callable[[SpecT, ResourceCache], pd.DataFrame]:
        for name in type_names:
            _REGISTRY[name] = loader
        return loader

    return decorator


def loader_for(type_name: str) -> Loader:
    return _REGISTRY[type_name]


def registered_types() -> list[str]:
    return sorted(_REGISTRY)
