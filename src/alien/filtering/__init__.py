"""Term filtering: size and broad-term filters, then redundancy removal."""

from __future__ import annotations

from .redundancy import FilterResult, RedundancyResult, filter_namespaces, redundancy_groups, remove_redundancy
from .size import SizeResult, apply_size_filters

__all__ = [
    "FilterResult",
    "RedundancyResult",
    "SizeResult",
    "apply_size_filters",
    "filter_namespaces",
    "redundancy_groups",
    "remove_redundancy",
]
