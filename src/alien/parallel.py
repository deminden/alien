"""Process-pool context shared by every parallel stage."""

from __future__ import annotations

import multiprocessing as mp
from multiprocessing.context import BaseContext


def process_pool_context() -> BaseContext:
    """forkserver (or spawn): never raw fork from parents that may already run helper threads."""
    for method in ("forkserver", "spawn"):
        try:
            return mp.get_context(method)
        except ValueError:
            continue
    return mp.get_context()
