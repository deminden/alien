#!/usr/bin/env python3
"""Redundancy-filter benchmark: N synthetic terms in one family, realistic sizes and hub genes.

Budget (CI, non-blocking): 50,000 terms in under 60 s and under 2 GB peak RSS.
Usage: python benchmarks/redundancy_benchmark.py --terms 50000 [--seed 1]
"""

from __future__ import annotations

import argparse
import resource
import sys
import time

import numpy as np

from alien.config.schema import RedundancyConfig
from alien.filtering.redundancy import remove_redundancy
from alien.model import TermRecord

TIME_BUDGET_SECONDS = 60.0
MEMORY_BUDGET_GB = 2.0
NEAR_DUPLICATE_RATE = 0.15


def synthetic_terms(n_terms: int, seed: int, universe: int = 20_000) -> dict[str, TermRecord]:
    rng = np.random.default_rng(seed)
    popularity = 1.0 / np.arange(1, universe + 1) ** 0.8  # Zipf-like: a few hub genes appear in many sets
    popularity /= popularity.sum()
    sizes = np.clip(rng.lognormal(mean=3.5, sigma=0.9, size=n_terms).astype(int), 10, 500)
    terms: dict[str, TermRecord] = {}
    base: list[set[str]] = []
    for index, size in enumerate(sizes):
        if base and rng.random() < NEAR_DUPLICATE_RATE:  # near-duplicates of an earlier term
            genes = set(base[rng.integers(len(base))])
            drop = rng.integers(0, max(1, len(genes) // 10) + 1)
            for gene in list(genes)[:drop]:
                genes.discard(gene)
        else:
            genes = {f"G{g}" for g in rng.choice(universe, size=size, replace=False, p=popularity)}
            base.append(genes)
        term_id = f"T{index:06d}"
        terms[term_id] = TermRecord(term_id, "", {"family": "f", "source_tag": "S", "display_name": term_id}, genes)
    return terms


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--terms", type=int, default=50_000)
    parser.add_argument("--seed", type=int, default=1)
    parser.add_argument("--cutoff", type=float, default=0.85)
    args = parser.parse_args(argv)
    terms = synthetic_terms(args.terms, args.seed)
    start = time.perf_counter()
    result = remove_redundancy(terms, "bench", "f", RedundancyConfig(jaccard_cutoff=args.cutoff), {})
    elapsed = time.perf_counter() - start
    peak_kb = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    peak_gb = peak_kb / (1024**3 if sys.platform == "darwin" else 1024**2)
    print(
        f"terms={args.terms} kept={len(result.kept)} removed={len(result.removed)} "
        f"seconds={elapsed:.1f} peak_rss_gb={peak_gb:.2f}"
    )
    within = elapsed < TIME_BUDGET_SECONDS and peak_gb < MEMORY_BUDGET_GB
    print("within budget" if within else "OVER BUDGET")
    return 0 if within else 1


if __name__ == "__main__":
    raise SystemExit(main())
