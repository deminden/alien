# ALIEN

**ALIEN: Audited Library Integration for External Namespaces** builds namespace-specific GMT libraries for human gene-set workflows.

ALIEN 0.1.0 is centered on one job: take local/cache source libraries, normalize their gene memberships into a canonical table, project them into configured target namespaces, and write combined GMT files with audit metadata.

## Install

```bash
python -m pip install -e ".[dev]"
```

## Quick Start

Prepare local source caches and a target annotation, then run:

```bash
alien build --config configs/production.yml --outdir data/alien_gmt --workers 4
```

The primary outputs are:

- `gmt/<target_namespace>.gmt`: combined GMT for each namespace, including `symbols.gmt` when enabled.
- `metadata/`: term manifest, gene mapping tables, removed terms, unmapped genes, ambiguity logs, and provenance.
- `qc/`: collection, mapping, redundancy, target-universe, and warning summaries.

## Source Inputs

ALIEN 0.1.0 reads local/cache sources. It supports:

- `msigdb_cache`: a directory of `msigdbr_<SOURCE_TAG>.tsv.gz` files.
- `msigdb_tsv`: a MSigDB-like TSV with term names, gene symbols, and optional source Ensembl IDs.
- `symbol_gmt`: a GMT file whose members are gene symbols.
- `canonical_tsv`: ALIEN’s normalized membership schema.

Canonical tables require:

```text
term_id    gene_symbol
```

Recommended optional fields are `original_name`, `display_name`, `description`, `source`, `source_tag`, `collection`, `family`, `aspect`, `gene_id`, and `gene_id_namespace`. If `gene_id_namespace` identifies Ensembl IDs, ALIEN audits those IDs before falling back through symbol mapping.

## Python API

```python
from alien import build

result = build("configs/production.yml", outdir="data/alien_gmt", workers=4)
print(result.namespaces)
```

Only the small facade API is stable for 0.1.0. Lower-level modules are importable for experimentation but may change while the package grows.

## MSigDB Cache Preparation

ALIEN 0.1.0 does not download MSigDB from Python. Use existing local caches or the optional helper:

```bash
Rscript scripts/fetch_msigdb.R data/alien_sources/msigdb FALSE 26.1.0
```

A Python-based MSigDB downloader is planned for 0.1.1.

## Scope

The 0.1.0 release officially supports human gene sets using HGNC symbols and Ensembl target namespaces. The code is organized so broader namespace integrations can be added later without tying the package to any single downstream analysis project.
