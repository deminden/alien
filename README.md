# ALIEN

**ALIEN: Audited Library Integration for External Namespaces** builds namespace-specific GMT libraries for human gene-set workflows.

ALIEN is centered on one job: take configured source libraries, normalize their gene memberships into a canonical table, project them into configured target namespaces, and write combined GMT files with audit metadata.

## Install

ALIEN is not on PyPI yet. Install it from the repository:

```bash
git clone https://github.com/deminden/alien.git
cd alien
python -m pip install -e .
```

## Quick Start

Prepare a config with sources and a target annotation, then run:

```bash
alien build --config examples/cancer_dependency.yml --workers 16
```

The primary outputs are:

- `gmt/<target_namespace>.gmt`: combined GMT for each namespace, including `symbols.gmt` when enabled.
- `metadata/`: term manifest, gene mapping tables, removed terms, unmapped genes, ambiguity logs, and provenance.
- `qc/`: collection, mapping, redundancy, target coverage, and warning summaries.

The two filesystem roots are configured in YAML: `project.source_dir` is for downloaded/cached source and mapping resources, while `project.outdir` is the output root containing `gmt/`, `metadata/`, and `qc/`.

## How It Works

ALIEN builds one canonical membership table from configured sources, audits source gene symbols against human mapping resources, then projects each term into the requested output namespaces.

The target namespace is defined in the config. For an Ensembl-style namespace, provide a target name and a GTF annotation:

```yaml
outputs:
  include_symbols: true

targets:
  - name: human_gencode49
    type: ensembl_gtf
    annotation:
      source: GENCODE
      version: "49"
```

This writes `gmt/human_gencode49.gmt`. When `include_symbols` is enabled, ALIEN also writes `gmt/symbols.gmt`.

For human GENCODE releases, any numeric `version` is enough; ALIEN builds the official FTP URL and caches the GTF under `data/alien_sources/gencode/` if it is missing. You can still provide `annotation.path` to pin a local file explicitly.

`source: GENCODE` is not special to the config shape; other Ensembl-style GTF origins can use the same adapter by providing a local path or URL. Fully different target ID systems, such as Entrez or UniProt GMT output, are planned as future target adapters.

Targets can also define `gene_universe` when a dataset file supplies the output ID namespace, or `gene_filter` when the annotation namespace should be intersected with a file/list of allowed Ensembl IDs. In those cases the GTF still acts as the symbol/metadata helper, not necessarily as the owner of the final namespace.

The build keeps audit outputs beside the GMTs so each namespace projection can be traced back to source terms, symbol repairs, unmapped genes, filters, redundancy decisions, and provenance.

`term_id` values must identify one source term unambiguously. If two source libraries reuse the same `term_id` for different term metadata, ALIEN fails by default and writes `metadata/term_id_collisions.tsv` so the IDs can be renamed or prefixed before rebuilding.

Redundancy filtering removes exact duplicate terms and then clusters highly overlapping terms by Jaccard similarity within each namespace and family. The default cutoff is `0.85`:

```yaml
redundancy:
  jaccard_cutoff: 0.85
```

Within each redundant cluster, ALIEN keeps one representative using the configured source priority, then term size and name-based tie-breaks.

Source priority is configured per term family and controls which library wins when redundant terms overlap. For example, pathway terms can prefer `REACTOME` over broader ontology-derived terms:

```yaml
source_priority:
  biology_process_pathway: [REACTOME, WIKIPATHWAYS, KEGG_MEDICUS, GOBP]
```

See [docs/usage.md](docs/usage.md) for the full configuration reference.

## Source Inputs

ALIEN 0.1.3 supports managed MSigDB and Enrichr download/cache sources and local file sources:

- `msigdb_remote`: a Python downloader/reader for the current `msigdbr` Zenodo release cache, with normalized Parquet caches for repeated builds.
- `enrichr_remote`: a Python downloader/reader for Enrichr libraries by public library name.
- `symbol_gmt`: a GMT file whose members are gene symbols.

Additional source metadata fields are documented in [docs/usage.md](docs/usage.md).

## Python API

```python
from alien import build

result = build("examples/pathways.yml", workers=16)
print(result.namespaces)
```

You can also pass the same configuration as a Python dictionary, including target gene universes such as `gene_universe: {"path": "expression.tsv.gz", "column": "feature_id"}` for a known Ensembl ID field.

## MSigDB

```yaml
sources:
  - type: msigdb_remote
    version: "2026.1"
    db_species: HS
    collection: C2
```

This downloads and caches the configured `msigdbr` release archive under `data/alien_sources/msigdb_remote/`, verifies it by MD5, and converts matching RDS files into ALIEN's normalized Parquet cache for later builds.

Remote caches are reused by default. ALIEN also caches normalized MSigDB memberships and prepared NCBI rescue maps, so repeated large builds avoid expensive source-format conversion. Use `alien build --force-download` or per-source `force: true` to refresh managed downloads such as MSigDB, Enrichr, GENCODE, HGNC, and NCBI resources.

## Scope

The 0.1.3 release officially supports human gene sets using HGNC symbols, Python MSigDB and Enrichr cache integration, optimized repeated-build caches, and Ensembl-style target namespaces. The code is organized so broader namespace integrations can be added later without tying the package to any single downstream analysis project.

## Contributing

Contributions are welcome. Useful areas include additional tests, documentation, source adapters, target namespace adapters, mapping-audit improvements, curated filtering/source-priority defaults, and validation against established gene-set resources. See [docs/development.md](docs/development.md) for development setup and current future plans.
