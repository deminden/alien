# Architecture

ALIEN is a pipeline of stages over a small set of shared primitives. Each module has one job.

## Stages (`alien.pipeline`)

1. **Configure** (`alien.config`): defaults on frozen pydantic models; layers merged from the `extends:` chain, the config file, overrides and explicit options; relative paths resolved against the declaring file; unknown keys rejected with their YAML path.
2. **Load sources** (`alien.sources`): every enabled source spec goes to the loader registered for its `type` and returns the canonical membership table (`alien.schema.CANONICAL_COLUMNS`). Term-ID collisions are checked here and, by default, renamed apart.
3. **Load resources** (`alien.resources`): the HGNC complete set, NCBI Gene rescue maps, prepared targets (annotation GTF, metadata fallbacks, output-gene lists) and the Ensembl archive resolver.
4. **Map** (`alien.mapping`): the `symbols` namespace, then each target. For a target, every distinct key `(symbol, source Ensembl ID, NCBI allowed)` is resolved once by `resolver.resolve`, a pure ordered chain of steps; outcomes are joined back to memberships to build term gene sets, audit tables and status counts.
5. **Filter** (`alien.filtering`): size and broad-term filters per namespace, then redundancy removal per (namespace, family) group, serially or in a process pool.
6. **Report** (`alien.reports`): term and source manifests, QC summaries and provenance.
7. **Publish** (`alien.outputs`): the build is written into a hidden staging directory, pruned to the output mode, given a `MANIFEST.json`, and swapped into place; on failure the staging tree is kept under `<outdir>/failed/`.

`fetch()` runs stages 1 to 3 plus the archive lookups and writes the lock file instead of building.

## Shared primitives

| Primitive | Module | Used for |
|---|---|---|
| `fetch`, `get_json`, `post_json`, `FetchRecord` | `alien.io.http` | Every HTTP call: one session with retries that honour Retry-After, a bioalien User-Agent, streaming to a unique temporary file, hashing while streaming, `os.replace`. |
| `atomic_write_text/bytes/json`, `atomic_dir_swap` | `alien.io.atomic` | Every file ALIEN writes, and the output-directory swap. |
| `write_table`, `read_table`, `with_columns` | `alien.io.tables` | Every TSV / TSV.GZ, with declared column order and deterministic gzip. |
| `parse_gmt`, `write_gmt` / `parse_gtf_genes` | `alien.io.gmt` / `alien.io.gtf` | GMT and GTF formats. |
| Column lists, `OUTPUT_FILES`, schema versions | `alien.schema` | One declaration per table; `docs/outputs.md` is generated from it. |
| Configuration options, defaults and their descriptions | `alien.config` models (attribute docstrings) | `docs/config.md` and `docs/config.schema.json` are generated from them. |
| `MappingStatus`, `RemovalReason`, `TermRecord`, `Target`, `ensembl_stable_id`, `sanitize_term_id` | `alien.model` | Domain types and the single implementation of ID normalization. `Target.effective_ids` is the one rule for which IDs a target may emit. |
| `ResourceCache`, `LockFile`, `derived_cache_path` | `alien.cache` | Cache layout, per-resource file locks, offline and frozen modes, the lock file, and cache keys for derived data. |
| `register`, `loader_for` | `alien.sources.registry` | Source types, including plugins from the `alien.sources` entry-point group. |
| `AlienError` and subclasses | `alien.errors` | Errors the CLI maps to exit codes. |

## Human-specific code

Three places assume human data: the `ENSG` prefix in `alien.model.ensembl_stable_id`, `alien.resources.hgnc`, and the taxon filter in `alien.resources.ncbi`. A species profile would touch only these.

## Tests

`tests/` mirrors the package. `tests/golden/` builds an offline fixture that exercises every mapping status in all output modes with 1 and 4 workers and compares every file byte for byte. `tests/mapping/test_equivalence.py` and `tests/filtering/test_redundancy_equivalence.py` are hypothesis property tests against frozen copies of the 0.1.x algorithms in `tests/legacy/`. `tests/regressions/` holds one test per fixed bug. Real network access fails the suite unless a test is marked `network`.
