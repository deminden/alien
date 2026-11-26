# Development

Install ALIEN with development dependencies from a local checkout and enable the pre-commit hooks:

```bash
git clone https://github.com/deminden/alien.git
cd alien
python -m pip install -e ".[dev]"
pre-commit install
```

## Checks

```bash
ruff check src tests scripts benchmarks     # lint (complexity <= 12, at most 6 arguments)
ruff format --check src tests scripts benchmarks
mypy                                         # strict, configured in pyproject.toml
python -m pytest -q --cov=alien              # coverage must stay at or above 90%
python -m pytest -q -m network               # optional: real downloads (run nightly in CI)
python benchmarks/redundancy_benchmark.py --terms 50000
```

The suite blocks real network connections; tests serve HTTP from the `fake_http` fixture in `tests/conftest.py`.

## Golden outputs

`tests/golden/` builds an offline fixture covering every mapping status and compares every output file with the files under `tests/golden/expected/`. When a change is meant to alter outputs, regenerate them and explain the difference in the commit message:

```bash
ALIEN_UPDATE_GOLDEN=1 python -m pytest -q tests/golden
git diff tests/golden/expected
```

## Generated files

- `README-PYPI.md` is generated from `README.md` with links pinned to the release tag: `python scripts/generate_pypi_readme.py` (CI runs it with `--check`; the release workflow regenerates it for the tag).
- `docs/outputs.md` is generated from `src/alien/schema.py`: `python scripts/generate_outputs_doc.py` (a test fails when it is stale).
- `docs/config.md` is generated from the configuration models in `src/alien/config/`: `python scripts/generate_config_doc.py`. Every option needs an attribute docstring (a test fails otherwise), and `docs/config.schema.json` is refreshed with `alien config schema > docs/config.schema.json`.
- `docs/config.schema.json` is the configuration JSON Schema: `alien config schema > docs/config.schema.json`.

## Continuous integration

`.github/workflows/ci.yml` runs lint and type checks, the suite on Python 3.12 and 3.13 (Ubuntu) and 3.12 (macOS), the suite against the lowest supported dependency versions (`uv pip install --resolution lowest-direct`), a non-blocking redundancy benchmark, a nightly network integration test, and a package build check. `.github/workflows/release.yml` publishes tagged releases to PyPI with trusted publishing, attaches build provenance attestations, and creates a GitHub Release from the matching `CHANGELOG.md` section.

## Releasing

1. Update `version` in `pyproject.toml` (the only place it is written) and add a `CHANGELOG.md` section for it.
2. Regenerate `README-PYPI.md` for the tag: `python scripts/generate_pypi_readme.py --ref vX.Y.Z`.
3. Tag `vX.Y.Z` and push the tag.

See [CONTRIBUTING.md](../CONTRIBUTING.md) for contribution guidelines and [architecture.md](architecture.md) for how the code is organized.

## Future plans

- Target adapters beyond Ensembl-style GTFs, such as Entrez or UniProt output.
- Species profiles (the human assumptions are confined to three modules).
- Community-reviewed filtering presets for broad terms, non-gene tokens and source priority.
- More cross-checks against upstream MSigDB/msigdbr behaviour.
