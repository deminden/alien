# Contributing to ALIEN

Contributions are welcome: bug reports, documentation, tests, source adapters, target adapters, mapping-audit improvements and curated defaults. Please open an issue to discuss larger changes before a pull request.

## Ground rules

- **Outputs are a contract.** The golden tests in `tests/golden/` compare every output byte for byte. If your change alters outputs on purpose, regenerate them (`ALIEN_UPDATE_GOLDEN=1 python -m pytest tests/golden`) in the same commit and explain why in the commit message. Changes to output files or columns bump `OUTPUT_SCHEMA_VERSION` in `src/alien/schema.py` and need a `CHANGELOG.md` entry.
- **Bugs get a failing test first.** Add a regression test that fails before the fix.
- **One concern per commit**, each leaving the suite green. Do not mix file moves with logic changes.
- **Keep it small.** Modules up to about 400 lines, functions up to about 60 lines with complexity at most 12 and at most 6 parameters (ruff enforces the last two). Reuse the primitives in [docs/architecture.md](docs/architecture.md) instead of writing a second HTTP client, table writer or ID normalizer.
- **No new runtime dependency** without a sentence of justification in the pull request.
- **Scientific defaults** (mapping order, size limits, Jaccard cutoff, representative choice) change only after discussion in an issue.
- Keep downloaded data, caches and large artifacts out of git; tiny fixtures are preferred for tests.

## Setup and checks

See [docs/development.md](docs/development.md). Before opening a pull request:

```bash
pre-commit run --all-files
python -m pytest -q --cov=alien
```
