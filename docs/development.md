# Development

Install ALIEN with development dependencies from a local checkout:

```bash
git clone https://github.com/deminden/alien.git
cd alien
python -m pip install -e ".[dev]"
```

Run the test suite:

```bash
python -m pytest -q
```

Build source and wheel distributions:

```bash
python -m build
```

Check package metadata before publishing:

```bash
python -m twine check dist/*
```

The GitHub Actions workflow in `.github/workflows/ci.yml` runs two checks on pushes and pull requests:

- `Tests`: installs ALIEN with development extras and runs `python -m pytest -q`.
- `Package build check`: builds the source distribution and wheel, runs `twine check`, installs the built wheel, and verifies the CLI/import smoke path.

## Contributing

Contributions are very welcome. If you would like to improve ALIEN, open an issue to discuss bugs, data-source behavior, target namespace needs, or design changes before large pull requests.

Pull requests are especially useful for:

- Additional source adapters for commonly used gene-set libraries.
- New target adapters beyond Ensembl-style GTFs, such as Entrez or UniProt output.
- Better curated defaults for filtering lists, broad-term patterns, non-gene tokens, and source priorities.
- Tests against small deterministic fixtures and established external resources.
- Mapping-audit improvements for ambiguous symbols, retired IDs, and source-provided identifiers.
- Documentation, worked examples, and reproducibility notes.
- Performance improvements for large source collections.

Please keep generated data, downloaded caches, and large release artifacts out of git. Tiny fixtures are preferred for tests.

## Future Plans

Near-term directions:

- Broader target adapter system for non-Ensembl namespaces.
- More source readers with explicit provenance and licensing notes.
- Better community-reviewed filtering presets for broad terms, non-gene tokens, redundancy behavior, and source priority by term family.
- Clearer cache management commands for managed resources.
- Richer QC summaries for mapping, redundancy, and source coverage.
- More cross-checks against upstream MSigDB/msigdbr behavior.

Longer-term directions:

- Support for non-human organisms where mapping resources are sufficiently auditable.
- Curated Enrichr source presets for common pathway, disease, and perturbation libraries.
- More stable lower-level APIs once the package has real-world usage feedback.
