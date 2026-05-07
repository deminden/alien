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
