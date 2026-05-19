import subprocess
import sys
from pathlib import Path


def test_pypi_readme_is_generated_from_readme():
    root = Path(__file__).resolve().parents[1]
    result = subprocess.run(
        [sys.executable, "scripts/generate_pypi_readme.py", "--check"],
        cwd=root,
        capture_output=True,
        text=True,
    )

    assert result.returncode == 0, result.stderr


def test_pypi_readme_uses_tagged_github_example_links():
    root = Path(__file__).resolve().parents[1]
    text = (root / "README-PYPI.md").read_text(encoding="utf-8")

    assert "https://github.com/deminden/alien/tree/v0.1.6/examples/" in text
    assert "https://github.com/deminden/alien/blob/v0.1.6/docs/usage.md" in text
