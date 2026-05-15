#!/usr/bin/env python3
"""Generate the PyPI README from README.md with GitHub tag links."""

from __future__ import annotations

import argparse
import re
import sys
import tomllib
from pathlib import Path
from urllib.parse import quote, urldefrag


REPO_URL = "https://github.com/deminden/alien"
RAW_URL = "https://raw.githubusercontent.com/deminden/alien"
ROOT = Path(__file__).resolve().parents[1]
README = ROOT / "README.md"
PYPI_README = ROOT / "README-PYPI.md"
PYPROJECT = ROOT / "pyproject.toml"
LINK_RE = re.compile(r"(!?)\[([^\]]*)\]\(([^)\s]+)(?:\s+\"[^\"]*\")?\)")


def project_version() -> str:
    metadata = tomllib.loads(PYPROJECT.read_text(encoding="utf-8"))
    return metadata["project"]["version"]


def should_rewrite(target: str) -> bool:
    if target.startswith(("#", "/", "http://", "https://", "mailto:", "tel:")):
        return False
    return True


def github_url(target: str, ref: str, *, raw: bool) -> str:
    path_part, fragment = urldefrag(target)
    path = Path(path_part)
    if path.is_absolute() or ".." in path.parts:
        return target

    encoded_path = quote(path_part, safe="/")
    if raw:
        url = f"{RAW_URL}/{quote(ref, safe='')}/{encoded_path}"
    else:
        route = "tree" if (ROOT / path_part).is_dir() or path_part.endswith("/") else "blob"
        url = f"{REPO_URL}/{route}/{quote(ref, safe='')}/{encoded_path}"
    if fragment:
        url = f"{url}#{fragment}"
    return url


def rewrite_links(markdown: str, ref: str) -> str:
    blocks = re.split(r"(^```.*?^```)", markdown, flags=re.MULTILINE | re.DOTALL)
    rewritten: list[str] = []

    for block in blocks:
        if block.startswith("```"):
            rewritten.append(block)
            continue

        def replace(match: re.Match[str]) -> str:
            image_marker, text, target = match.groups()
            if not should_rewrite(target):
                return match.group(0)
            url = github_url(target, ref, raw=bool(image_marker))
            return f"{image_marker}[{text}]({url})"

        rewritten.append(LINK_RE.sub(replace, block))

    return "".join(rewritten)


def generated_readme(ref: str) -> str:
    markdown = README.read_text(encoding="utf-8")
    body = rewrite_links(markdown, ref).rstrip()
    return "<!-- Generated from README.md by scripts/generate_pypi_readme.py. Do not edit directly. -->\n\n" + body + "\n"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ref", default=f"v{project_version()}", help="Git ref used in generated GitHub links.")
    parser.add_argument("--check", action="store_true", help="Fail if README-PYPI.md is not up to date.")
    args = parser.parse_args()

    generated = generated_readme(args.ref)
    if args.check:
        current = PYPI_README.read_text(encoding="utf-8") if PYPI_README.exists() else ""
        if current != generated:
            print("README-PYPI.md is out of date. Run scripts/generate_pypi_readme.py.", file=sys.stderr)
            return 1
        return 0

    PYPI_README.write_text(generated, encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
