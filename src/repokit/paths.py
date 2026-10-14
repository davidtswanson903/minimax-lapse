"""Where things live. The root is found from the working directory, so the same
code serves any checkout of the repository."""
from __future__ import annotations

import os
import pathlib
import tomllib


def _find_root() -> pathlib.Path:
    env = os.environ.get("REPOKIT_ROOT")
    if env:
        return pathlib.Path(env).resolve()
    here = pathlib.Path.cwd().resolve()
    for candidate in (here, *here.parents):
        pyproject = candidate / "pyproject.toml"
        if pyproject.is_file() and "[tool.repokit]" in pyproject.read_text(encoding="utf-8"):
            return candidate
    return pathlib.Path(__file__).resolve().parents[2]


ROOT = _find_root()
DATA = ROOT / "data"
RESULTS = ROOT / "results"
GOLDEN = RESULTS / "golden"
BUILD = ROOT / "build"
LATEX = BUILD / "latex"
DOCS_SRC = ROOT / "docs-src"
DOCS = ROOT / "docs"
PAPER = ROOT / "paper" / "draft" / "paper.tex"
TESTS = ROOT / "tests"
PREVIEW_CUT = 190  # characters a link preview keeps of the description
ARXIV_LIMIT = 1920  # characters arXiv accepts in its abstract field


def config() -> dict:
    text = (ROOT / "pyproject.toml").read_text(encoding="utf-8")
    return tomllib.loads(text).get("tool", {}).get("repokit", {})


def package() -> str:
    return config().get("package", "work")


def project_version() -> str:
    text = (ROOT / "pyproject.toml").read_text(encoding="utf-8")
    return tomllib.loads(text)["project"]["version"]


def rel(path: pathlib.Path) -> str:
    return path.resolve().relative_to(ROOT).as_posix()
