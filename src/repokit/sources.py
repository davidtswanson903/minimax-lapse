"""Reading the reader's tree: the claims, the stored results, and the values the work's own package computes
from them. The composer reads the registries that drive the paper's text from the private root."""
from __future__ import annotations

import importlib
import json
import pathlib
import re

import yaml

from . import paths
from .values import Value


def _yaml(path: pathlib.Path):
    with open(path, encoding="utf-8") as handle:
        return yaml.safe_load(handle)


def claims() -> dict:
    return _yaml(paths.DATA / "claims.yaml") or {}


def numbers_from() -> dict[str, str]:
    """Which grade's stored results each experiment's printed numbers come from: `numbers_from` under
    [tool.repokit] in pyproject.toml."""
    return {str(k): str(v) for k, v in (paths.config().get("numbers_from") or {}).items()}


def results() -> dict[str, dict]:
    chosen = numbers_from()
    if not chosen:
        return {p.stem: json.loads(p.read_text(encoding="utf-8")) for p in sorted((paths.RESULTS / "paper").glob("*.json"))}
    out = {}
    for name, grade in chosen.items():
        path = paths.RESULTS / grade / f"{name}.json"
        if path.exists():
            out[name] = json.loads(path.read_text(encoding="utf-8"))
    return out


def values() -> dict[str, Value]:
    module = importlib.import_module(f"{paths.package()}.values")
    out: dict[str, Value] = {}
    for value in module.values(results()):
        if value.key in out:
            raise ValueError(f"value {value.key} is defined twice")
        out[value.key] = value
    return out


def text_values(vals: dict[str, Value]) -> dict[str, str]:
    return {key: value.text() for key, value in vals.items()}


def paper() -> str:
    """The paper's source with its own included parts expanded, so every check reads them; generated
    inputs, named by a path, stay as they are."""
    text = paths.PAPER.read_text(encoding="utf-8")
    def expand(match):
        part = paths.PAPER.parent / f"{match.group(1)}.tex"
        return part.read_text(encoding="utf-8") if part.exists() else match.group(0)
    return re.sub(r"\\input\{([A-Za-z0-9_-]+)\}", expand, text)
