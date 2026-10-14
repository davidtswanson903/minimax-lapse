"""Data from the world: published datasets brought into the paper with their provenance, their
role, the premises they add, and the equations that carry them into the paper's quantities.

`data/datasets.yaml` lists each dataset, and its file sits under `data/external/` with its
checksum. A dataset is fixed by comparison: imported, cited, a property of its source, never a
claim. A placeholder is a computed stand-in with the schema of the published data it stands for:
it may be committed so the whole route from file to printed number can run, it is marked on every
surface that uses it, and no release is made while one remains.

Roles:
    consistency     a derived inequality held against measurement
    estimate        a quantity the paper defines, estimated from measurement
    correspondence  a measured structure with a reading's form, its closeness stated
    witness         something the paper's object forbids, observed where it is observed
"""
from __future__ import annotations

import csv
import hashlib
import pathlib

import yaml

from . import paths
from .tables import Sym, Table

ROLES = ("consistency", "estimate", "correspondence", "witness")


def registry() -> list[dict]:
    path = paths.DATA / "datasets.yaml"
    return list(yaml.safe_load(path.read_text(encoding="utf-8")) or []) if path.exists() else []


def entry(key: str) -> dict:
    for item in registry():
        if item["key"] == key:
            return item
    raise KeyError(f"no dataset {key!r} in data/datasets.yaml")


def checksum(path: pathlib.Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load(key: str) -> list[dict]:
    """A dataset's rows, after its file is checked against its recorded checksum."""
    item = entry(key)
    path = paths.ROOT / item["file"]
    if checksum(path) != item["sha256"]:
        raise ValueError(f"{key}: {item['file']} does not match its recorded checksum")
    with path.open(encoding="utf-8", newline="") as handle:
        return [{k: _number(v) for k, v in row.items()} for row in csv.DictReader(handle)]


def _number(value: str):
    try:
        return float(value)
    except ValueError:
        return value


def placeholders() -> list[str]:
    return [item["key"] for item in registry() if item.get("placeholder")]


def table() -> Table | None:
    items = registry()
    if not items:
        return None
    rows = []
    for item in items:
        source = item.get("source") or item.get("stands_in_for_ref")
        stands = item.get("title", item["key"])
        cite = [Sym(f" ({source})", r"~\cite{" + source + "}")] if source else []
        standing = "placeholder" if item.get("placeholder") else "fixed by comparison"
        bears = ", ".join(item.get("bears_on", []))
        rows.append(([stands] + cite, ", ".join(item.get("role", [])), bears, standing))
    return Table(name="datasets", columns=("Data", "Role", "Bears on", "Standing"),
                 align=(r">{\raggedright\arraybackslash}p{2.7in}", "l", r">{\raggedright\arraybackslash}p{2.1in}", "l"),
                 rows=tuple(rows))
