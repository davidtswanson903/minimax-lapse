"""Shared base for a paper repository: named values and their one formatter,
generation of every derived file from one source, golden checks for experiments,
and the checks that keep the published surfaces in agreement."""
from __future__ import annotations

import pathlib
import re
import sys

__version__ = "1.0.0"

PYPROJECT = pathlib.Path(__file__).resolve().parents[2] / "pyproject.toml"


def required_python(pyproject: pathlib.Path = PYPROJECT) -> tuple[int, int] | None:
    """The least Python version `requires-python` names, read without tomllib so an older interpreter can still say what it lacks."""
    try:
        text = pyproject.read_text(encoding="utf-8")
    except FileNotFoundError:
        return None  # an installed copy, away from its pyproject.toml
    found = re.search(r'^requires-python\s*=\s*">=\s*(\d+)\.(\d+)', text, re.M)
    return (int(found.group(1)), int(found.group(2))) if found else None


def require_python(version: tuple = tuple(sys.version_info), pyproject: pathlib.Path = PYPROJECT) -> None:
    """Refuse an interpreter older than the one the work requires, naming both."""
    least = required_python(pyproject)
    if least and tuple(version[:2]) < least:
        raise SystemExit(f"this repository needs Python {least[0]}.{least[1]} or later (requires-python in pyproject.toml); "
                         f"this is Python {version[0]}.{version[1]}")


require_python()
