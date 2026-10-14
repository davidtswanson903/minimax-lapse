"""Equations as code. A function that computes an equation of the paper carries the
equation's label and states the equation in its docstring, so a reader can go from the
equation to the code and back; every such function has a check that reaches the same
result by another route, and says which route in the first line of its docstring.

    @equation("eq:density")
    def densities(...):
        '''D_x = h_x + 1/2 sum of the link pieces at x (eq:density).'''

    @checks("eq:density")
    def test_densities_split_the_hamiltonian():
        '''By summing every piece once ...'''

A labelled equation the code does not compute is marked in the paper by a comment
`% not computed: <why>` on the lines after its label.
"""
from __future__ import annotations

import ast
import importlib
import inspect
import pathlib
import pkgutil
import re
from dataclasses import dataclass

from . import paths
from .tables import Sym, Table


def equation(*labels: str):
    def mark(function):
        function.__equations__ = tuple(labels)
        return function
    return mark


def checks(*labels: str):
    def mark(function):
        function.__checks__ = tuple(labels)
        return function
    return mark


@dataclass(frozen=True)
class Entry:
    label: str
    function: str
    path: str
    line: int
    doc: str


def functions() -> list[Entry]:
    """Every equation function in the work's package, found by importing its modules."""
    package = importlib.import_module(paths.package())
    found = set()
    for info in pkgutil.walk_packages(package.__path__, package.__name__ + "."):
        if info.name.endswith(".figures"):
            continue
        module = importlib.import_module(info.name)
        for obj in vars(module).values():
            labels = getattr(obj, "__equations__", None)
            if labels and getattr(obj, "__module__", None) == module.__name__:
                code = inspect.unwrap(obj).__code__
                for label in labels:
                    found.add(Entry(label, f"{module.__name__}.{obj.__qualname__}",
                                    paths.rel(pathlib.Path(code.co_filename)), code.co_firstlineno,
                                    inspect.getdoc(obj) or ""))
    return sorted(found, key=lambda e: (e.label, e.function))


def check_functions() -> dict[str, list[tuple[str, str]]]:
    """Label -> [(test, the first line of its docstring)], read from the work's own tests."""
    out: dict[str, list[tuple[str, str]]] = {}
    for path in sorted((paths.ROOT / "tests").glob("test_*.py")):
        for node in ast.parse(path.read_text(encoding="utf-8")).body:
            if not isinstance(node, ast.FunctionDef):
                continue
            for decorator in node.decorator_list:
                name = getattr(getattr(decorator, "func", None), "id", None) or getattr(getattr(decorator, "func", None), "attr", None)
                if isinstance(decorator, ast.Call) and name == "checks":
                    doc = (ast.get_docstring(node) or "").split("\n")[0]
                    for arg in decorator.args:
                        if isinstance(arg, ast.Constant):
                            out.setdefault(arg.value, []).append((f"{paths.rel(path)}::{node.name}", doc))
    return out


def paper_labels(paper: str) -> dict[str, bool]:
    """Each labelled equation of the paper, in order, and whether it is marked as not computed."""
    out = {}
    for match in re.finditer(r"\\label\{(eq:[^}]+)\}", paper):
        out[match.group(1)] = bool(re.match(r"\s*(?:\n\s*)?%\s*not computed:", paper[match.end():match.end() + 200]))
    return out


def _tex_code(text: str) -> str:
    return r"\texttt{" + text.replace("_", r"\_\allowbreak{}").replace(".", r".\allowbreak{}") + "}"


def table(paper: str) -> Table | None:
    """Each equation, the function that computes it, and the check that reaches it another way."""
    entries = functions()
    if not entries:
        return None
    order = {label: i for i, label in enumerate(paper_labels(paper))}
    tested = check_functions()
    rows = []
    for label in sorted({e.label for e in entries}, key=lambda l: order.get(l, len(order))):
        here = [e for e in entries if e.label == label]
        code = [Sym(f"[`{e.function}`]({e.path}#L{e.line})", _tex_code(e.function)) for e in here]
        tests = [Sym(f"`{t.split('::')[-1]}`", _tex_code(t.split("::")[-1])) for t, _ in tested.get(label, [])]
        rows.append((Sym(label, r"Eq.~\eqref{" + label + "}"), _joined(code), _joined(tests)))
    return Table(name="equations", columns=("Equation", "Computed by", "Checked by"),
                 align=("l", r">{\raggedright\arraybackslash}p{2.55in}", r">{\raggedright\arraybackslash}p{3.15in}"), rows=tuple(rows))


def _joined(parts: list) -> list:
    out: list = []
    for i, part in enumerate(parts):
        out += ([Sym("; ", "; ")] if i else []) + [part]
    return out or ["none"]
