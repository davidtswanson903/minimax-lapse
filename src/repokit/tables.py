"""Tables, built by the work's own `tables` module from stored results.

A table is rows of cells; a cell is text, a Value, or a sequence of them, so a
table prints every number through the one formatter, as the text does. Each table
is written as a LaTeX tabular for the paper to place, build/latex/table-<name>.tex,
and as a Markdown block a template may place, {{table_<name>}}.
"""
from __future__ import annotations

import importlib
from dataclasses import dataclass

from . import paths
from .values import Value


@dataclass(frozen=True)
class Sym:
    """Text written differently in TeX: a symbol, a dash, a little mathematics."""
    text: str
    tex: str


@dataclass(frozen=True)
class Table:
    name: str
    columns: tuple   # header cells
    align: tuple     # one LaTeX column spec per column: l, c, r or p{width}
    rows: tuple      # tuples of cells; a row of None draws a rule


def collect(results: dict, values: dict) -> list[Table]:
    try:
        module = importlib.import_module(f"{paths.package()}.tables")
    except ModuleNotFoundError:
        return []
    return list(module.tables(results, values))


# a bar in text is an em-dash in the default font encoding, so it is written as the text bar
_TEX_SPECIAL = {"&": r"\&", "%": r"\%", "#": r"\#", "_": r"\_", "–": "--", "—": "---", "|": r"\textbar{}"}


def _parts(cell) -> tuple:
    return tuple(cell) if isinstance(cell, (list, tuple)) else (cell,)


def _tex(part) -> str:
    if isinstance(part, Value):
        return f"${part.tex()}$"
    if isinstance(part, Sym):
        return part.tex
    return "".join(_TEX_SPECIAL.get(c, c) for c in str(part))


def _text(part) -> str:
    if isinstance(part, Value):
        return part.text()
    if isinstance(part, Sym):
        return part.text
    return str(part)


RAGGED = r">{\raggedright\arraybackslash}"


def tex(table: Table) -> str:
    # A ragged-right paragraph column is written as a plain p column whose cells open with \raggedright, its rows
    # ended by \tabularnewline, so the table needs no array package: under the LaTeX kernel of 2025-11-01, array and
    # revtex4-2 together break every tabular.
    ragged = [a.startswith(RAGGED) for a in table.align]
    align = [a[len(RAGGED):] if r else a for a, r in zip(table.align, ragged)]

    def row(cells):
        if cells is None:
            return r"\midrule"
        return " & ".join((r"\raggedright " if r else "") + "".join(_tex(p) for p in _parts(c))
                          for c, r in zip(cells, ragged)) + r" \tabularnewline"
    # booktabs rules, which the presentation layer loads
    lines = [rf"\begin{{tabular}}{{{''.join(align)}}}", r"\toprule", row(table.columns), r"\midrule"]
    lines += [row(r) for r in table.rows]
    lines += [r"\bottomrule", r"\end{tabular}"]
    return "\n".join(lines) + "\n"


def markdown(table: Table) -> str:
    def row(cells):
        return "| " + " | ".join("".join(_text(p) for p in _parts(c)).replace("|", "\\|") for c in cells) + " |"
    rule = "| " + " | ".join({"r": "---:", "c": ":---:"}.get(a[0], "---") for a in table.align) + " |"
    return "\n".join([row(table.columns), rule, *[row(r) for r in table.rows if r is not None]]) + "\n"
