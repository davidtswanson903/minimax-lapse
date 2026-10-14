"""Writes the files a reader regenerates from the stored results: the values, the tables the work's own package
builds from its results, the table of the equations beside their code, and the table of the datasets. With --check,
lists the files that differ from what it would write and exits non-zero. The paper's other generated parts are
written by the composer, at the private root, and stay here as files.

    python -m repokit.generate            write
    python -m repokit.generate --check    compare only
"""
from __future__ import annotations

import argparse
import json
import pathlib
import sys

from . import datasets, equations, latex, paths, render, sources, tables

VALUES_JSON = paths.BUILD / "values.json"


def built() -> list[tables.Table]:
    """The tables a reader regenerates: the work's own, the equations' and the datasets'."""
    out = tables.collect(sources.results(), sources.values())
    listed = equations.table(sources.paper())
    if listed is not None:
        out.append(listed)
    held = datasets.table()
    if held is not None:
        out.append(held)
    return out


def outputs() -> dict[pathlib.Path, str]:
    vals = sources.values()
    return {
        VALUES_JSON: json.dumps(
            {key: {"kind": v.kind, "text": v.text(), "tex": v.tex()} for key, v in sorted(vals.items())},
            indent=2, ensure_ascii=False,
        ) + "\n",
        paths.LATEX / "values.tex": latex.values_tex(vals),
        **{paths.LATEX / f"table-{table.name}.tex":
           render.banner(f"src/{paths.package()}/tables.py", "tex") + tables.tex(table) for table in built()},
    }


def stale(out: dict[pathlib.Path, str] | None = None) -> list[str]:
    out = outputs() if out is None else out
    return [paths.rel(p) for p, text in out.items() if not p.exists() or p.read_text(encoding="utf-8") != text]


def write(out: dict[pathlib.Path, str]) -> None:
    """Write each output whose text changed; an unchanged file keeps its time, so what is built from it stays fresh."""
    for path, text in out.items():
        if path.is_file() and path.read_text(encoding="utf-8") == text:
            continue
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(text, encoding="utf-8", newline="\n")


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(prog="python -m repokit.generate")
    parser.add_argument("--check", action="store_true", help="compare with the files on disk; write nothing")
    args = parser.parse_args(argv)
    out = outputs()
    if args.check:
        found = stale(out)
        if found:
            print("generated files differ from the stored results; run `make generate`:")
            for name in found:
                print(f"  {name}")
            return 1
        print(f"{len(out)} generated files match the stored results")
        return 0
    write(out)
    print(f"wrote {len(out)} generated files")
    return 0


if __name__ == "__main__":
    sys.exit(main())
