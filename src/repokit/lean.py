"""Proofs in Lean: the project's pinning, the stored report of what each declaration rests on, and the
blueprint linking each derived claim to its proposition and its Lean declarations.

    python -m repokit.lean report     rebuild the stored report from the built project (run by `make lean`)

The checks read the stored report, `results/lean.json`, so `make check` needs no Lean. The report records a
hash of the Lean sources it was built from, so a report older than the proofs is caught without building them.
"""
from __future__ import annotations

import hashlib
import json
import pathlib
import re
import subprocess
import sys
import tempfile
import tomllib

from . import paths

STANDARD_AXIOMS = ("propext", "Classical.choice", "Quot.sound")
REPORT = paths.RESULTS / "lean.json"
PINNED_FILES = ("lean-toolchain", "lakefile.toml", "lake-manifest.json")
DECLARATION = re.compile(
    r"^(?:@\[[^\]]*\]\s*)?(?:(?:private|protected|noncomputable)\s+)*"
    r"(theorem|lemma|def|abbrev|structure|class|inductive|instance)\s+([^\s:({\[]+)", re.M)
LABEL = re.compile(r"`((?:prop|lem|cor|thm|def|eq|app):[A-Za-z0-9_.-]+)`")
STATEMENT_KINDS = ("prop:", "lem:", "thm:", "cor:")  # the labels of statements a declaration can prove
TEX_LABEL = re.compile(r"\\label\{([^}]+)\}")
IMPORT = re.compile(r"^import\s+(\S+)", re.M)
# Mathlib's conventions, as far as they can be read from the source: a theorem is named in snake case
# whose words are lower camel case (`weighted_isHermitian`); a definition carries no underscore.
THEOREM_NAME = re.compile(r"^[a-z][A-Za-z0-9']*(?:_[a-zA-Z0-9][A-Za-z0-9']*)*$")
DEFINITION_NAME = re.compile(r"^[A-Za-z][A-Za-z0-9']*$")
TOOLCHAIN = re.compile(r"^leanprover/lean4:(?:v\d+\.\d+\.\d+(?:-rc\d+)?|nightly-\d{4}-\d{2}-\d{2})$")
COMMIT = re.compile(r"^[0-9a-f]{40}$")


def project() -> pathlib.Path | None:
    """The Lean project's folder, named by `lean` in `[tool.repokit]`; None when the work has none."""
    name = paths.config().get("lean")
    return paths.ROOT / name if name else None


def library(root: pathlib.Path) -> str:
    """The name of the project's Lean library, from its lakefile."""
    lake = tomllib.loads((root / "lakefile.toml").read_text(encoding="utf-8"))
    return lake["lean_lib"][0]["name"]


def lean_files(root: pathlib.Path) -> list[pathlib.Path]:
    return sorted(p for p in root.rglob("*.lean") if ".lake" not in p.relative_to(root).parts)


def module(root: pathlib.Path, path: pathlib.Path) -> str:
    """The module a Lean file defines: its path inside the project, dotted."""
    return ".".join(path.relative_to(root).with_suffix("").parts)


def sources_hash(root: pathlib.Path) -> str:
    """A hash of everything the proofs are built from: the Lean files and the files that pin their versions."""
    digest = hashlib.sha256()
    for path in lean_files(root) + [root / name for name in PINNED_FILES if (root / name).is_file()]:
        digest.update(path.relative_to(root).as_posix().encode() + b"\0" + path.read_bytes().replace(b"\r\n", b"\n"))
    return digest.hexdigest()


def declarations(path: pathlib.Path) -> list[dict]:
    """The declarations a Lean file states, each with its kind, name, docstring and line."""
    text = path.read_text(encoding="utf-8")
    found = []
    for match in DECLARATION.finditer(text):
        before = text[:match.start()].rstrip()
        doc = before[before.rfind("/--") + 3:-2].strip() if before.endswith("-/") and "/--" in before else None
        found.append({"kind": match.group(1), "name": match.group(2), "doc": doc,
                      "line": text.count("\n", 0, match.start()) + 1})
    return found


def references(claims: dict) -> list[dict]:
    """Every Lean declaration a claim names, as evidence or as owed."""
    out = []
    for key, claim in claims.items():
        for field in ("evidence", "owed"):
            for ref in evidence_list(claim, field):
                if ref.startswith("lean:"):
                    path, _, name = ref[len("lean:"):].partition("::")
                    out.append({"claim": key, "field": field, "path": path, "name": name})
    return out


def reported(report: dict, root: pathlib.Path, ref: dict) -> dict | None:
    """The report's entry for a referenced declaration, matched by its module and the end of its name."""
    wanted = module(root, paths.ROOT / ref["path"])
    for name, entry in (report.get("declarations") or {}).items():
        if entry.get("module") == wanted and (name == ref["name"] or name.endswith("." + ref["name"])):
            return entry
    return None


def is_proved(entry: dict | None) -> bool:
    return entry is not None and set(entry.get("axioms", [])) <= set(STANDARD_AXIOMS)


def load_report() -> dict:
    return json.loads(REPORT.read_text(encoding="utf-8")) if REPORT.is_file() else {}


def paper_labels(paper: str) -> set[str]:
    return set(TEX_LABEL.findall(paper))


REPORT_SCRIPT = """import {library}
open Lean Elab Command

#eval show CommandElabM Unit from do
  let env ← getEnv
  for (name, info) in env.constants.map₁.toList do
    if name.isInternal then continue
    match env.getModuleIdxFor? name with
    | none => pure ()
    | some idx =>
      let mod := env.header.moduleNames[idx.toNat]!
      if (`{library}).isPrefixOf mod then
        let kind := match info with
          | .thmInfo _ => "theorem"
          | .defnInfo _ => "def"
          | _ => "other"
        let axioms ← collectAxioms name
        IO.println s!"{{name}}\\t{{mod}}\\t{{kind}}\\t{{axioms.toList}}"
"""


def build_report(root: pathlib.Path) -> dict:
    """Ask Lean, through the built project, what every declaration of the library rests on."""
    lib = library(root)
    with tempfile.TemporaryDirectory() as tmp:
        script = pathlib.Path(tmp) / "Report.lean"
        script.write_text(REPORT_SCRIPT.format(library=lib), encoding="utf-8", newline="\n")
        run = subprocess.run(["lake", "env", "lean", str(script)], cwd=root, capture_output=True, text=True)
    if run.returncode != 0:
        raise SystemExit(f"Lean could not report on the library:\n{run.stdout}{run.stderr}")
    stated = {d["name"] for p in lean_files(root) for d in declarations(p)}
    found = {}
    for line in run.stdout.splitlines():
        parts = line.split("\t")
        # a declaration may be stated under a dotted name (`HodgeSplit.gaugeFixedAction_one`): match any suffix
        words = parts[0].split(".")
        if len(parts) != 4 or not any(".".join(words[i:]) in stated for i in range(len(words))):
            continue
        axioms = [a.strip() for a in parts[3].strip("[]").split(",") if a.strip()]
        found[parts[0]] = {"module": parts[1], "kind": parts[2], "axioms": sorted(axioms)}
    toolchain = (root / "lean-toolchain").read_text(encoding="utf-8").strip()
    manifest = json.loads((root / "lake-manifest.json").read_text(encoding="utf-8"))
    return {"toolchain": toolchain, "packages": {p["name"]: p["rev"] for p in manifest.get("packages", [])},
            "sources": sources_hash(root), "declarations": dict(sorted(found.items()))}


def blueprint_md(claims: dict, paper: str) -> str:
    """The blueprint: each derived claim with the proposition proving it, its Lean declarations and whether they are proved."""
    root = project()
    report = load_report()
    docs = {}
    if root is not None and root.is_dir():
        docs = {(module(root, p), d["name"]): d["doc"] or "" for p in lean_files(root) for d in declarations(p)}
    rows = []
    for key, claim in claims.items():
        if claim.get("standing") != "derived":
            continue
        refs = [r for r in references({key: claim})]
        cells, labels, states = [], set(), []
        for ref in refs:
            doc = docs.get((module(root, paths.ROOT / ref["path"]), ref["name"]), "") if root else ""
            labels |= {label for label in LABEL.findall(doc) if label.startswith(STATEMENT_KINDS)}
            proved = ref["field"] == "evidence" and is_proved(reported(report, root, ref)) if root else False
            states.append("proved" if proved else "owed")
            cells.append(f"`{ref['name']}` ({ref['path']})")
        status = "not formalized" if not refs else ("proved" if all(s == "proved" for s in states) else "owed")
        if not refs and claim.get("unformalized"):
            status += f": {claim['unformalized']}"
        proposition = ", ".join(sorted(labels)) or str(claim.get("proof") or "")
        rows.append(f"| {key} | {proposition} | {'; '.join(cells) or '—'} | {status} |")
    if not rows:
        return "No claim is derived."
    return "| Claim | Proved in the paper by | Lean declaration | In Lean |\n| --- | --- | --- | --- |\n" + "\n".join(rows)


def main(argv: list[str] | None = None) -> int:
    args = sys.argv[1:] if argv is None else argv
    root = project()
    if args != ["report"] or root is None:
        print("usage: python -m repokit.lean report   (the work names its Lean project as `lean` in [tool.repokit])")
        return 2
    report = build_report(root)
    REPORT.parent.mkdir(parents=True, exist_ok=True)
    REPORT.write_text(json.dumps(report, indent=2, ensure_ascii=False) + "\n", encoding="utf-8", newline="\n")
    proved = sum(1 for e in report["declarations"].values() if is_proved(e))
    print(f"wrote {paths.rel(REPORT)}: {len(report['declarations'])} declarations, {proved} resting only on the standard axioms")
    return 0


if __name__ == "__main__":
    sys.exit(main())


def evidence_list(claim: dict, field: str = "evidence") -> list[str]:
    value = claim.get(field) or []
    return [value] if isinstance(value, str) else [str(v) for v in value]
