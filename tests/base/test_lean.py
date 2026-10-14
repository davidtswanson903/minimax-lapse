"""Proofs in Lean: pinned, reported, proved without `sorry` or extra axioms, named for their propositions,
imported in the order the claims depend on each other, named as Mathlib names, and linked in the blueprint.
A work without a Lean project has nothing to check."""
import json

import pytest

from repokit import lean, paths

ROOT = lean.project()
pytestmark = pytest.mark.skipif(ROOT is None or not ROOT.is_dir(), reason="the work has no Lean project")


def test_lean_is_pinned():
    toolchain = (ROOT / "lean-toolchain").read_text(encoding="utf-8").strip()
    assert lean.TOOLCHAIN.match(toolchain), f"lean-toolchain names {toolchain}, not one release"
    lake = lean.tomllib.loads((ROOT / "lakefile.toml").read_text(encoding="utf-8"))
    loose = [r["name"] for r in lake.get("require", []) if r.get("rev") in (None, "master", "main")]
    assert not loose, f"requirements without a pinned revision: {loose}"
    manifest = json.loads((ROOT / "lake-manifest.json").read_text(encoding="utf-8"))
    unpinned = [p["name"] for p in manifest.get("packages", []) if not lean.COMMIT.match(p.get("rev", ""))]
    assert manifest.get("packages") and not unpinned, f"lake-manifest.json pins no commit for {unpinned}"


def test_lean_report_is_current():
    report = lean.load_report()
    assert report, "no stored report; run `make lean`"
    assert report.get("sources") == lean.sources_hash(ROOT), "the Lean sources changed since the report; run `make lean`"


def test_cited_declarations_are_proved(claims):
    report = lean.load_report()
    unproved = [f"{r['claim']}: {r['name']}" for r in lean.references(claims)
                if r["field"] == "evidence" and not lean.is_proved(lean.reported(report, ROOT, r))]
    assert not unproved, f"evidence that is not proved without sorry or extra axioms (list it as owed): {unproved}"


def test_cited_declarations_name_their_proposition(claims, paper):
    labels = lean.paper_labels(paper)
    docs = {(lean.module(ROOT, p), d["name"]): d["doc"] for p in lean.lean_files(ROOT) for d in lean.declarations(p)}
    problems = []
    for ref in lean.references(claims):
        doc = docs.get((lean.module(ROOT, paths.ROOT / ref["path"]), ref["name"]))
        named = set(lean.LABEL.findall(doc or ""))
        statements = {label for label in named if label.startswith(lean.STATEMENT_KINDS)}
        if not statements or not named <= labels:
            problems.append(f"{ref['name']}: its docstring names {sorted(named) or 'no label'}; it must name the label of "
                            "a proposition, lemma, theorem or corollary, and every label it names must be the paper's")
    assert not problems, problems


def test_every_lean_theorem_has_a_docstring():
    bare = [f"{lean.module(ROOT, p)}.{d['name']}" for p in lean.lean_files(ROOT) for d in lean.declarations(p)
            if d["kind"] in ("theorem", "lemma", "def", "abbrev") and not d["doc"]]
    assert not bare, f"declarations without a docstring: {bare}"


def test_lean_names_follow_mathlib():
    wrong = []
    for path in lean.lean_files(ROOT):
        for d in lean.declarations(path):
            pattern = lean.THEOREM_NAME if d["kind"] in ("theorem", "lemma") else lean.DEFINITION_NAME
            if d["kind"] != "instance" and not pattern.match(d["name"].rsplit(".", 1)[-1]):
                wrong.append(f"{lean.module(ROOT, path)}: {d['kind']} {d['name']}")
    assert not wrong, wrong


def test_lean_imports_follow_the_dependencies(claims):
    """A module proving claims imports another that proves claims only when its claims use them; a module proving
    none imports none that does. The library's root module, which only gathers the others, is exempt."""
    proved_in: dict[str, set[str]] = {}
    for ref in lean.references(claims):
        proved_in.setdefault(lean.module(ROOT, paths.ROOT / ref["path"]), set()).add(ref["claim"])

    def uses(key, seen=None):
        seen = set() if seen is None else seen
        for used in claims.get(key, {}).get("uses") or []:
            if used not in seen:
                seen.add(used)
                uses(used, seen)
        return seen

    wrong = []
    for path in lean.lean_files(ROOT):
        here = lean.module(ROOT, path)
        if here == lean.library(ROOT):
            continue
        reach = set().union(*(uses(k) | {k} for k in proved_in.get(here, set()))) if here in proved_in else set()
        for imported in lean.IMPORT.findall(path.read_text(encoding="utf-8")):
            theirs = proved_in.get(imported, set())
            if theirs and not theirs <= reach:
                wrong.append(f"{here} imports {imported}, whose claims {sorted(theirs - reach)} its own do not use")
    assert not wrong, wrong


def test_every_unformalized_claim_gives_its_reason(claims):
    """By reading the registry: a derived claim with no Lean declaration says why, and one with a declaration does not."""
    named = {r["claim"] for r in lean.references(claims)}
    silent = [k for k, c in claims.items() if c.get("standing") == "derived" and k not in named
              and not str(c.get("unformalized") or "").strip()]
    stale = [k for k, c in claims.items() if k in named and c.get("unformalized")]
    assert not silent, f"derived claims not proved in Lean that give no reason (`unformalized`): {silent}"
    assert not stale, f"claims proved in Lean that still give a reason for not being: {stale}"


def test_blueprint_lists_every_derived_claim(claims):
    blueprint = (paths.ROOT / "docs" / "blueprint.md").read_text(encoding="utf-8")
    missing = [k for k, c in claims.items() if c.get("standing") == "derived" and f"| {k} |" not in blueprint]
    assert not missing, f"derived claims missing from the blueprint: {missing}"
