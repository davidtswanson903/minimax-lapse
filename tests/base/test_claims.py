"""Every claim key in the paper resolves, and every claim names evidence that exists:
a test, or a Lean declaration. Evidence still owed is listed as owed, not as evidence."""
import ast
import re

from repokit import latex, paths
from repokit.lean import evidence_list

STATUSES = {"checked", "stated", "open"}
STANDINGS = {"derived", "stated", "measured", "fixed by comparison"}


def _exists(reference: str) -> bool:
    if reference.startswith("lean:"):
        path, _, declaration = reference[len("lean:"):].partition("::")
        file = paths.ROOT / path
        return file.is_file() and re.search(rf"\b(?:theorem|lemma|def)\s+{re.escape(declaration)}\b",
                                             file.read_text(encoding="utf-8")) is not None
    path, _, function = reference.partition("::")
    file = paths.ROOT / path
    if not file.is_file():
        return False
    return function in {n.name for n in ast.parse(file.read_text(encoding="utf-8")).body if isinstance(n, ast.FunctionDef)}


def test_paper_claim_keys_resolve(paper, claims):
    unknown = sorted(set(latex.claim_keys(paper)) - set(claims))
    assert not unknown, f"claim keys in the paper with no row in data/claims.yaml: {unknown}"


def test_derived_claims_carry_their_proofs(paper, claims):
    in_appendix = latex.appendix_labels(paper)
    missing = [f"{key}: {claim.get('proof') or 'no proof named'}" for key, claim in claims.items()
               if claim.get("standing") == "derived" and claim.get("proof") not in in_appendix]
    assert not missing, f"a derived claim names the appendix that proves it in full: {missing}"


def test_claim_evidence_exists(claims):
    problems = []
    for key, claim in claims.items():
        if claim.get("status") not in STATUSES:
            problems.append(f"{key}: status {claim.get('status')!r}")
        if claim.get("standing") and claim["standing"] not in STANDINGS:
            problems.append(f"{key}: standing {claim['standing']!r}")
        if claim.get("status") == "checked" and not evidence_list(claim):
            problems.append(f"{key}: checked, with no evidence")
        problems += [f"{key}: no such evidence {e}" for e in evidence_list(claim) if not _exists(e)]
    assert not problems, problems
