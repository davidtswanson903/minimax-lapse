"""Golden outputs pin each experiment's logic, and the stored results follow it."""
from types import SimpleNamespace

from repokit import golden, sources


def test_golden_outputs_hold():
    problems = [p for name, module in golden.experiments().items() for p in golden.check_golden(name, module)]
    assert not problems, problems


def test_stored_results_follow_current_logic():
    problems = [p for name in golden.experiments() for grade in golden.stored_grades(name)
                for p in golden.check_provenance(name, grade)]
    assert not problems, problems


def test_every_grade_shares_the_golden_keys():
    wrong = [name for name, module in golden.experiments().items()
             if not golden.grades_share_keys(module) or set(module.CONFIGS) - set(golden.GRADES)]
    assert not wrong, f"every grade runs the golden configuration's code path, under a known grade name: {wrong}"


def test_printed_grades_are_stored(paper):
    chosen = sources.numbers_from()
    missing = [f"{name}: {grade}" for name, grade in chosen.items() if not golden.grade_path(name, grade).exists()]
    unnamed = [grade for grade in set(chosen.values()) if grade not in paper]
    assert not missing, f"the surfaces print grades that are not stored: {missing}"
    assert not unnamed, f"the paper does not name the grade its numbers come from: {unnamed}"


def test_reproduce_compares_and_writes_nothing(tmp_path, monkeypatch):
    monkeypatch.setattr(golden.paths, "RESULTS", tmp_path)
    experiment = SimpleNamespace(CONFIGS={"full": {"x": 2.0}, "golden": {"x": 1.0}},
                                 run=lambda config: {"y": config["x"] ** 2, "residual": 1e-9}, NOISE=("residual",))
    experiment.CONFIGS["quick"] = experiment.CONFIGS["full"]
    assert golden.reproduction_differences("toy", experiment, {"y": 4.0, "residual": 5e-9, "_provenance": {}}) == []
    assert golden.reproduction_differences("toy", experiment, {"y": 4.1, "residual": 0.0})
    assert list(tmp_path.iterdir()) == []


def test_a_run_told_where_hands_over_its_result_with_a_receipt(tmp_path, monkeypatch):
    """With REPOKIT_RESULTS_FROM naming a variable that holds a directory, a run leaves its result there with a
    receipt: the result's hash, the commit named by REPOKIT_COMMIT_FROM, and the hashes of the code it ran."""
    import hashlib
    import json
    name, module = next((n, m) for n, m in sorted(golden.experiments().items()) if golden.grade_path(n, "quick").exists())
    monkeypatch.setenv("WHERE_RESULTS_GO", str(tmp_path))
    monkeypatch.setenv("WHICH_COMMIT", "0123abc")
    monkeypatch.setenv("REPOKIT_RESULTS_FROM", "WHERE_RESULTS_GO")
    monkeypatch.setenv("REPOKIT_COMMIT_FROM", "WHICH_COMMIT")
    monkeypatch.setenv("REPOKIT_TAKEN", "remote")
    assert golden.hand_over(name, "quick", module) == tmp_path
    result = (tmp_path / f"{name}.quick.json").read_bytes()
    receipt = json.loads((tmp_path / f"{name}.quick.receipt.json").read_text(encoding="utf-8"))
    assert receipt["result_sha256"] == hashlib.sha256(result).hexdigest()
    assert receipt["commit"] == "0123abc" and receipt["taken"] == "remote"
    assert any(path.endswith(f"experiments/{name}.py") for path in receipt["sources"])
    monkeypatch.delenv("REPOKIT_RESULTS_FROM")
    assert golden.hand_over(name, "quick", module) is None


def test_a_named_key_takes_its_own_tolerance():
    """A key with its own tolerance is compared under it, with everything beneath it; the rest keep the module's."""
    old = {"fits": [{"limit": 1.0}], "depth": 1.0}
    new = {"fits": [{"limit": 1.0005}], "depth": 1.0005}
    found = golden.differences(new, old, 1e-6, 1e-9, "x", frozenset(), {"fits": {"rtol": 1e-3, "atol": 1e-9}})
    assert found == ["x.depth: 1.0005, stored 1.0"], found
