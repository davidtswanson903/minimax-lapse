"""Experiments, grades, golden checks and provenance.

An experiment is a module in the work's `experiments` package with CONFIGS, one
configuration per grade, every one with the golden configuration's keys, and
run(config) -> dict. Grades run the same code path at widening breadth:

    golden          pins the logic; runs on every check
    quick           what a reader reruns after downloading, in seconds
    paper           the numbers the paper prints
    survey          a sample of the certification catalogue, on one machine
    certification   the full catalogue, routed by the private plan

The golden output pins the logic: while it holds, stored results of every grade
stay valid without being rerun. Each grade's results, in results/<grade>/, record
the golden output they were made under, the machine, and the measured run time.

    python -m repokit.golden check               golden outputs, and every stored grade's provenance
    python -m repokit.golden pin [NAME...]       store new golden outputs, deliberately
    python -m repokit.golden run GRADE [NAME...] run a grade; refuses while a golden check fails
    python -m repokit.golden reproduce [GRADE]   rerun a grade (default quick) and compare; writes nothing
"""
from __future__ import annotations

import argparse
import hashlib
import importlib
import json
import os
import pathlib
import pkgutil
import platform
import shutil
import sys
import time

from . import paths

GRADES = ("golden", "quick", "paper", "survey", "certification")
DEFAULT_TOLERANCE = {"rtol": 1e-12, "atol": 1e-14}


def experiments() -> dict:
    package = importlib.import_module(f"{paths.package()}.experiments")
    found = {}
    for info in pkgutil.iter_modules(package.__path__):
        if info.name.startswith("_"):
            continue
        module = importlib.import_module(f"{package.__name__}.{info.name}")
        if hasattr(module, "run") and hasattr(module, "CONFIGS"):
            found[info.name] = module
    return found


def canonical(obj) -> str:
    return json.dumps(obj, sort_keys=True, separators=(",", ":"), ensure_ascii=False)


def digest(obj) -> str:
    return hashlib.sha256(canonical(obj).encode("utf-8")).hexdigest()


def _normal(obj):
    return json.loads(json.dumps(obj))


def _load(path):
    return json.loads(path.read_text(encoding="utf-8"))


def _write(path, obj) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(obj, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8", newline="\n")


def grade_path(name: str, grade: str):
    return paths.RESULTS / grade / f"{name}.json"


def golden_path(name: str):
    return grade_path(name, "golden")


def tolerance(module) -> dict:
    return {**DEFAULT_TOLERANCE, **getattr(module, "TOLERANCE", {})}


def key_tolerances(module) -> dict[str, dict]:
    """Looser tolerances for named keys, wherever they occur, each stated with its reason in the experiment: a fit
    that amplifies the noise in what it fits, or an iterate an optimizer stops short of exact."""
    return {key: {**DEFAULT_TOLERANCE, **tol} for key, tol in getattr(module, "TOLERANCE_BY_KEY", {}).items()}


def noise(module) -> frozenset:
    return frozenset(getattr(module, "NOISE", ()))


def differences(new, old, rtol: float, atol: float, where: str = "", noise=frozenset(),
                by_key: dict | None = None) -> list[str]:
    """Differences beyond tolerance. Keys named as noise (solver residuals, timings,
    weights that need not be unique) are skipped wherever they occur, and a key with its own
    tolerance is compared under it, with everything beneath it."""
    if isinstance(old, dict) and isinstance(new, dict):
        out: list[str] = []
        for key in sorted(set(old) | set(new)):
            if key in noise:
                continue
            if key not in old:
                out.append(f"{where}.{key}: not in the stored output")
            elif key not in new:
                out.append(f"{where}.{key}: no longer produced")
            else:
                own = (by_key or {}).get(key)
                r, a = (own["rtol"], own["atol"]) if own else (rtol, atol)
                out += differences(new[key], old[key], r, a, f"{where}.{key}", noise, by_key)
        return out
    if isinstance(old, list) and isinstance(new, list):
        if len(old) != len(new):
            return [f"{where}: {len(new)} items, stored {len(old)}"]
        out = []
        for i, (a, b) in enumerate(zip(new, old)):
            out += differences(a, b, rtol, atol, f"{where}[{i}]", noise, by_key)
        return out
    numbers = (int, float)
    if isinstance(old, numbers) and isinstance(new, numbers) and not isinstance(old, bool) and not isinstance(new, bool):
        return [] if abs(new - old) <= atol + rtol * abs(old) else [f"{where}: {new!r}, stored {old!r}"]
    return [] if new == old else [f"{where}: {new!r}, stored {old!r}"]


def grades_share_keys(module) -> bool:
    keys = set(module.CONFIGS["golden"])
    return all(set(config) == keys for config in module.CONFIGS.values())


def run_config(module, grade: str) -> dict:
    return _normal(module.run(dict(module.CONFIGS[grade])))


def check_golden(name: str, module) -> list[str]:
    path = golden_path(name)
    if not path.exists():
        return [f"{name}: no stored golden output; pin it with `make golden EXP={name}`"]
    tol = tolerance(module)
    return differences(run_config(module, "golden"), _load(path), tol["rtol"], tol["atol"], name, noise(module),
                       key_tolerances(module))


def stored_grades(name: str) -> list[str]:
    return [g for g in GRADES if g != "golden" and grade_path(name, g).exists()]


def check_provenance(name: str, grade: str) -> list[str]:
    path, gold = grade_path(name, grade), golden_path(name)
    if not path.exists():
        return [f"{name}: no stored {grade} results; run `make run GRADE={grade} EXP={name}`"]
    if not gold.exists():
        return [f"{name}: no stored golden output"]
    if _load(path).get("_provenance", {}).get("golden_sha256") != digest(_load(gold)):
        return [f"{name}: the stored {grade} results were made under a different golden output; "
                f"rerun with `make run GRADE={grade} EXP={name}`"]
    return []


def machine() -> dict:
    return {"platform": platform.platform(), "processor": platform.processor() or platform.machine(),
            "cpus": os.cpu_count(), "workers": int(os.environ.get("REPOKIT_WORKERS", "1")),
            "python": platform.python_version()}


def _select(names: list[str], grade: str | None = None) -> dict:
    found = experiments()
    unknown = [n for n in names if n not in found]
    if unknown:
        raise SystemExit(f"unknown experiments: {unknown}; known: {sorted(found)}")
    chosen = {n: found[n] for n in (names or sorted(found))}
    if grade is not None:
        chosen = {n: m for n, m in chosen.items() if grade in m.CONFIGS}
    return chosen


def pin(names: list[str]) -> None:
    for name, module in _select(names).items():
        _write(golden_path(name), run_config(module, "golden"))
        print(f"pinned the golden output of {name}")


def run_grade(grade: str, names: list[str]) -> None:
    if grade not in GRADES or grade == "golden":
        raise SystemExit(f"grade {grade!r} is not one of {GRADES[1:]}")
    for name, module in _select(names, grade).items():
        problems = check_golden(name, module)
        if problems:
            raise SystemExit(
                f"{name}: its golden check fails, so its logic changed. Pin the new golden output "
                f"deliberately (`make golden EXP={name}`), then run it.\n  " + "\n  ".join(problems)
            )
        start = time.time()
        out = run_config(module, grade)
        out["_provenance"] = {"grade": grade, "config": module.CONFIGS[grade],
                              "golden_sha256": digest(_load(golden_path(name))),
                              "seconds": round(time.time() - start, 3), "machine": machine()}
        _write(grade_path(name, grade), out)
        handed = hand_over(name, grade, module)
        print(f"ran {name} at the {grade} grade in {out['_provenance']['seconds']:.1f} s" + (f"; left in {handed}" if handed else ""))


def sources_digest(module) -> dict:
    """The hashes of the code a result depends on: the experiment's module and every module of the shared base
    and of the work's package, by their paths in the public tree."""
    files = sorted({pathlib.Path(module.__file__), *(paths.ROOT / "src").rglob("*.py")})
    return {paths.rel(f): hashlib.sha256(f.read_bytes().replace(b"\r\n", b"\n")).hexdigest() for f in files
            if "__pycache__" not in f.parts}


def hand_over(name: str, grade: str, module) -> pathlib.Path | None:
    """Where a run is told to leave its result for collection, copy it there with a receipt.

    REPOKIT_RESULTS_FROM names the environment variable that holds the directory, and REPOKIT_COMMIT_FROM the one
    that holds the commit the run was made at; passing names rather than values keeps a server's own variables
    out of the committed tree. Unset, the run is local and nothing is handed over."""
    holder = os.environ.get("REPOKIT_RESULTS_FROM")
    if not holder or not os.environ.get(holder):
        return None
    where = pathlib.Path(os.environ[holder])
    where.mkdir(parents=True, exist_ok=True)
    result = where / f"{name}.{grade}.json"
    shutil.copyfile(grade_path(name, grade), result)
    receipt = {"experiment": name, "grade": grade, "taken": os.environ.get("REPOKIT_TAKEN", "local"),
               "commit": os.environ.get(os.environ.get("REPOKIT_COMMIT_FROM", ""), None),
               "result_sha256": hashlib.sha256(result.read_bytes()).hexdigest(), "sources": sources_digest(module)}
    _write(where / f"{name}.{grade}.receipt.json", receipt)
    return where


def reproduction_differences(name: str, module, stored: dict, grade: str = "quick") -> list[str]:
    """Rerun one grade and compare it with its stored results; writes nothing."""
    kept = {k: v for k, v in stored.items() if not k.startswith("_")}
    tol = tolerance(module)
    return differences(run_config(module, grade), kept, tol["rtol"], tol["atol"], name, noise(module),
                       key_tolerances(module))


def reproduce(grade: str, names: list[str]) -> int:
    failures = 0
    for name, module in _select(names, grade).items():
        path = grade_path(name, grade)
        if not path.exists():
            print(f"{name}: no stored {grade} results to compare with")
            failures += 1
            continue
        found = reproduction_differences(name, module, _load(path), grade)
        if found:
            failures += 1
            print(f"{name}: does not reproduce its {grade} results\n  " + "\n  ".join(found[:20]))
        else:
            print(f"{name}: reproduces its stored {grade} results, {path.relative_to(paths.ROOT).as_posix()}")
    return 1 if failures else 0


def check(names: list[str]) -> int:
    problems: list[str] = []
    for name, module in _select(names).items():
        problems += check_golden(name, module)
        for grade in stored_grades(name):
            problems += check_provenance(name, grade)
    for line in problems:
        print(line)
    if not problems:
        print("golden outputs hold, and every stored grade follows the current logic")
    return 1 if problems else 0


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(prog="python -m repokit.golden")
    sub = parser.add_subparsers(dest="action", required=True)
    sub.add_parser("check").add_argument("names", nargs="*")
    sub.add_parser("pin").add_argument("names", nargs="*")
    runner = sub.add_parser("run")
    runner.add_argument("grade")
    runner.add_argument("names", nargs="*")
    again = sub.add_parser("reproduce")
    again.add_argument("grade", nargs="?", default="quick")
    again.add_argument("names", nargs="*")
    args = parser.parse_args(argv)
    if args.action == "pin":
        pin(args.names)
        return 0
    if args.action == "run":
        run_grade(args.grade, args.names)
        return 0
    if args.action == "reproduce":
        return reproduce(args.grade, args.names)
    return check(args.names)


if __name__ == "__main__":
    sys.exit(main())
