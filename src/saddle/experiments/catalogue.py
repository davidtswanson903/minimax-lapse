"""The catalogue: every general claim of the paper, checked on many kinds of model.

Each model is a family of pieces (transverse-field Ising, XXZ, random) on a lattice
(ring, open chain, two-leg ladder, random graph) at one size and one choice of
parameters. On each model the experiment computes the saddle and checks: that the
two programs agree and the four conditions hold; that every place clearly below the common
value carries no weight; that the one-sided rates of the common value bracket the
first place's weight; on uniform rings, that the common value is the ground energy
per cell; and whether the first place, deepened, loses its weight.

Grades differ in breadth, not in precision: more families' parameters, more sizes,
more seeds. Models run in parallel when REPOKIT_WORKERS is set; the result does
not depend on how many workers run it."""
from __future__ import annotations

import os
import time

try:
    import resource  # Unix only: the heavy grades run there, and peak memory is recorded with them
except ImportError:
    resource = None
from concurrent.futures import ProcessPoolExecutor

import numpy as np

from saddle import model, programs

ENTRIES_QUICK = [
    {"family": "ising", "lattice": "ring", "sizes": [3, 4], "params": [{"g": 1.0, "J": 1.0}]},
    {"family": "random", "lattice": "random", "sizes": [4], "params": [{"seed": 0}]},
]


def _entries(sizes: list[int], g: list[float], delta: list[float], seeds: list[int]) -> list[dict]:
    """The kinds of model, at the given sizes. Breadth comes from the kinds: families of
    pieces, lattices, signs of the coupling, a longitudinal field, and which cell is
    deepened. Size is not breadth; the scaling experiment takes it up separately. Each lattice starts at the least
    size at which it is not another lattice relabelled: the complete graph on three cells is the three-cell ring, the
    three-cell star the three-cell chain, and the two-leg ladder of four cells the four-cell ring, so complete graphs
    and stars start at four cells and the ladder at six."""
    ising = lambda **extra: [{"g": v, "J": 1.0, **extra} for v in g]
    four = [n for n in sizes if n >= 4]
    return [
        {"family": "ising", "lattice": "ring", "sizes": sizes, "params": ising()},
        {"family": "ising", "lattice": "ring", "sizes": sizes,
         "params": [{"g": v, "J": -1.0} for v in g] + [{"g": v, "J": 1.0, "h": 0.5} for v in g]},
        {"family": "ising", "lattice": "chain", "sizes": sizes, "params": ising() + ising(deep_at="middle")},
        {"family": "ising", "lattice": "ladder", "sizes": [n for n in sizes if n >= 6 and n % 2 == 0] or [6], "params": ising()},
        {"family": "ising", "lattice": "complete", "sizes": four, "params": ising() + [{"g": v, "J": -1.0} for v in g]},
        {"family": "ising", "lattice": "star", "sizes": four, "params": ising() + ising(deep_at="last")},
        {"family": "xxz", "lattice": "ring", "sizes": sizes, "params": [{"J": 1.0, "delta": d, "h": 0.2} for d in delta]},
        {"family": "xxz", "lattice": "chain", "sizes": sizes,
         "params": [{"J": 1.0, "delta": d, "h": 0.2, "deep_at": at} for d in delta for at in (None, "middle")]},
        {"family": "xxz", "lattice": "star", "sizes": four,
         "params": [{"J": 1.0, "delta": d, "h": 0.2, "deep_at": at} for d in delta for at in (None, "last") if at or d]},
        {"family": "random", "lattice": "random", "sizes": sizes, "params": [{"seed": s} for s in seeds]},
    ]


CONFIGS = {
    "golden": {"models": [
        {"family": "ising", "lattice": "ring", "sizes": [3], "params": [{"g": 1.0, "J": 1.0}]},
        {"family": "ising", "lattice": "complete", "sizes": [4], "params": [{"g": 1.0, "J": -1.0, "h": 0.3}]},
        {"family": "ising", "lattice": "star", "sizes": [4], "params": [{"g": 1.0, "J": 1.0, "deep_at": "last"}]},
        {"family": "ising", "lattice": "ladder", "sizes": [6], "params": [{"g": 1.0, "J": 1.0}]},
        {"family": "ising", "lattice": "chain", "sizes": [3], "params": [{"g": 1.0, "J": 1.0, "deep_at": "middle"}]},
        {"family": "xxz", "lattice": "chain", "sizes": [3], "params": [{"J": 1.0, "delta": 0.5, "h": 0.2}]},
        {"family": "random", "lattice": "random", "sizes": [3], "params": [{"seed": 0}]},
    ], "deep": 2.0, "shift": 1e-3},
    "quick": {"models": ENTRIES_QUICK, "deep": 2.0, "shift": 1e-3},
    "survey": {"models": _entries([3, 4, 5], [0.5, 1.0, 2.0], [-0.5, 0.5, 1.5], [0, 1, 2]), "deep": 2.0, "shift": 1e-3},
    "certification": {"models": _entries([3, 4, 5], [0.25, 0.5, 1.0, 2.0, 4.0], [-1.5, -0.5, 0.0, 0.5, 1.0, 1.5],
                                         list(range(12))), "deep": 2.0, "shift": 1e-3},
}
TOLERANCE = {"rtol": 1e-5, "atol": 1e-7}
NOISE = ("gap", "constraint", "above", "slackness", "normalization", "forward", "backward",
         "weights", "deep_weight", "deep_margin", "below", "below_weight_max", "seconds", "rss_mb", "worst")
BOUND_BY = "double-precision dense linear algebra"  # what the heavy grades are limited by
# A place counts as clearly below the common value at a margin of 1e-3, well clear of
# flat optima; an interior-point solver leaves weights of order gap/margin there, so
# no weight means below 1e-5.
NO_WEIGHT = 1e-5
BELOW = 1e-3
RATE_SLACK = 1e-4


def specs(config: dict) -> list[dict]:
    return [{"family": e["family"], "lattice": e["lattice"], "size": size,
             **{k: v for k, v in params.items() if v is not None}}
            for e in config["models"] for size in e["sizes"] for params in e["params"]]


def _shifted_level(dens, step):
    moved = list(dens)
    moved[0] = moved[0] + step * np.eye(moved[0].shape[0])
    return programs.saddle(moved)["level"]


def check_model(spec: dict, deep: float, shift: float) -> dict:
    start = time.time()
    dens = model.build(spec)
    result = programs.saddle(dens)
    level, weights = result["level"], result["weights"]
    below = result["densities"] < level - BELOW
    forward = (_shifted_level(dens, shift) - level) / shift
    backward = (level - _shifted_level(dens, -shift)) / shift
    uniform = None
    if model.translation_invariant(spec):
        uniform = bool(abs(level - np.linalg.eigvalsh(sum(dens))[0] / spec["size"]) < 1e-6)
    at = model.deep_index(spec)
    deep_result = programs.saddle(model.build({**spec, "depth": deep}))
    margin = deep_result["level"] - deep_result["densities"][at]
    return {
        "spec": spec,
        "level": level,
        "programs": bool(result["gap"] < 1e-7 and result["constraint"] < 1e-4
                         and max(result["above"], result["slackness"], result["normalization"]) < 1e-7),
        "forbids": bool(np.all(weights[below] < NO_WEIGHT)),
        "rates": bool(forward >= weights[0] - RATE_SLACK and backward <= weights[0] + RATE_SLACK),
        "uniform": uniform,
        "dropped": bool(deep_result["weights"][at] < NO_WEIGHT and margin > BELOW),
        "weighted": bool(weights[at] >= NO_WEIGHT),
        "below": int(below.sum()),
        "below_weight_max": float(np.max(weights[below])) if below.any() else 0.0,
        "weights": weights.tolist(),
        "forward": forward, "backward": backward,
        "deep_weight": float(deep_result["weights"][at]), "deep_margin": float(margin),
        **{k: result[k] for k in ("gap", "constraint", "above", "slackness", "normalization")},
        "seconds": time.time() - start,
        "rss_mb": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024 if resource else None,
    }


def undeepened(spec: dict) -> tuple:
    """A model before deepening: the variants that differ only in which cell is deepened are one model."""
    return tuple(sorted((k, v) for k, v in spec.items() if k != "deep_at"))


def distinct_below(rows: list) -> int:
    """The places below the common value, counted once for each distinct model before deepening."""
    return sum({undeepened(r["spec"]): r["below"] for r in rows}.values())


def _check(args):
    return check_model(*args)


def run(config: dict) -> dict:
    jobs = [(spec, config["deep"], config["shift"]) for spec in specs(config)]
    workers = int(os.environ.get("REPOKIT_WORKERS", "1"))
    if workers > 1:
        with ProcessPoolExecutor(workers) as pool:
            rows = list(pool.map(_check, jobs))
    else:
        rows = [_check(job) for job in jobs]
    uniform = [r["uniform"] for r in rows if r["uniform"] is not None]
    return {
        "models": rows,
        "summary": {
            "models": len(rows),
            "families": len({r["spec"]["family"] for r in rows}),
            "lattices": len({r["spec"]["lattice"] for r in rows}),
            "largest": max(r["spec"]["size"] for r in rows),
            "programs": sum(r["programs"] for r in rows),
            "forbids": sum(r["forbids"] for r in rows),
            "below_places": distinct_below(rows),
            "rates": sum(r["rates"] for r in rows),
            "uniform_checked": len(uniform),
            "uniform": sum(uniform),
            "dropped": sum(r["dropped"] for r in rows),
            "worst": {k: max(r[k] for r in rows)
                      for k in ("gap", "constraint", "above", "slackness", "normalization", "below_weight_max")},
        },
    }
