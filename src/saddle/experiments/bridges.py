"""Data set beside the paper (Sec. IX). A measured entanglement-Hamiltonian profile is set beside the
readings, and the coldest equilibria a tightening bound allows are traced toward the saddle weights."""
from __future__ import annotations

import numpy as np

from repokit import datasets
from saddle import entropy, model, programs

CONFIGS = {
    "golden": {"ring_cells": 3, "ring_depth": 0.2, "gaps": [0.3, 0.03], "near": 2.0, "chain_cells": 3},
    "quick": {"ring_cells": 4, "ring_depth": 0.2, "gaps": [0.3, 0.03, 0.003], "near": 2.0, "chain_cells": 4},
    "paper": {"ring_cells": 6, "ring_depth": 0.2, "gaps": [0.3, 0.1, 0.03, 0.01, 0.003, 0.001], "near": 2.0,
              "chain_cells": 6},
}
TOLERANCE = {"rtol": 1e-4, "atol": 1e-7}
NOISE = ()
BOUND_BY = "double-precision dense linear algebra"


def _profile(config: dict) -> dict:
    rows = sorted(datasets.load("entanglement_profile"), key=lambda r: r["distance"])
    near = [r for r in rows if r["distance"] <= config["near"]]
    slope, intercept = np.polyfit([r["distance"] for r in near], [r["weight"] for r in near], 1)
    largest = max(r["weight"] for r in rows)
    return {
        "distances": [r["distance"] for r in rows], "weights": [r["weight"] for r in rows],
        "slope": float(slope), "intercept": float(intercept), "cut_fraction": float(intercept / largest),
        "nearest_over_largest": float(rows[0]["weight"] / largest),
        "increasing": bool(all(a["weight"] < b["weight"] for a, b in zip(rows, rows[1:]))),
    }


def _coldest(config: dict) -> dict:
    cells = config["ring_cells"]
    dens = model.ring(cells, 1.0, 1.0, [config["ring_depth"]] + [0.0] * (cells - 1))
    found = programs.saddle(dens)
    steps, multipliers = [], None
    for gap in config["gaps"]:
        multipliers = entropy.coldest_multipliers(dens, found["level"] + gap, start=multipliers)
        steps.append({"gap": gap, "total": float(multipliers.sum()),
                      "distance": float(np.abs(multipliers / multipliers.sum() - found["weights"]).max())})
    return {"cells": cells, "depth": config["ring_depth"], "weights": [float(w) for w in found["weights"]], "steps": steps}


def _coldest_chain(config: dict) -> dict:
    """The coldest equilibria of the open chain, whose optimal weightings form a segment: the proportions at each
    level, their two ends, and their largest distance from the even split, the segment's one symmetric point."""
    cells = config["chain_cells"]
    dens = model.build({"family": "ising", "lattice": "chain", "size": cells, "g": 1.0, "J": 1.0})
    level = programs.saddle(dens)["level"]
    even = np.zeros(cells)
    even[[0, -1]] = 0.5
    steps, multipliers = [], None
    for gap in config["gaps"]:
        multipliers = entropy.coldest_multipliers(dens, level + gap, start=multipliers)
        shape = multipliers / multipliers.sum()
        steps.append({"gap": gap, "total": float(multipliers.sum()), "shape": shape.tolist(),
                      "ends": [float(shape[0]), float(shape[-1])], "distance": float(np.abs(shape - even).max())})
    return {"cells": cells, "level": level, "steps": steps}


def run(config: dict) -> dict:
    return {"profile": _profile(config), "coldest": _coldest(config), "coldest_chain": _coldest_chain(config)}
