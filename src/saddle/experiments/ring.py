"""A ring of qubits, uniform and then with one cell made deeper.

One code path for both configurations: the saddle of the uniform ring set against
its exact ground energy per cell; a sweep of the deep cell's depth; a bisection for
the depth past which that cell carries no weight; and, on the uniform ring, a shallow one and a deep one, the
rate at which the common value moves when each cell's density is shifted, set against that cell's weight, with the
least and greatest weight each cell takes over the weightings near the common value, whose spread bounds how far the
optimal weighting can be from unique."""
from __future__ import annotations

import numpy as np

from saddle import model, programs

CONFIGS = {
    "paper": {
        "cells": 6, "field": 1.0, "coupling": 1.0,
        "depths": [0.0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.5, 0.75, 1.0], "shallow": 0.2, "deep": 1.0, "shift": 1e-3,
    },
    "quick": {
        "cells": 4, "field": 1.0, "coupling": 1.0,
        "depths": [0.0, 0.25, 0.5, 1.0], "shallow": 0.2, "deep": 1.0, "shift": 1e-3,
    },
    "golden": {
        "cells": 3, "field": 1.0, "coupling": 1.0,
        "depths": [0.0, 0.5, 1.0], "shallow": 0.2, "deep": 1.0, "shift": 1e-3,
    },
}
TOLERANCE = {"rtol": 1e-4, "atol": 1e-7}
# Solver residuals vary with the platform's linear algebra; the tests check them
# against stated bounds, so the golden output does not pin them.
NOISE = ("coarse", "gap", "constraint", "above", "slackness", "normalization")
NO_WEIGHT = 1e-6  # a weight below this counts as none
NEAR = 1e-7  # the weightings within this of the common value, over which each cell's weight is bounded
BOUND_BY = "double-precision dense linear algebra"


def _deep(config: dict, depth: float) -> list[np.ndarray]:
    cells = config["cells"]
    return model.ring(cells, config["field"], config["coupling"], [depth] + [0.0] * (cells - 1))


def _summary(result: dict) -> dict:
    keys = ("level", "gap", "constraint", "above", "slackness", "normalization")
    return {**{k: result[k] for k in keys},
            "weights": result["weights"].tolist(), "densities": result["densities"].tolist()}


def _rate(dens: list[np.ndarray], cell: int, step: float) -> float:
    """How fast the common value moves when one cell's density is shifted by a constant."""
    def shifted(amount: float) -> float:
        moved = list(dens)
        moved[cell] = moved[cell] + amount * np.eye(moved[cell].shape[0])
        return programs.saddle(moved)["level"]
    return (shifted(step) - shifted(-step)) / (2 * step)


def _rates(dens: list[np.ndarray], step: float) -> list[dict]:
    """Each cell's rate and weight, and the least and greatest weight it takes over the weightings within NEAR of
    the common value (Proposition rates): their spread bounds how far the optimal weighting is from unique."""
    weights = programs.saddle(dens)["weights"]
    out = []
    for cell in range(len(dens)):
        least, greatest = programs.rate_bounds(dens, cell, NEAR)
        out.append({"rate": _rate(dens, cell, step), "weight": float(weights[cell]), "least": least, "greatest": greatest})
    return out


def run(config: dict) -> dict:
    cells = config["cells"]
    uniform_dens = model.ring(cells, config["field"], config["coupling"])
    uniform = programs.saddle(uniform_dens)
    ground = float(np.linalg.eigvalsh(sum(uniform_dens))[0]) / cells
    sweep = [{"depth": depth, **_summary(programs.saddle(_deep(config, depth)))} for depth in config["depths"]]
    spin = model.on_site(model.PAULI_Z, 0, config["cells"])
    dropout = {"depth": programs.drop_depth(_deep(config, 0.0), 0, spin, tol=1e-9),
               "coarse": programs.drop_depth(_deep(config, 0.0), 0, spin, tol=1e-8)}
    deep_dens = _deep(config, config["deep"])
    deep = programs.saddle(deep_dens)
    return {
        "cells": cells,
        "ground_per_cell": ground,
        "uniform": _summary(uniform),
        "sweep": sweep,
        "dropout": dropout,
        "deep": {"depth": config["deep"], **_summary(deep)},
        "rates": {"uniform": _rates(uniform_dens, config["shift"]),
                  "shallow": _rates(_deep(config, config["shallow"]), config["shift"]),
                  "deep": _rates(deep_dens, config["shift"])},
    }
