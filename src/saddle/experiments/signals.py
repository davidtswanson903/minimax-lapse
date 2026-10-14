"""Signals under the weighted evolution: how long a disturbance takes to cross the open chain.

The state the Hamiltonian reweighted by place holds fixed is kicked at the first cell inside the chain's ends,
and evolves under that Hamiltonian; each place's density is watched until it moves from its value just after
the kick. Along weightings running from uniform to the saddle's, every interior cell carries the same weight, and
the arrival at the last interior cell is timed against it. At the saddle the interior carries no weight, the
reweighted Hamiltonian has no term across the interior's links, and the disturbance is followed over a long
window to measure how far it gets."""
from __future__ import annotations

import numpy as np
from scipy.linalg import expm

from saddle import model, programs

PAULI_Y = np.array([[0.0, -1j], [1j, 0.0]])
KICK = np.pi / 8  # the kick: a rotation of the first interior cell about Y by this angle
KICKED = 1  # the kicked cell, the first inside the chain's ends
THRESHOLD = 1e-3  # a density has moved once it departs from its stationary value by this much
WINDOW = 8.0  # how long a density is watched, in units of one over the interior's weight
SADDLE_WINDOW = 1e4  # at the saddle, where the interior's weight is zero, the window is this long
NO_WEIGHT = 1e-6  # a weight below this is the solver's zero

CONFIGS = {
    "golden": {"cells": 5, "field": 1.0, "coupling": 1.0, "approach": [0.0, 0.9, 0.99, 1.0], "steps": 200},
    "quick": {"cells": 5, "field": 1.0, "coupling": 1.0, "approach": [0.0, 0.5, 0.9, 0.99, 1.0], "steps": 400},
    "paper": {"cells": 6, "field": 1.0, "coupling": 1.0,
              "approach": [0.0, 0.5, 0.8, 0.9, 0.95, 0.99, 0.999, 1.0], "steps": 800},
}
TOLERANCE = {"rtol": 1e-3, "atol": 1e-6}
NOISE = ("wall",)
BOUND_BY = "the saddle's semidefinite program, whose state has dimension two to the number of cells"


def chain(config: dict) -> list[np.ndarray]:
    """The densities of the open transverse-field chain."""
    return model.build({"family": "ising", "lattice": "chain", "size": config["cells"],
                        "g": config["field"], "J": config["coupling"]})


def interior_weight(weights: np.ndarray) -> float:
    """The mean weight of the cells between the two ends."""
    return float(np.mean(weights[1:-1]))


def watch(dens: list[np.ndarray], weights: np.ndarray, span: float, steps: int) -> dict:
    """Kick the first interior cell of the state the weighted Hamiltonian holds fixed, evolve, and return for
    each cell the first time its density moves by the threshold from its value just after the kick (None if
    never within the span), and the largest such move over the span."""
    cells = len(dens)
    energies, vectors = np.linalg.eigh(programs.weighted(dens, weights))
    held = vectors[:, 0]
    kicked = expm(-1j * KICK * model.on_site(PAULI_Y, KICKED, cells)) @ held
    start = vectors.conj().T @ kicked
    rest = np.array([np.real(np.vdot(kicked, d @ kicked)) for d in dens])
    first, largest = [None] * cells, np.zeros(cells)
    for t in np.linspace(0.0, span, steps):
        psi = vectors @ (np.exp(-1j * energies * t) * start)
        moved = np.abs(np.array([np.real(np.vdot(psi, d @ psi)) for d in dens]) - rest)
        largest = np.maximum(largest, moved)
        for j in range(cells):
            if first[j] is None and moved[j] > THRESHOLD:
                first[j] = float(t)
    return {"arrivals": first, "largest": [float(v) for v in largest]}


def run(config: dict) -> dict:
    dens = chain(config)
    cells = config["cells"]
    found = programs.saddle(dens)
    saddle = np.array(found["weights"])
    uniform = np.full(cells, 1.0 / cells)
    path = []
    for s in config["approach"]:
        weights = (1.0 - s) * uniform + s * saddle
        inner = interior_weight(weights)
        span = WINDOW / inner if inner > NO_WEIGHT else SADDLE_WINDOW
        seen = watch(dens, weights, span, config["steps"])
        path.append({"approach": s, "interior_weight": inner, "span": span, **seen,
                     "far": seen["arrivals"][cells - 2]})
    fitted = [p for p in path if p["far"] is not None and p["interior_weight"] > NO_WEIGHT]
    exponent = None
    if len(fitted) >= 2:
        exponent = float(np.polyfit(np.log([p["interior_weight"] for p in fitted]),
                                    np.log([p["far"] for p in fitted]), 1)[0])
    at_saddle = path[-1]
    return {
        "cells": cells,
        "saddle_weights": [float(w) for w in saddle],
        "interior_saddle_weight": interior_weight(saddle),
        "path": path,
        "exponent": exponent,
        "reached_at_saddle": max([j for j, a in enumerate(at_saddle["arrivals"]) if a is not None], default=0),
        "wall": max(at_saddle["largest"][KICKED + 1:]),
        "uniform_far": path[0]["far"],
    }
