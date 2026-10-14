"""Exits: instances where dropping a hypothesis of a result breaks it, or attempts where it did
not. The rotation of a uniform ring, dropped on an open chain; the nondegeneracy of the coldest
equilibrium, dropped with densities that are multiples of the identity; and the conservation of
the magnetization on an even XXZ ring, broken by a transverse field."""
from __future__ import annotations

import numpy as np

from saddle import entropy, model, programs

CONFIGS = {grade: {"chain_cells": 6, "identity_levels": [0.1, 0.2, 0.3], "identity_dim": 4,
                   "identity_gaps": [0.3, 0.03], "fields": [0.1, 0.3], "xxz_cells": 4}
           for grade in ("golden", "quick", "paper")}
TOLERANCE = {"rtol": 1e-4, "atol": 1e-7}
BOUND_BY = "double-precision dense linear algebra"


def run(config: dict) -> dict:
    cells = config["chain_cells"]
    chain = model.build({"family": "ising", "lattice": "chain", "size": cells, "g": 1.0, "J": 1.0})
    rotation = {"floor": programs.saddle(chain)["level"], "ground_per_cell": programs.uniform_value(chain)}
    flat = [c * np.eye(config["identity_dim"]) for c in config["identity_levels"]]
    floor = programs.saddle(flat)["level"]
    largest = max(float(entropy.coldest_multipliers(flat, floor + gap).max()) for gap in config["identity_gaps"])
    n = config["xxz_cells"]
    xxz = model.build({"family": "xxz", "lattice": "ring", "size": n, "J": 1.0, "delta": 0.5, "h": 0.2})
    spin = model.on_site(model.PAULI_Z, 0, n)
    broken = [{"field": f, "drop_depth": programs.drop_depth(
        [d - 0.5 * f * model.on_site(model.PAULI_X, x, n) for x, d in enumerate(xxz)], 0, spin)} for f in config["fields"]]
    return {"rotation": rotation, "degenerate": {"floor": floor, "largest_multiplier": largest}, "conservation": broken}
