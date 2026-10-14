"""Scaling: one kind of model across sizes, for the one finding that depends on size.

The general claims hold at every size, so size adds nothing to them. What size can
change is a measured, family-specific finding: whether a deepened cell drops out,
and at what depth. For Ising rings this experiment brackets that depth at each size;
for XXZ rings it records whether the deepened cell drops out at a fixed depth."""
from __future__ import annotations

import itertools

import numpy as np

from saddle import laws, model, programs
from saddle.experiments.catalogue import BELOW, NO_WEIGHT

CONFIGS = {
    "golden": {"ising_sizes": [3], "xxz_sizes": [3], "deep": 2.0,
               "xxz_depths": [0.0, 2.0], "xxz_fields": [0.2], "anisotropies": [1.5]},
    "quick": {"ising_sizes": [3, 4], "xxz_sizes": [3, 4], "deep": 2.0,
              "xxz_depths": [0.0, 1.0, 2.0], "xxz_fields": [0.1, 0.2], "anisotropies": [0.5, 1.5]},
    "survey": {"ising_sizes": [3, 4, 5, 6], "xxz_sizes": [3, 4, 5, 6], "deep": 2.0,
               "xxz_depths": [0.0, 0.1, 0.2, 0.3, 0.4, 0.6, 0.8, 1.0, 1.5, 2.0, 3.0, 4.0, 6.0, 8.0],
               "xxz_fields": [0.1, 0.2, 0.4], "anisotropies": [-0.5, 0.5, 1.0, 1.5]},
    "certification": {"ising_sizes": [3, 4, 5, 6, 7], "xxz_sizes": [3, 4, 5, 6, 7], "deep": 2.0,
                      "xxz_depths": [0.0, 0.1, 0.2, 0.3, 0.4, 0.6, 0.8, 1.0, 1.5, 2.0, 3.0, 4.0, 6.0, 8.0, 12.0, 16.0],
                      "xxz_fields": [0.05, 0.1, 0.2, 0.4, 0.8], "anisotropies": [-0.5, 0.5, 1.0, 1.5]},
}
TOLERANCE = {"rtol": 1e-4, "atol": 1e-7}
COUPLING, ANISOTROPY = 1.0, 0.5  # the XXZ rings' J and Delta
CURVE_FIELD = 0.2  # the field of the XXZ drop-out depths and curves
PAIR_BOUND = 1e-7  # how far above e' a density may sit and still count as held at it
PAIR_LEVEL = 1e-6  # how near e' a sector's lowest level must be to belong to the others' ground level
PAIR_FLAT = 1e-4  # a density that changes less than this across the pair counts as flat: on three cells the two
# others are equal by reflection, up to the solver's weights, which a floor flat to second order leaves this loose
PAIR_FLOOR = 1e-4  # how near e' the pair's own least largest density must be: on seven cells the solver's weights
# leave it up to 1e-5 above, two of the bounds' lines being nearly flat
NOISE = ("weight", "margin", "product", "other", "predicted", "coarse", "lift")
BOUND_BY = "double-precision dense linear algebra"


def _ising(cells: int, depth: float) -> list:
    return model.build({"family": "ising", "lattice": "ring", "size": cells, "g": 1.0, "J": 1.0, "depth": depth})


def _xxz(cells: int, depth: float, field: float = CURVE_FIELD, delta: float = ANISOTROPY) -> list:
    return model.build({"family": "xxz", "lattice": "ring", "size": cells, "J": COUPLING, "delta": delta, "h": field,
                        "depth": depth})


def _pair(dens: list) -> dict | None:
    """The other cells' optimal states on an XXZ ring whose deepened cell is the first, when the ground level of their
    optimal weighted Hamiltonian is one state of magnetization M = +1, one of M = -1, or a pair of the two. Every
    optimal state of theirs lies in that level, and every density conserves M, so a state's densities there depend
    only on its share p of M = +1; the bounds Tr(rho D_x) <= e' fix an interval of p, and the deepened spin over that
    interval decides the drop-out: its least depth is the least ratio of Eq. dropout at the interval's ends with
    positive spin, and it never comes when the spin is nowhere positive (Proposition dropout). None when the level is
    not of that kind."""
    others = dens[1:]
    upper = programs.greatest_least(others)
    level = upper["value"]
    weighted = sum(w * d for w, d in zip(upper["weights"], others))
    cells = int(round(np.log2(dens[0].shape[0])))
    magnetization = np.diag(sum(model.on_site(model.PAULI_Z, x, cells) for x in range(cells)))
    at_level, states = [], {}
    for sector in np.unique(magnetization):
        keep = np.isclose(magnetization, sector)
        values, vectors = np.linalg.eigh(weighted[np.ix_(keep, keep)])
        if values[0] - level < PAIR_LEVEL:
            if len(values) > 1 and values[1] - level < 1e-3:
                return None
            at_level.append(float(sector))
        states[float(sector)] = np.zeros(len(keep))
        states[float(sector)][keep] = vectors[:, 0]
    if not at_level or not set(at_level) <= {1.0, -1.0}:
        return None
    up, down = states[1.0], states[-1.0]
    held = lambda d, p: p * (up @ d @ up) + (1 - p) * (down @ d @ down)
    least, most = (0.0 if -1.0 in at_level else 1.0), (1.0 if 1.0 in at_level else 0.0)
    # The bounds' lines all pass through e' at the optimal shares, up to rounding, so the interval is taken around
    # the pair's own least largest density rather than at e' itself.
    lines = [(held(d, 0.0), held(d, 1.0) - held(d, 0.0)) for d in others] if len(at_level) == 2 else []
    lines = [(b, s) for b, s in lines if abs(s) > PAIR_FLAT]
    if lines:
        peak = lambda p: max(b + s * p for b, s in lines)
        crossings = [(b2 - b1) / (s1 - s2) for (b1, s1), (b2, s2) in itertools.combinations(lines, 2) if s1 != s2]
        floor = min(peak(p) for p in [0.0, 1.0] + [p for p in crossings if 0.0 <= p <= 1.0])
        if abs(floor - level) > PAIR_FLOOR:
            return None
        for b, s in lines:
            edge = (floor + PAIR_BOUND - b) / s
            least, most = (least, min(most, edge)) if s > 0 else (max(least, edge), most)
    spin = model.on_site(model.PAULI_Z, 0, cells)
    ends = [least, most]
    spins = [float(held(spin, p)) for p in ends]
    ratios = [max(0.0, float(held(dens[0], p) - level) / s) for p, s in zip(ends, spins) if s > 0]
    return {"share": ends, "spin": spins, "depth": min(ratios) if ratios else None, "others": float(level)}


def _dropped(dens: list) -> tuple[bool, float, float]:
    result = programs.saddle(dens)
    margin = result["level"] - result["densities"][0]
    return bool(result["weights"][0] < NO_WEIGHT and margin > BELOW), float(result["weights"][0]), float(margin)


def run(config: dict) -> dict:
    ising = []
    for cells in config["ising_sizes"]:
        spin = model.on_site(model.PAULI_Z, 0, cells)
        ising.append({"cells": cells, "depth": programs.drop_depth(_ising(cells, 0.0), 0, spin, tol=1e-9),
                      "coarse": programs.drop_depth(_ising(cells, 0.0), 0, spin, tol=1e-8)})
    xxz = []
    for cells in config["xxz_sizes"]:
        dropped, weight, margin = _dropped(_xxz(cells, config["deep"]))
        spin = model.on_site(model.PAULI_Z, 0, cells)
        exact = programs.drop_depth(_xxz(cells, 0.0), 0, spin, tol=1e-9)
        curve = []
        for depth in config["xxz_depths"]:
            down, w, m = _dropped(_xxz(cells, depth))
            curve.append({"depth": depth, "dropped": down, "weight": w, "margin": m})
        xxz.append({"cells": cells, "dropped": dropped, "weight": weight, "margin": margin,
                    "exact": exact, "curve": curve})
    # On even rings the weight falls as the inverse of the depth; the product of the two at
    # the deepest depth, against the weak field the other cells carry and the derived limit.
    deepest = max(config["xxz_depths"])
    fields = []
    for cells in [c for c in config["xxz_sizes"] if c % 2 == 0]:
        for field in config["xxz_fields"]:
            weights = programs.saddle(_xxz(cells, deepest, field))["weights"]
            fields.append({"cells": cells, "field": field, "product": float(weights[0] * deepest),
                           "other": float(field * weights[1:].mean()), "predicted": laws.even_ring_limit(cells, field, weights)})
    # On odd rings the outcome is the drop-out criterion's, read on the other cells' optimal pair at each anisotropy,
    # with the full saddle at the deepest depth beside it.
    odd = []
    for delta in config["anisotropies"]:
        for cells in [c for c in config["xxz_sizes"] if c % 2 == 1]:
            pair = _pair(_xxz(cells, 0.0, delta=delta))
            deep = programs.saddle(_xxz(cells, deepest, delta=delta))
            entry = {"delta": delta, "cells": cells, "pair": pair, "product": float(deep["weights"][0] * deepest)}
            if pair is not None:
                entry["lift"] = float(deep["level"] - pair["others"])
            odd.append(entry)
    return {"ising": ising, "xxz": xxz, "fields": fields, "odd": odd}
