"""The common value's two parts on uniform critical rings, and the central charge read two ways.

On a ring the shift around it maps into itself, the uniform weighting is optimal and the common value is the ground
energy per cell (Proposition 2). At the critical point, field equal to coupling, the transverse-field ring is free
fermions: after the Jordan-Wigner map its Hamiltonian is quadratic in Majorana operators, antiperiodic across the
ring's seam in the sector that holds the ground state, and is solved exactly at sizes dense algebra cannot reach.
The common value is split into a part every place carries, the same at every size, and a part set by the ring's
size, set against the conformal term -pi c v / (6 n^2) (Eq. casimir). The speed v is read from the ring's lowest
levels; the central charge c from the size-set part and, independently, from the entanglement entropy of an
interval of the largest ring (Eq. chord)."""
from __future__ import annotations

import numpy as np

from repokit.equations import equation

CONFIGS = {
    "golden": {"field": 1.0, "coupling": 1.0, "sizes": [8, 12, 16, 20], "shortest": 3},
    "quick": {"field": 1.0, "coupling": 1.0, "sizes": [8, 12, 16, 20, 24, 28, 32], "shortest": 3},
    "paper": {"field": 1.0, "coupling": 1.0, "sizes": [8, 12, 16, 20, 24, 28, 32, 36, 40, 44, 48, 52, 56, 60],
              "shortest": 3},
}
TOLERANCE = {"rtol": 1e-9, "atol": 1e-12}
NOISE = ()
BOUND_BY = "dense algebra over 2n Majorana modes"


def majorana_generator(cells: int, field: float, coupling: float) -> np.ndarray:
    """The real antisymmetric A with H = (i/4) sum_ab A_ab g_a g_b for the ring's Hamiltonian
    -g sum_x X_x - J sum_x Z_x Z_{x+1}, in the sector of even parity: X_x = i g_{2x} g_{2x+1},
    Z_x Z_{x+1} = i g_{2x+1} g_{2x+2}, and across the seam the sign the sector's parity sets."""
    a = np.zeros((2 * cells, 2 * cells))
    for x in range(cells):
        a[2 * x, 2 * x + 1] = -2 * field
    for x in range(cells - 1):
        a[2 * x + 1, 2 * x + 2] = -2 * coupling
    a[2 * cells - 1, 0] = 2 * coupling
    return a - a.T


def solve(cells: int, field: float, coupling: float) -> dict:
    """The ground energy, the single-particle levels, and the ground state's Majorana correlations sign(iA)."""
    levels, modes = np.linalg.eigh(1j * majorana_generator(cells, field, coupling))
    return {"energy": float(-0.25 * np.abs(levels).sum()),
            "levels": np.sort(levels[levels > 0]),
            "correlations": (modes * np.sign(levels)) @ modes.conj().T}


def interval_entropy(correlations: np.ndarray, length: int) -> float:
    """The entanglement entropy of the first `length` cells, from the eigenvalues +-nu of their Majorana
    correlations: S = sum over pairs of the binary entropy of (1 + nu)/2."""
    block = correlations[: 2 * length, : 2 * length]
    nu = np.clip(np.abs(np.linalg.eigvalsh(block)), 0.0, 1.0)
    p = np.clip((1 + nu) / 2, 1e-300, 1.0)
    q = np.clip((1 - nu) / 2, 1e-300, 1.0)
    return float(0.5 * np.sum(-p * np.log(p) - q * np.log(q)))


@equation("eq:casimir")
def size_set_part(cells: int, charge: float, speed: float) -> float:
    """The part of the common value set by the ring's size, to leading order: -pi c v / (6 n^2) (eq:casimir)."""
    return -np.pi * charge * speed / (6 * cells ** 2)


@equation("eq:chord")
def chord(cells: int, length: int) -> float:
    """The coefficient of c in an interval's entanglement entropy on a ring:
    S = (c/3) ln[(n/pi) sin(pi l/n)] + const (eq:chord)."""
    return float(np.log(cells / np.pi * np.sin(np.pi * length / cells)) / 3)


def run(config: dict) -> dict:
    sizes = config["sizes"]
    rings = []
    for cells in sizes:
        solved = solve(cells, config["field"], config["coupling"])
        rings.append({"cells": cells, "common": solved["energy"] / cells,
                      "speed": float(solved["levels"][0] * cells / np.pi)})
    n = np.array(sizes, dtype=float)
    common = np.array([r["common"] for r in rings])
    # the common value: a per-place part, the size-set part's 1/n^2 term and the next order
    per_place, size_coefficient, _ = np.linalg.lstsq(np.c_[np.ones_like(n), n ** -2, n ** -4], common, rcond=None)[0]
    # the speed of the lowest excitations, read from each ring's lowest level, to the same order
    speed = float(np.linalg.lstsq(np.c_[np.ones_like(n), n ** -2], [r["speed"] for r in rings], rcond=None)[0][0])
    charge_size = float(size_coefficient / size_set_part(1, 1.0, speed))

    largest = sizes[-1]
    correlations = solve(largest, config["field"], config["coupling"])["correlations"]
    lengths = list(range(config["shortest"], largest // 2 + 1))
    entropies = [interval_entropy(correlations, length) for length in lengths]
    slope = np.polyfit([chord(largest, length) for length in lengths], entropies, 1)[0]
    return {"rings": rings, "per_place": float(per_place), "size_coefficient": float(size_coefficient),
            "speed": speed, "charge_size": charge_size,
            "entropy": {"cells": largest, "lengths": lengths, "entropies": entropies, "charge": float(slope)}}
