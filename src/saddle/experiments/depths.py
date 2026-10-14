"""The Ising ring's drop-out depth from the weighting side, at sizes the state side cannot reach.

A cell deepened by mu Z drops out exactly when the other cells, optimized without it, hold it at or below their
common value e' (Proposition dropout), so its depth needs only the others' optimum, which no depth enters. On the
transverse-field ring the others' weighted Hamiltonian has no term of the deep cell's own: it conserves that cell's
spin, and the global spin flip maps the spin-up sector onto the spin-down one, so its ground space is a pair,
psi_+ and its flip. Over that pair, the least depth at which some optimal state holds the deep cell at or below e'
is sqrt((l - e')^2 - (g r)^2) (Eq. depth), with l the deep cell's link energy in psi_+ and r the others' flip in
psi_+. So the others' optimum is found in the spin-up sector alone, 2^(n-1) amplitudes, by ascending the weighted
ground energy with a sparse eigensolver, the weights taken symmetric under the ring's reflection through the deep
cell."""
from __future__ import annotations

import numpy as np
import scipy.sparse as sparse
from scipy.optimize import curve_fit, minimize
from scipy.sparse.linalg import eigsh

from repokit.equations import equation

CONFIGS = {
    "golden": {"field": 1.0, "coupling": 1.0, "sizes": [3, 4, 5, 6, 7], "fit_from": [3, 4]},
    "quick": {"field": 1.0, "coupling": 1.0, "sizes": list(range(3, 13)), "fit_from": [6, 8]},
    "paper": {"field": 1.0, "coupling": 1.0, "sizes": list(range(3, 21)), "fit_from": [8, 12]},
}
TOLERANCE = {"rtol": 1e-6, "atol": 1e-9}
NOISE = ()
BOUND_BY = "a sparse eigensolver on the deep cell's spin-up sector, 2^(n-1) amplitudes, at each step of the ascent"

_X = sparse.csr_matrix(np.array([[0.0, 1.0], [1.0, 0.0]]))
_Z = sparse.csr_matrix(np.diag([1.0, -1.0]))


def _on(single, qubit: int, count: int):
    out = sparse.identity(1, format="csr")
    for q in range(count):
        out = sparse.kron(out, single if q == qubit else sparse.identity(2, format="csr"), format="csr")
    return out


def sector(cells: int, field: float, coupling: float) -> dict:
    """The other cells' densities with the deep cell 0's spin up, its half links' energy there, and the others' flip.
    Cell c > 0 is qubit c - 1; a link to cell 0 becomes a field on its neighbour."""
    count = cells - 1
    z = [_on(_Z, q, count) for q in range(count)]
    x = [_on(_X, q, count) for q in range(count)]
    one = sparse.identity(2 ** count, format="csr")
    spin = lambda c: one if c % cells == 0 else z[c % cells - 1]
    densities = [(-field * x[c - 1] - 0.5 * coupling * (spin(c - 1) @ spin(c) + spin(c) @ spin(c + 1))).tocsr()
                 for c in range(1, cells)]
    flip = x[0]
    for q in range(1, count):
        flip = flip @ x[q]
    return {"densities": densities, "link": (-0.5 * coupling * (z[0] + z[count - 1])).tocsr(), "flip": flip.tocsr()}


def _lowest(h) -> tuple[float, np.ndarray]:
    if h.shape[0] <= 256:
        levels, states = np.linalg.eigh(h.toarray())
        return float(levels[0]), states[:, 0]
    levels, states = eigsh(h, k=1, which="SA", tol=1e-12, v0=np.ones(h.shape[0]))
    return float(levels[0]), states[:, 0]


@equation("eq:depth")
def depth_from_others(link: float, common: float, flip: float, field: float) -> float:
    """mu_c = sqrt((d - e')^2 - (g xi)^2): the least depth at which a state optimal for the others, in the pair
    psi_+ and its flip, holds the deep cell at or below their common value e' (eq:depth)."""
    return float(np.sqrt(max((link - common) ** 2 - (field * flip) ** 2, 0.0)))


def deepest_density(link: float, flip: float, field: float, depth: float) -> float:
    """d - sqrt((g xi)^2 + mu^2): the least the deep cell's density reaches over the others' optimal pair at depth
    mu, the pair's lowest level of D_0 - mu Z_0, whose diagonal is d -/+ mu and whose off-diagonal is -g xi. Past
    the drop-out every density from there up to e' is held by some optimal state."""
    return float(link - np.sqrt((field * flip) ** 2 + depth ** 2))


def others_optimum(cells: int, field: float, coupling: float) -> dict:
    """Ascend the others' weighted ground energy over weightings symmetric under the reflection through cell 0."""
    held = sector(cells, field, coupling)
    count = cells - 1
    pairs = sorted({tuple(sorted({q, count - 1 - q})) for q in range(count)})
    sizes = np.array([len(p) for p in pairs], dtype=float)

    def weights(u):
        w = np.zeros(count)
        for k, pair in enumerate(pairs):
            w[list(pair)] = u[k]
        return w

    def negative(u):
        level, state = _lowest(sum(w * d for w, d in zip(weights(u), held["densities"])))
        slope = [sum(float(state @ (held["densities"][q] @ state)) for q in pair) for pair in pairs]
        return -level, -np.array(slope)

    ascent = minimize(negative, np.ones(len(pairs)) / count, jac=True, method="SLSQP", bounds=[(0.0, 1.0)] * len(pairs),
                      constraints=[{"type": "eq", "fun": lambda u: sizes @ u - 1.0, "jac": lambda u: sizes}],
                      options={"ftol": 1e-14, "maxiter": 500})
    w = weights(ascent.x)
    common, state = _lowest(sum(wi * d for wi, d in zip(w, held["densities"])))
    return {"weights": w.tolist(), "common": common, "link": float(state @ (held["link"] @ state)),
            "flip": float(state @ (held["flip"] @ state)), "iterations": int(ascent.nit)}


def run(config: dict) -> dict:
    rings = []
    for cells in config["sizes"]:
        found = others_optimum(cells, config["field"], config["coupling"])
        rings.append({"cells": cells, **found,
                      "depth": depth_from_others(found["link"], found["common"], found["flip"], config["field"])})
    n = np.array([r["cells"] for r in rings], dtype=float)
    depth = np.array([r["depth"] for r in rings])
    fits = []
    for start in config["fit_from"]:
        keep = n >= start
        (limit, scale, power), _ = curve_fit(lambda m, a, b, p: a - b * m ** (-p), n[keep], depth[keep],
                                             p0=(depth[-1] + 0.02, 0.2, 1.0), maxfev=20000)
        fits.append({"from": int(start), "limit": float(limit), "scale": float(scale), "power": float(power)})
    return {"rings": rings, "fits": fits}
