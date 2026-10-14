"""The Ising drop-out at size: the other cells' optimum on the transverse-field ring by free fermions.

With the deepened cell's weight at zero it carries no field of its own, so its spin is conserved. The weighted
Hamiltonian of the other cells is still quadratic in Majorana operators, the transverse-field ring's even-odd block
with the deepened cell's diagonal entry zero, and its lowest level is the pair of Eq. depth: the ground state of the
even fermion parity sector and its image under the deepened spin, which carries one sector to the other at the same
energy. So the floor is the even sector's free-fermion ground energy (Eq. freefloor), and the pair's two numbers come
from that sector's correlations: d, the deepened cell's half links, and xi = <X_0>, since the pair state of spin up is
(1 + Z_0) times the even ground state over the square root of two. The other cells' optimal weighting is found by
ascent and Newton's method on the weightings that give the deepened cell none, certified by the bracket, and checked
unique by the floor's Hessian there; Eq. depth then gives the drop-out depth at hundreds of cells."""
from __future__ import annotations

import numpy as np
from scipy.optimize import minimize
from scipy.sparse.linalg import LinearOperator, cg

from saddle.experiments.bump import _hessian_vector, block, densities_from_polar
from saddle.experiments.depths import depth_from_others

CONFIGS = {
    # Every fit takes three sizes at least, one for each of its three parameters, so it is determined.
    "golden": {"sizes": [6, 8, 10], "fit_from": [6]},
    "quick": {"sizes": [6, 10, 20, 30, 40], "fit_from": [10, 20]},
    "paper": {"sizes": [6, 8, 10, 12, 14, 16, 18, 20, 30, 40, 60, 80, 120, 160, 240, 320],
              "fit_from": [40, 80]},
}
# The ring's optimum is an iterate the optimizer stops at a gap near 1e-7, so its depths and weights agree across
# platforms to about 1e-5; the gap itself is the solver's residual. A fit picks its power on a grid of step 1e-3,
# and its limit and scale move with that choice.
TOLERANCE = {"rtol": 1e-4, "atol": 1e-8}
NOISE = ("gap",)
TOLERANCE_BY_KEY = {"fits": {"rtol": 1e-2, "atol": 1e-6}}
BOUND_BY = "the singular value decomposition of the n-by-n block, once per step"
FIELD = COUPLING = 1.0  # the ring's g and J, the critical point


def _floor(weights: np.ndarray):
    """The even sector's floor and its polar factor at a weighting of the whole ring, the deepened cell's included."""
    U, s, Vt = np.linalg.svd(block(np.full(len(weights), FIELD), weights))
    return -0.5 * float(s.sum()), U, s, Vt


def _others(cells: int, start: np.ndarray | None = None) -> np.ndarray:
    """The other cells' optimal weighting, the deepened cell 0 given none: L-BFGS over softmax weights, from a smaller
    ring's profile where one is given, then Newton's method on the condition that every other cell's density sits at
    the floor."""
    c = np.full(cells, FIELD)

    def full(z):
        w = np.exp(z - z.max())
        return np.concatenate([[0.0], w / w.sum()])

    def objective(z):
        N = full(z)
        f, U, s, Vt = _floor(N)
        d = densities_from_polar(c, U @ Vt)[1:]
        return -f, -N[1:] * (d - N[1:] @ d)

    z0 = np.zeros(cells - 1)
    if start is not None:  # the smaller ring's profile, read at the same fractions of the way round
        profile = np.interp(np.arange(1, cells) / cells, np.arange(1, len(start) + 1) / (len(start) + 1), start)
        z0 = np.log(profile)
    N = full(minimize(objective, z0, jac=True, method="L-BFGS-B",
                      options={"maxiter": 5000, "gtol": 1e-14, "ftol": 1e-16}).x)
    project = lambda v: v - v.mean()
    for _ in range(30):
        f, U, s, Vt = _floor(N)
        d = densities_from_polar(c, U @ Vt)[1:]
        if d.max() - f < 1e-13:
            break
        hvp = lambda v: project(_hessian_vector(c, U, s, Vt, np.concatenate([[0.0], project(v)]))[1:])
        with np.errstate(divide="ignore", invalid="ignore"):
            step = project(cg(LinearOperator((cells - 1, cells - 1), matvec=lambda v: -hvp(v)), project(d),
                              rtol=1e-12, maxiter=4000)[0])
        if not np.all(np.isfinite(step)):  # the gradient already lies at the solver's floor: nothing left to take
            break
        alpha, negative = 1.0, step < 0
        if negative.any():
            alpha = min(1.0, 0.9 * float(np.min(-N[1:][negative] / step[negative])))
        while alpha > 1e-10:
            trial = N.copy()
            trial[1:] = N[1:] + alpha * step
            trial /= trial.sum()
            if _floor(trial)[0] >= f - 1e-15:
                N = trial
                break
            alpha *= 0.5
    return N


def _curvature(N: np.ndarray, U, s, Vt) -> float:
    """The largest eigenvalue of the floor's Hessian on the weightings that leave the deepened cell at zero and sum to
    zero: negative, the other cells' optimal weighting is their floor's only maximizer."""
    n = len(N)
    c = np.full(n, FIELD)
    H = np.empty((n - 1, n - 1))
    for k in range(n - 1):
        e = np.zeros(n)
        e[k + 1] = 1.0
        H[:, k] = _hessian_vector(c, U, s, Vt, e)[1:]
    basis = np.linalg.qr(np.vstack([np.ones(n - 1), np.eye(n - 1)[:-1]]).T)[0][:, 1:]
    return float(np.linalg.eigvalsh(basis.T @ ((H + H.T) / 2) @ basis).max())


def ring(cells: int, start: np.ndarray | None = None) -> dict:
    N = _others(cells, start)
    f, U, s, Vt = _floor(N)
    P = U @ Vt
    densities = densities_from_polar(np.full(cells, FIELD), P)
    # <Z_x Z_{x+1}> = P[x+1, x] in the bulk and <Z_{n-1} Z_0> = -P[0, n-1] across the closing bond; <X_x> = -P[x, x].
    link = -0.5 * COUPLING * (-P[0, cells - 1] + P[1, 0])
    flip = -float(P[0, 0])
    return {"cells": cells, "common": f, "link": float(link), "flip": flip,
            "depth": depth_from_others(float(link), f, flip, FIELD), "gap": float(densities[1:].max() - f),
            "single": float(s.min()), "curvature": _curvature(N, U, s, Vt), "weights": N[1:].tolist()}


def _fit(rows: list, start: int) -> dict:
    """mu_c(n) = mu_inf - a n^-p over the sizes from `start`, by least squares in p on a grid and linearly in the rest."""
    use = [r for r in rows if r["cells"] >= start]
    n = np.array([r["cells"] for r in use], dtype=float)
    mu = np.array([r["depth"] for r in use])
    best = None
    for p in np.linspace(0.2, 3.0, 2801):
        design = np.column_stack([np.ones(len(n)), -n ** -p])
        coef, *_ = np.linalg.lstsq(design, mu, rcond=None)
        error = float(np.sum((design @ coef - mu) ** 2))
        if best is None or error < best[0]:
            best = (error, p, coef)
    return {"from": start, "limit": float(best[2][0]), "scale": float(best[2][1]), "power": float(best[1]),
            "residual": best[0]}


def run(config: dict) -> dict:
    rings, start = [], None
    for n in sorted(config["sizes"]):
        rings.append(ring(n, start))
        start = np.array(rings[-1]["weights"])
    return {"rings": rings, "fits": [_fit(rings, start) for start in config["fit_from"]]}
