"""The saddle's own weights at size: critical transverse-field rings whose couplings vary smoothly, solved from the
weighting side, and the spatial reading tried on those weights.

On a transverse-field ring the weighted Hamiltonian is quadratic in Majorana operators whose bonds join even sites to
odd ones, so its lowest eigenvalue, the floor a weighting certifies, is minus half the sum of the singular values of
the n-by-n block C of even against odd sites (Eq. freefloor). The ring's ground state lies in the sector of even
fermion parity, with antiperiodic fermions: every coupling is positive, so by Perron-Frobenius the ground state is
unique and in that sector. The places' densities in that ground state are the floor's gradient, read from C's polar
factor. The weights that maximize the floor are found by L-BFGS over softmax weights, multiplicative steps, and
Newton's method on the full support with the Hessian applied through the polar factor's derivative; Corollary cor:gap
certifies them, the common value lying between the floor reached and the largest density. The Hessian along the
simplex, negative throughout, makes them the only maximizer; on the uniform ring of the same size it shows how flat the
floor is around uniform weights, which is why couplings varied so slightly move the weights so far.

On those weights, a long-wave packet's travel time is set against the optical distance (Eq. optical), and the saddle
state's entanglement across intervals inside the bump against the conformal vacuum's on the metric whose lapse is the
weights (Eq. interval)."""
from __future__ import annotations

import numpy as np
from scipy.optimize import minimize
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import LinearOperator, cg, expm_multiply

from repokit import workers
from repokit.equations import equation
from saddle.experiments.space import SPEED, optical_time

CONFIGS = {
    "golden": {"rings": [{"cells": 60, "eps": 1e-3}], "packet": 6, "levels": [0.9, 0.7, 0.5],
               "intervals": [[-15, 15], [-20, 10], [-10, 20]],
               "sweep": [{"cells": 40, "eps": 3e-3}, {"cells": 60, "eps": 3e-3}]},
    "quick": {"rings": [{"cells": 120, "eps": 1e-3}], "packet": 8, "levels": [0.9, 0.7, 0.5, 0.35],
              "intervals": [[-20, 20], [-30, 10], [-10, 30]],
              "sweep": [{"cells": 60, "eps": 3e-3}, {"cells": 120, "eps": 3e-4}]},
    "paper": {"rings": [{"cells": 300, "eps": 1e-5}, {"cells": 600, "eps": 1e-6}], "packet": 20,
              "levels": [0.9, 0.7, 0.5, 0.35, 0.25],
              "intervals": [[-30, 30], [-60, 20], [-20, 70], [-80, 80], [10, 90], [-100, -20]],
              "sweep": [{"cells": 150, "eps": 1e-5}, {"cells": 300, "eps": 1e-4}, {"cells": 300, "eps": 1e-6},
                        {"cells": 600, "eps": 1e-5}]},
}
# The golden ring's weights converge only to the bracket's gap there, about 1e-6, so the threads of the linear algebra
# move them in the fifth figure; the tolerance is set to that, and the gap itself is the solver's noise.
TOLERANCE = {"rtol": 1e-4, "atol": 1e-7}
RESPONSE = 3.0  # the variation over the uniform ring's flattest curvature above which the width law is fitted
NOISE = ("iterations", "gap")
BOUND_BY = "the singular value decomposition of the n-by-n block, once per step of the ascent"


def couplings(cells: int, eps: float) -> np.ndarray:
    """g_x = J_x = 1 + eps cos(2 pi x / n): critical at every place, the couplings weakest at the middle cell."""
    return 1 + eps * np.cos(2 * np.pi * np.arange(cells) / cells)


def block(c: np.ndarray, weights: np.ndarray) -> np.ndarray:
    """C, the even-odd block of the weighted ring's Majorana generator; linear in the weights."""
    n = len(weights)
    C = np.zeros((n, n))
    i = np.arange(n)
    C[i, i] = -2 * c * weights
    C[i[1:], i[:-1]] = c[:-1] * (weights[:-1] + weights[1:])
    C[0, n - 1] = -c[-1] * (weights[-1] + weights[0])
    return C


def densities_from_polar(c: np.ndarray, P: np.ndarray) -> np.ndarray:
    n = len(c)
    i = np.arange(n)
    d = c * np.diag(P)
    inter = -P[i[1:], i[:-1]]
    d[:-1] += 0.5 * c[:-1] * inter
    d[1:] += 0.5 * c[:-1] * inter
    d[-1] += 0.5 * c[-1] * P[0, n - 1]
    d[0] += 0.5 * c[-1] * P[0, n - 1]
    return d


@equation("eq:freefloor")
def free_floor(c: np.ndarray, weights: np.ndarray) -> tuple[float, np.ndarray]:
    """f(N) = -(1/2) sum_k sigma_k(C(N)), the weighted ring's lowest eigenvalue, and the densities in its ground state,
    from C's polar factor U V^T (eq:freefloor)."""
    U, s, Vt = np.linalg.svd(block(c, weights))
    return -0.5 * float(s.sum()), densities_from_polar(c, U @ Vt)


def _softmax(z):
    w = np.exp(z - z.max())
    return w / w.sum()


def _ascend(c: np.ndarray) -> np.ndarray:
    n = len(c)

    def objective(z):
        N = _softmax(z)
        f, d = free_floor(c, N)
        return -1e6 * f, -1e6 * N * (d - N @ d)

    z = minimize(objective, np.zeros(n), jac=True, method="L-BFGS-B",
                 options={"maxiter": 20000, "maxfun": 40000, "gtol": 1e-16, "ftol": 1e-18, "maxcor": 50}).x
    N = _softmax(z)
    f, d = free_floor(c, N)
    eta = 1e3
    for _ in range(4000):  # multiplicative steps move the small weights the softmax gradient barely moves
        if d.max() - f < 1e-9 or eta < 1e-6:
            break
        trial = N * np.exp(np.clip(eta * (d - f), -50, 50))
        trial /= trial.sum()
        ft, dt = free_floor(c, trial)
        if ft > f:
            N, f, d, eta = trial, ft, dt, eta * 1.5
        else:
            eta *= 0.3
    return N


def _hessian_vector(c, U, s, Vt, delta):
    Et = U.T @ block(c, delta) @ Vt.T
    X = (Et - Et.T) / (s[:, None] + s[None, :])
    return densities_from_polar(c, U @ X @ Vt)


def _newton(c: np.ndarray, N: np.ndarray, tol: float = 1e-12) -> np.ndarray:
    n = len(c)
    project = lambda v: v - v.mean()
    for _ in range(30):
        U, s, Vt = np.linalg.svd(block(c, N))
        f = -0.5 * s.sum()
        d = densities_from_polar(c, U @ Vt)
        if d.max() - f < tol:
            break
        op = LinearOperator((n, n), matvec=lambda v: -project(_hessian_vector(c, U, s, Vt, project(v))))
        step = project(cg(op, project(d), rtol=1e-10, maxiter=2000)[0])
        negative = step < 0
        alpha = 1.0 if not negative.any() else min(1.0, 0.9 * float(np.min(-N[negative] / step[negative])))
        while alpha > 1e-8:
            trial = N + alpha * step
            trial /= trial.sum()  # the floor scales with the weights, so the simplex is kept exactly
            if -0.5 * np.linalg.svd(block(c, trial), compute_uv=False).sum() >= f - 1e-15:
                N = trial
                break
            alpha *= 0.5
    return N


def saddle(c: np.ndarray) -> dict:
    N = _newton(c, _ascend(c))
    f, d = free_floor(c, N)
    return {"weights": N, "floor": f, "gap": float(d.max() - f)}


def curvature(c: np.ndarray, weights: np.ndarray) -> dict:
    """The floor's second derivatives along the simplex at the weights: the eigenvalues of its Hessian on the weightings
    that sum to zero, the Hessian built column by column from the polar factor's derivative. All negative, at weights
    that are all positive, make the weights the floor's only maximizer, since the floor is concave (Proposition
    weights). Also the share of the flattest direction that lies in the ring's longest waves, cos and sin of 2 pi x / n,
    and the least single-particle energy, whose being positive makes the floor smooth there."""
    n = len(weights)
    U, s, Vt = np.linalg.svd(block(c, weights))
    H = np.empty((n, n))
    for y in range(n):
        H[:, y] = _hessian_vector(c, U, s, Vt, np.eye(n)[y])
    basis = np.linalg.qr(np.vstack([np.ones(n), np.eye(n)[:-1]]).T)[0][:, 1:]
    values, vectors = np.linalg.eigh(basis.T @ ((H + H.T) / 2) @ basis)
    flattest = basis @ vectors[:, -2:]
    x = 2 * np.pi * np.arange(n) / n
    longest = np.linalg.qr(np.vstack([np.cos(x), np.sin(x)]).T)[0]
    return {"most": float(values[0]), "least": float(values[-1]), "single": float(s.min()),
            "longest": float(np.linalg.norm(longest.T @ flattest) ** 2 / 2)}


def _generator(c: np.ndarray, weights: np.ndarray):
    n2 = 2 * len(weights)
    t = np.empty(n2)
    t[0:n2:2] = c * weights
    t[1:n2 - 1:2] = c[:-1] * (weights[:-1] + weights[1:]) / 2
    t[n2 - 1] = c[-1] * (weights[-1] + weights[0]) / 2
    m = np.arange(n2 - 1)
    rows = np.concatenate([m + 1, m, [0, n2 - 1]])
    cols = np.concatenate([m, m + 1, [n2 - 1, 0]])
    vals = np.concatenate([2 * t[:-1], -2 * t[:-1], [-2 * t[-1], 2 * t[-1]]])
    return csr_matrix((vals, (rows, cols)), shape=(n2, n2))


def waves(c: np.ndarray, weights: np.ndarray, width: float, levels: list) -> list:
    """A packet from the bump's centre down its flank, timed between the cells where the weight falls to each level
    of its peak, against the optical time with the lapse c N scaled to one at the peak."""
    n = len(weights)
    lapse = c * weights / weights.max()
    centre = float(np.argmax(weights))
    x = np.arange(n)
    marks = [float(np.flatnonzero((x > centre) & (weights / weights.max() <= level)).min()) for level in levels]
    t_end = optical_time(lapse, centre, marks[-1] + 5) * 1.3
    position = np.arange(2 * n) / 2
    u0 = np.exp(-((position - centre) ** 2) / (2 * width ** 2))
    states = expm_multiply(_generator(c, weights / weights.max()), u0 / np.linalg.norm(u0),
                           start=0.0, stop=t_end, num=4000, endpoint=True)
    times = np.linspace(0.0, t_end, 4000)
    w = states ** 2
    centres = (w @ position) / w.sum(axis=1)
    t0 = np.interp(marks[0], centres, times)
    return [{"level": level, "cell": at, "measured": float(np.interp(at, centres, times) - t0),
             "optical": optical_time(lapse, marks[0], at)} for level, at in zip(levels[1:], marks[1:])]


@equation("eq:interval")
def interval_entropy(lapse: np.ndarray, start: int, stop: int, power: float = 1.0, wave_speed: float = SPEED) -> float:
    """S = (c/6) ln[l~^2 (v N_a)^power (v N_b)^power], c = 1/2, up to a constant: the entanglement of the cells from
    start to stop in the conformal vacuum on the metric whose lapse is N, l~ their optical length and N_a, N_b the
    lapse at the interval's two cuts; power one is a cutoff of one cell (eq:interval)."""
    length = float(np.sum(1.0 / (wave_speed * lapse[start:stop])))
    cut = lambda k: wave_speed * (lapse[k - 1] + lapse[k]) / 2
    return (1 / 12) * float(np.log(length ** 2 * (cut(start) * cut(stop)) ** power))


def _correlations(c: np.ndarray, weights: np.ndarray) -> np.ndarray:
    C = block(c, weights)
    n2 = 2 * len(weights)
    A = np.zeros((n2, n2))
    A[0::2, 1::2] = C
    A[1::2, 0::2] = -C.T
    e, W = np.linalg.eigh(1j * A)
    return (W * np.sign(e)) @ W.conj().T


def _entropy(G: np.ndarray, start: int, stop: int) -> float:
    p = np.clip((1 + np.linalg.eigvalsh(G[2 * start:2 * stop, 2 * start:2 * stop])) / 2, 1e-15, 1)
    return float(-np.sum(p * np.log(p)))


def entanglement(c: np.ndarray, weights: np.ndarray, intervals: list) -> dict:
    """The saddle state's entanglement across intervals inside the bump, with the constant from a uniform ring of
    the same size, where S = (1/6) ln[(n / pi) sin(pi l / n)] + C0."""
    n = len(weights)
    uniform = _correlations(np.ones(n), np.ones(n) / n)
    constants = [_entropy(uniform, 0, l) - (1 / 6) * np.log((n / np.pi) * np.sin(np.pi * l / n)) for l in (20, 30, 40)]
    C0 = float(np.mean(constants))
    G = _correlations(c, weights)
    lapse = c * weights / weights.max()
    centre = int(np.argmax(weights))
    rows = []
    for a, b in intervals:
        start, stop = centre + a, centre + b
        S = _entropy(G, start, stop)
        rows.append({"start": start, "stop": stop, "entropy": S,
                     "lapse_ends": [float((lapse[start - 1] + lapse[start]) / 2), float((lapse[stop - 1] + lapse[stop]) / 2)],
                     "miss": S - (interval_entropy(lapse, start, stop) + C0),
                     "optical_miss": S - (interval_entropy(lapse, start, stop, power=0.0) + C0)})
    powers = np.linspace(0.8, 1.2, 4001)
    errors = [np.mean([(r["entropy"] - interval_entropy(lapse, r["start"], r["stop"], p) - C0) ** 2 for r in rows])
              for p in powers]
    return {"constant": C0, "constant_spread": float(np.ptp(constants)), "intervals": rows,
            "power": float(powers[int(np.argmin(errors))])}


def _ring(job: tuple) -> dict:
    """One ring: its saddle, certified, with its curvatures; the reading's two tests where `full`, else the width."""
    ring, packet, levels, intervals, full = job
    n = ring["cells"]
    c = couplings(n, ring["eps"])
    found = saddle(c)
    N = found["weights"]
    half = np.flatnonzero(N > N.max() / 2)
    out = {"cells": n, "eps": ring["eps"], "floor": found["floor"], "gap": found["gap"],
           "width": int(half.max() - half.min() + 1), "uniform_curvature": curvature(np.ones(n), np.ones(n) / n)}
    if full:
        out.update({"weights": N.tolist(), "antipode": float(N[0] / N.max()), "curvature": curvature(c, N),
                    "waves": waves(c, N, packet, levels), "entanglement": entanglement(c, N, intervals)})
    return out


def run(config: dict) -> dict:
    """Every ring is independent, so they run one worker each when REPOKIT_WORKERS is set, sharing the cores'
    threads (repokit.workers); the result does not depend on how many workers run it."""
    jobs = [(ring, config["packet"], config["levels"], config["intervals"], True) for ring in config["rings"]]
    jobs += [(ring, None, None, None, False) for ring in config["sweep"]]
    done = workers.map(_ring, jobs)  # one worker per ring, the cores' threads shared among them

    keys = ("cells", "eps", "floor", "gap", "width", "weights", "antipode", "curvature", "uniform_curvature", "waves",
            "entanglement")
    rings = [{k: row[k] for k in keys} for row in done[:len(config["rings"])]]
    # The width alone across the couplings' variation and the ring's size, each saddle certified by the bracket, with
    # the variation over the uniform ring's flattest curvature: the response, which sets the regime.
    sweep = []
    for row in done:
        flat = -row["uniform_curvature"]["least"]
        sweep.append({"cells": row["cells"], "eps": row["eps"], "width": row["width"], "gap": row["gap"], "flat": flat,
                      "response": row["eps"] / flat})
    bent = [r for r in sweep if r["response"] >= RESPONSE]
    scale = np.array([(r["cells"] ** 2 / r["eps"]) ** 0.25 for r in bent])
    widths = np.array([r["width"] for r in bent], dtype=float)
    design = np.column_stack([np.ones(len(bent)), np.log([r["cells"] for r in bent]), np.log([r["eps"] for r in bent])])
    fit = np.linalg.lstsq(design, np.log(widths), rcond=None)[0]
    return {"rings": rings, "sweep": sweep, "law": {"fitted": len(bent), "cells_power": float(fit[1]),
                                                    "eps_power": float(fit[2]),
                                                    "ratio_least": float((widths / scale).min()),
                                                    "ratio_most": float((widths / scale).max())}}
