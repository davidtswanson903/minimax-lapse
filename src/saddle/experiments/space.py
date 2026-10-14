"""The spatial part, on weightings chosen by hand: long waves and the ground state's entanglement on the weighted
critical transverse-field chain.

The weighted chain, each cell's field and each link's coupling multiplied by the weights as the weighted Hamiltonian
multiplies them, is a Majorana hopping chain over 2n sites: the bond inside cell x carries t = N_x and the bond
between cells x and x + 1 carries t = (N_x + N_{x+1}) / 2, at field equal to coupling. A Majorana amplitude evolves by
dv_m/dt = 2 t_{m-1} v_{m-1} - 2 t_m v_{m+1}, so a long wave moves at v N cells per unit time, v = 2. The experiment
sends packets of long waves along smooth weightings and times them against the optical distance, the sum of
1 / (v N) over the cells crossed (Eq. optical); follows them toward a zero of the weights, where that time grows as
ln(d1 / d2) / (v kappa); measures the share a sharp step turns back; and compares the ground state's entanglement
across every cut with the conformal vacuum's on the metric -(v N)^2 dt^2 + dx^2, with one cell as the cutoff
(Eq. vacuum)."""
from __future__ import annotations

import numpy as np
from scipy.sparse import diags
from scipy.sparse.linalg import expm_multiply

from repokit.equations import equation

CONFIGS = {
    "golden": {"cells": 600, "widths": [10, 20], "start": 200, "stop": 560, "zero_from": 400, "zero_to": [200, 100],
               "step": 0.3, "entropy_cells": 120, "times": 1500},
    "quick": {"cells": 1500, "widths": [20, 40], "start": 300, "stop": 1350, "zero_from": 1000, "zero_to": [400, 100],
              "step": 0.3, "entropy_cells": 300, "times": 3000},
    "paper": {"cells": 3000, "widths": [10, 20, 40], "start": 300, "stop": 2700, "zero_from": 2000, "zero_to": [400, 100],
              "step": 0.3, "entropy_cells": 600, "times": 6000},
}
TOLERANCE = {"rtol": 1e-3, "atol": 1e-6}
NOISE = ()
BOUND_BY = "the lattice's dispersion, which the packets' widths set"
SPEED = 2.0  # the critical chain's long waves at field and coupling one
FALL_MIN = 0.1  # the falling weighting's least weight, at the chain's far end
WELL_MIN = 0.2  # the well's least weight, at its centre
WELL_WIDTH = 0.1  # the well's width, as a share of the chain
STEPS_DROP = 0.35  # half the smoothed steps' fall, so the weight falls from one to 1 - 2 STEPS_DROP
STEPS_SMOOTH = 15.0  # the cells over which each smoothed step turns


def profiles(cells: int) -> dict:
    """The smooth weightings: uniform; falling from 1 to 0.1; a well from 1 down to 0.2 and back."""
    u = np.arange(cells) / (cells - 1)
    x = np.arange(cells, dtype=float)
    return {"uniform": np.ones(cells), "falling": FALL_MIN + (1 - FALL_MIN) * (1 - u) ** 2,
            "well": 1 - (1 - WELL_MIN) * np.exp(-((x - cells / 2) / (WELL_WIDTH * cells)) ** 2)}


def hoppings(weights: np.ndarray) -> np.ndarray:
    t = np.empty(2 * len(weights) - 1)
    t[0::2] = weights
    t[1::2] = (weights[:-1] + weights[1:]) / 2
    return t


def generator(weights: np.ndarray):
    t = hoppings(weights)
    return diags([2 * t, -2 * t], [-1, 1], format="csr")


@equation("eq:optical")
def optical_time(weights: np.ndarray, start: float, stop: float, wave_speed: float = SPEED) -> float:
    """t = sum over the cells crossed of 1 / (v N_x): the time a long wave takes from start to stop where the lapse is
    N, the metric's optical distance (eq:optical); the weights are read between cells by linear interpolation."""
    x = np.arange(len(weights))
    xs = np.linspace(start, stop, 20001)
    return float(np.trapezoid(1.0 / (wave_speed * np.interp(xs, x, weights)), xs))


def packet(cells: int, centre: float, width: float) -> np.ndarray:
    position = np.arange(2 * cells) / 2
    u = np.exp(-((position - centre) ** 2) / (2 * width ** 2))
    return u / np.linalg.norm(u)


def travel(weights: np.ndarray, centre: float, width: float, t_end: float, times: int):
    """The packet's centre, the weighted mean of its squared amplitudes over the cells, at each time."""
    states = expm_multiply(generator(weights), packet(len(weights), centre, width),
                           start=0.0, stop=t_end, num=times, endpoint=True)
    position = np.arange(2 * len(weights)) / 2
    w = states ** 2
    return np.linspace(0.0, t_end, times), (w @ position) / w.sum(axis=1)


def crossing(times: np.ndarray, centres: np.ndarray, x: float) -> float:
    k = int(np.argmax(centres >= x))
    return float(np.interp(x, centres[k - 1:k + 1], times[k - 1:k + 1])) if k > 0 else float("nan")


def waves(config: dict) -> dict:
    cells, a, b = config["cells"], float(config["start"]), float(config["stop"])
    out = {}
    for name, weights in profiles(cells).items():
        rows = []
        for width in config["widths"]:
            launch = a - 3 * width
            t_end = optical_time(weights, launch, b) * 1.05
            times, centres = travel(weights, launch, width, t_end, config["times"])
            measured = crossing(times, centres, b) - crossing(times, centres, a)
            row = {"width": width, "measured": measured, "optical": optical_time(weights, a, b)}
            if width == max(config["widths"]):
                speed = np.gradient(centres, times)
                local = np.interp(centres, np.arange(cells), weights)
                inside = (centres > a) & (centres < b)
                row["local_speed_gap"] = float(np.max(np.abs(speed[inside] / (SPEED * local[inside]) - 1)))
            rows.append(row)
        out[name] = rows
    return out


def toward_zero(config: dict) -> dict:
    """The weights falling linearly to zero at the chain's far end, N = kappa d; packets from d = zero_from."""
    cells = config["cells"]
    x0 = float(cells - 1)
    kappa = 1 / x0
    weights = np.clip((x0 - np.arange(cells)) * kappa, 1e-9, None)
    start = x0 - config["zero_from"]
    rows = []
    for width in config["widths"]:
        nearest = min(config["zero_to"])
        t_end = np.log(config["zero_from"] / (nearest * 0.5)) / (SPEED * kappa)
        times, centres = travel(weights, start - 3 * width, width, t_end, 3 * config["times"])
        t0 = crossing(times, centres, start)
        for d in config["zero_to"]:
            rows.append({"width": width, "distance": d, "measured": crossing(times, centres, x0 - d) - t0,
                         "law": float(np.log(config["zero_from"] / d) / (SPEED * kappa))})
        rows.append({"width": width, "nearest": float(x0 - centres.max())})
    return {"slope": kappa, "rows": rows}


def step(config: dict) -> dict:
    """A long wave meets a sharp step of the weights, from 1 to the step's weight and back up; the share left behind
    the step once the packet has passed, before anything reaches an end."""
    cells, low = config["cells"], config["step"]
    x = np.arange(cells, dtype=float)
    middle = cells / 2
    position = np.arange(2 * cells) / 2
    shares = {}
    width = max(config["widths"])
    start = middle - 8 * width
    for name, weights in (("down", np.where(x >= middle, low, 1.0)), ("up", np.where(x >= middle, 1.0, low))):
        t_end = optical_time(weights, start, middle) + 3 * width / (SPEED * low)
        state = expm_multiply(generator(weights), packet(cells, start, width), start=0.0, stop=t_end, num=2,
                              endpoint=True)[-1]
        shares[name] = float(np.sum(state[position < middle - 5] ** 2))
    return shares


def ground_correlations(weights: np.ndarray) -> np.ndarray:
    n2 = 2 * len(weights)
    t = hoppings(weights)
    A = np.zeros((n2, n2))
    i = np.arange(n2 - 1)
    A[i, i + 1] = -2 * t
    A[i + 1, i] = 2 * t
    e, W = np.linalg.eigh(1j * A)
    return (W * np.sign(e)) @ W.conj().T


def cut_entropies(weights: np.ndarray) -> np.ndarray:
    """The ground state's entanglement entropy of the first l cells, for every l, from its Majorana correlations."""
    G = ground_correlations(weights)
    out = np.empty(len(weights) - 1)
    for l in range(1, len(weights)):
        p = np.clip((1 + np.linalg.eigvalsh(G[:2 * l, :2 * l])) / 2, 1e-15, 1)
        out[l - 1] = float(-np.sum(p * np.log(p)))
    return out


@equation("eq:vacuum")
def vacuum_entropy(weights: np.ndarray, power: float = 1.0, wave_speed: float = SPEED) -> np.ndarray:
    """S(l) = (c/6) ln[(2 L~/pi) sin(pi x~/L~) (v N_cut)^power], c = 1/2, up to a constant: the entanglement of the
    first l cells in the conformal vacuum of an open chain on the metric -(v N)^2 dt^2 + dx^2, with x~ the optical
    position of the cut, L~ the chain's optical length and N_cut the weight at the cut; power one is a cutoff of one
    cell, power zero a cutoff of one unit of optical length (eq:vacuum)."""
    inverse = 1.0 / (wave_speed * weights)
    position = np.cumsum(inverse)[:-1]
    length = inverse.sum()
    at_cut = (weights[:-1] + weights[1:]) / 2
    return (1 / 12) * np.log((2 * length / np.pi) * np.sin(np.pi * position / length) * (wave_speed * at_cut) ** power)


def entanglement(config: dict) -> dict:
    cells = config["entropy_cells"]
    x = np.arange(cells, dtype=float)
    chosen = profiles(cells)
    chosen["steps"] = 1 - STEPS_DROP * (np.tanh((x - cells / 3) / STEPS_SMOOTH) - np.tanh((x - 2 * cells / 3) / STEPS_SMOOTH))
    bulk = slice(cells // 10, cells - cells // 10)
    uniform = cut_entropies(chosen["uniform"]) - vacuum_entropy(chosen["uniform"])
    out = {}
    for name in ("falling", "well", "steps"):
        weights = chosen[name]
        S = cut_entropies(weights)
        miss = S - (vacuum_entropy(weights) + uniform)
        optical_miss = S - (vacuum_entropy(weights, power=0.0) + uniform)
        powers = np.linspace(0.9, 1.1, 2001)
        errors = [np.mean((S[bulk] - (vacuum_entropy(weights, p)[bulk] + uniform[bulk])) ** 2) for p in powers]
        out[name] = {"miss": float(np.max(np.abs(miss[bulk]))), "optical_miss": float(np.max(np.abs(optical_miss[bulk]))),
                     "power": float(powers[int(np.argmin(errors))])}
    return out


def run(config: dict) -> dict:
    return {"waves": waves(config), "zero": toward_zero(config), "step": step(config),
            "entanglement": entanglement(config)}
