"""Time read from a timeless whole: the ring, an idle probe at every place, and a clock, solved once.

The whole is held to one condition and has no time parameter: it is stationary under the saddle's weighted
generator, with the probes attached, plus the clock's own energy. Its state is the part of a starting state that
the condition keeps: the ring at its floor, every probe in an equal superposition of its two levels, and the clock
at its first reading, projected onto the generator's levels the clock cancels to within its resolution. That one
state, conditioned on each of the clock's readings, gives the rest at that reading, and each probe's coherence
turns at some rate per unit of the reading (Proposition 9). The ring at its floor and the probes are diagonal in the
weighted generator's eigenbasis, so the whole is solved there, over the probes' levels and the clock's. As the
clock's levels grow its resolution shrinks, and the rates approach the time-stepped simulation's, its ideal-clock
limit."""
from __future__ import annotations

import itertools

import numpy as np

from repokit.equations import equation
from saddle import model, programs

CONFIGS = {
    "golden": {"cells": 3, "field": 1.0, "coupling": 1.0, "depths": [0.2, 1.0], "omega": 1.0,
               "levels": [8, 32], "readings": 32},
    "quick": {"cells": 4, "field": 1.0, "coupling": 1.0, "depths": [0.2, 1.0], "omega": 1.0,
              "levels": [8, 32, 128], "readings": 48},
    "paper": {"cells": 6, "field": 1.0, "coupling": 1.0, "depths": [0.2, 1.0], "omega": 1.0,
              "levels": [8, 32, 128, 512], "readings": 64},
}
TOLERANCE = {"rtol": 1e-4, "atol": 1e-7}
NOISE = ()
BOUND_BY = "the saddle's semidefinite program and dense linear algebra over the probes' and the clock's levels"
TURNS = 2.0  # the readings span this many turns of the fastest probe, or a quarter of the clock's period if shorter


def weights_at(config: dict, depth: float, drop: float | None) -> np.ndarray:
    """The saddle weights of the ring with its first cell at the depth; past the drop-out depth the deep cell's
    weight is zero, as in the time-stepped simulation, rather than the solver's noise."""
    cells = config["cells"]
    weights = np.array(programs.saddle(model.ring(cells, config["field"], config["coupling"],
                                                   [depth] + [0.0] * (cells - 1)))["weights"], dtype=float)
    if drop is not None and depth > drop + 1e-3:
        weights[0] = 0.0
    return weights


@equation("eq:conditioned")
def conditioned(whole: np.ndarray, clock_energy: np.ndarray, t: float) -> np.ndarray:
    """psi(t) = <tau(t)|Psi>, with |tau(t)> = exp(-i H_C t)|tau> and |tau> the clock's equal superposition of its
    levels: the rest conditioned on the clock at reading t (eq:conditioned). Psi is given by its amplitudes over
    the rest's levels and the clock's, in the basis that diagonalizes both."""
    reading = np.exp(1j * clock_energy * t) / np.sqrt(len(clock_energy))
    return whole @ reading


def conditioned_rates(weights: np.ndarray, omega: float, levels: int, readings: int) -> dict:
    """Solve the whole once and read each probe's rate per unit of the clock's reading."""
    strings = np.array(list(itertools.product((0, 1), repeat=len(weights))))
    probe_energy = omega * strings @ weights               # the probes' levels under the weighted generator
    span = float(probe_energy.max())
    resolution = span / (levels - 1)
    clock_energy = -resolution * np.arange(levels)         # the clock's levels, spanning the probes' range
    # the condition keeps the (probe level, clock level) pairs whose energies cancel to within half the resolution
    kept = np.abs(probe_energy[:, None] + clock_energy[None, :]) <= resolution / 2 + 1e-12
    whole = kept / np.sqrt(len(strings) * levels)          # the starting state's amplitudes, projected
    horizon = min(TURNS * 2 * np.pi / (omega * weights.max()), 0.25 * 2 * np.pi / resolution)
    times = np.linspace(0.0, horizon, readings)
    rates = []
    for x in range(len(weights)):
        up, down = strings[:, x] == 1, strings[:, x] == 0
        partner = [int(np.flatnonzero((strings == np.where(np.arange(len(weights)) == x, 0, s)).all(axis=1))[0])
                   for s in strings[up]]
        phases = []
        for t in times:
            rest = conditioned(whole, clock_energy, t)
            phases.append(np.angle(np.sum(rest[up] * np.conj(rest[partner]))))
        rates.append(float(-np.polyfit(times, np.unwrap(phases), 1)[0] / omega))
    return {"levels": levels, "resolution": resolution, "rates": rates,
            "gap": float(np.max(np.abs(np.array(rates) - weights)))}


def run(config: dict) -> dict:
    cells = config["cells"]
    spin = model.on_site(model.PAULI_Z, 0, cells)
    drop = programs.drop_depth(model.ring(cells, config["field"], config["coupling"]), 0, spin)
    rings = []
    for depth in config["depths"]:
        weights = weights_at(config, depth, drop)
        rings.append({"depth": depth, "weights": [float(w) for w in weights],
                      "clocks": [conditioned_rates(weights, config["omega"], levels, config["readings"])
                                 for levels in config["levels"]]})
    return {"cells": cells, "drop_depth": drop, "rings": rings}
