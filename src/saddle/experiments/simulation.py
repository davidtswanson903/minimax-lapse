"""The readings, sourced in simulation. Idle probe qubits are attached to every place as
clocks and thermometers (Appendix H). The evolution and the equilibrium the saddle weights
generate are simulated, and each reading's quantity is measured from them: every clock's
rate from its simulated phase, every thermometer's temperature from its simulated
populations, the energy an excitation adds to the weighted generator, and the rate of the
saddle state's global phase.

Past the exact drop-out depth the deep cell's weight is zero, by Proposition 7; it is set to
zero there rather than left at the solver's noise floor."""
from __future__ import annotations

import numpy as np

from saddle import model, probes, programs

CONFIGS = {
    "golden": {"cells": 3, "field": 1.0, "coupling": 1.0, "depths": [0.0, 1.0], "omega": 1.0, "beta": 2.0,
               "steps": 8, "span": 10.0, "xxz": {"cells": 4, "delta": 0.5, "h": 0.2, "depth": 8.0}},
    "quick": {"cells": 4, "field": 1.0, "coupling": 1.0, "depths": [0.0, 0.2, 1.0], "omega": 1.0, "beta": 2.0,
              "steps": 16, "span": 20.0, "xxz": {"cells": 4, "delta": 0.5, "h": 0.2, "depth": 8.0}},
    "paper": {"cells": 6, "field": 1.0, "coupling": 1.0, "depths": [0.0, 0.2, 1.0], "omega": 1.0, "beta": 2.0,
              "steps": 64, "span": 20.0, "xxz": {"cells": 4, "delta": 0.5, "h": 0.2, "depth": 8.0}},
}
TOLERANCE = {"rtol": 1e-4, "atol": 1e-7}
NOISE = ("gaps",)
BOUND_BY = "double-precision dense linear algebra"


def _sourced(dens: list, weights: np.ndarray, level: float, config: dict, times) -> dict:
    omega, beta = config["omega"], config["beta"]
    held = programs.weighted(dens, weights)
    energies, vectors = np.linalg.eigh(held)
    ground = vectors[:, 0]
    overlap = [np.vdot(ground, (vectors * np.exp(-1j * energies * t)) @ (vectors.conj().T @ ground)) for t in times]
    phases = [probes.simulate_clock(w, omega, times) for w in weights]
    rates = [probes.fitted_rate(times, p) / omega for p in phases]
    temperatures = [probes.simulate_thermometer(w, omega, beta) for w in weights]
    redshift = [float(w * (probes.probe_piece(omega)[0, 0] - probes.probe_piece(omega)[1, 1])) / omega for w in weights]
    phase_rate = -probes.fitted_rate(times, np.unwrap(np.angle(overlap)))
    finite = [(T, w) for T, w in zip(temperatures, weights) if T is not None]
    return {
        "weights": [float(w) for w in weights], "level": float(level),
        "phases": [[float(v) for v in p] for p in phases], "clock_rates": rates,
        "temperatures": [None if T is None else float(beta * T) for T in temperatures],
        "redshift": redshift, "phase_rate": float(phase_rate),
        "gaps": {"clock": max(abs(r - probes.clock_rate(w, omega) / omega) for r, w in zip(rates, weights)),
                 "tolman": max([abs(beta * T * w - 1.0) for T, w in finite] or [0.0]),
                 "redshift": max(abs(z - probes.global_energy(w, omega) / omega) for z, w in zip(redshift, weights)),
                 "phase": abs(phase_rate - probes.phase_rate(level))},
    }


def run(config: dict) -> dict:
    cells, times = config["cells"], np.linspace(0.0, config["span"], config["steps"])
    spin = model.on_site(model.PAULI_Z, 0, cells)
    drop = programs.drop_depth(model.ring(cells, config["field"], config["coupling"]), 0, spin)
    rings = []
    for depth in config["depths"]:
        depths = [depth] + [0.0] * (cells - 1)
        dens = model.ring(cells, config["field"], config["coupling"], depths)
        found = programs.saddle(dens)
        weights = np.clip(found["weights"], 0.0, None)
        past = drop is not None and depth > drop + 1e-3
        if past:  # zero by Proposition 7, not the solver's noise floor
            weights[0] = 0.0
            weights = weights / weights.sum()
        rings.append({"depth": depth, "past_drop": past, **_sourced(dens, weights, found["level"], config, times)})
    x = config["xxz"]
    spec = {"family": "xxz", "lattice": "ring", "size": x["cells"], "J": 1.0, "delta": x["delta"], "h": x["h"], "depth": x["depth"]}
    dens = model.build(spec)
    found = programs.saddle(dens)
    xxz = {"cells": x["cells"], "depth": x["depth"], **_sourced(dens, np.clip(found["weights"], 0.0, None), found["level"], config, times)}
    worst = {k: max(r["gaps"][k] for r in rings + [xxz]) for k in ("clock", "tolman", "redshift", "phase")}
    return {"drop_depth": drop, "times": [float(t) for t in times], "rings": rings, "xxz": xxz, "gaps": worst}
