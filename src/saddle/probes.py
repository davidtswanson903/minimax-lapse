"""Probes: idle test qubits attached to the places, as clocks and thermometers. A probe's
energy counts in its place's density and it couples to nothing else, so it reads the
evolution and the equilibrium the weights generate without moving the weights (Appendix H).
Each relation the readings use is one function here, named for its equation."""
from __future__ import annotations

import numpy as np

from repokit.equations import equation


def probe_piece(omega: float) -> np.ndarray:
    """K = (omega/2) Z, a probe qubit's own energy."""
    return 0.5 * omega * np.diag([1.0, -1.0])


@equation("eq:probes")
def probed_level(level: float, piece: np.ndarray) -> float:
    """lambda_min(A'_N) = lambda_min(A_N) + lambda_min(K): identical idle probes, one per place,
    shift the common value by their ground energy and leave every optimal weighting optimal (eq:probes)."""
    return level + float(np.linalg.eigvalsh(piece)[0])


@equation("eq:clockrate")
def clock_rate(weight: float, omega: float) -> float:
    """d phi_x / dt = N_x omega: a probe clock of frequency omega at place x advances its phase at its
    place's weight times omega per unit of the clock's reading (eq:clockrate)."""
    return weight * omega


@equation("eq:tolman")
def local_temperature(weight: float, beta: float) -> float:
    """T_x = 1 / (beta N_x): in equilibrium at global inverse temperature beta, a probe thermometer at
    place x reads this, so T_x N_x = 1/beta at every place, unbounded where N_x = 0 (eq:tolman)."""
    return float("inf") if weight <= 0 else 1.0 / (beta * weight)


@equation("eq:redshift")
def global_energy(weight: float, local: float) -> float:
    """Delta E = N_x Delta e_x: an excitation of local energy Delta e_x at place x raises the weighted
    generator's energy by its place's weight times that energy (eq:redshift)."""
    return weight * local


@equation("eq:phase")
def phase_rate(level: float) -> float:
    """The saddle state is an eigenstate of the weighted generator with eigenvalue e*, so its global
    phase turns at e* per unit of the clock's reading: e^{-i A t} psi = e^{-i e* t} psi (eq:phase)."""
    return level


def simulate_clock(weight: float, omega: float, times) -> np.ndarray:
    """A probe clock started in |+>, evolved under its weighted piece; its phase, from <X> and <Y>."""
    phases = []
    for t in times:
        state = np.exp(-1j * weight * np.diag(probe_piece(omega)) * t) * np.array([1.0, 1.0]) / np.sqrt(2)
        x = 2 * np.real(np.conj(state[0]) * state[1])
        y = 2 * np.imag(np.conj(state[0]) * state[1])
        phases.append(np.arctan2(y, x))
    return np.unwrap(np.array(phases))


def simulate_thermometer(weight: float, omega: float, beta: float) -> float | None:
    """A probe thermometer in equilibrium: its populations in the thermal state of its weighted
    piece, and the temperature at which they are thermal for its own energy; None if unbounded."""
    energies = weight * np.diag(probe_piece(omega))
    populations = np.exp(-beta * (energies - energies.min()))
    populations /= populations.sum()
    ratio = populations[1] / populations[0]  # ground over excited
    return None if ratio <= 1.0 + 1e-15 else float(omega / np.log(ratio))


def fitted_rate(times, phases) -> float:
    return float(np.polyfit(np.asarray(times, dtype=float), np.asarray(phases, dtype=float), 1)[0])
