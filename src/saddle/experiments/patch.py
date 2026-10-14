"""A part's horizon: the entanglement Hamiltonian of an interval of the infinite critical chain, and its weights.

An interval of the infinite transverse-field chain at its critical point, field equal to coupling, is in a Gaussian
state, so its entanglement Hamiltonian is quadratic in the interval's Majorana operators and follows exactly from
their correlations: K = (i/4) sum_ab kappa_ab g_a g_b, with kappa read from 2 artanh of the correlations' eigenvalues.
The eigenvalues crowd toward one, so they are taken at sixty digits, enough for intervals of up to forty cells. The interval's weights are K's coefficients in
units of the chain's own terms, at each place where a term of the chain sits: from the nearest terms alone, and with
the longer-range terms read by their first moment, which is what a term of range m contributes to the energy of a
long wave (Eq. moment). Both are set against the profile (2 pi / v) x (l - x) / l (Eq. patch), with the speed v read
from the chain's dispersion. The de Sitter temperature of the measured cosmological constant follows from the
Planck values (Eq. desitter), its uncertainty carried from theirs."""
from __future__ import annotations

import mpmath as mp
import numpy as np

from repokit.equations import equation

CONFIGS = {
    "golden": {"field": 1.0, "coupling": 1.0, "lengths": [8, 12], "digits": 60,
               "hubble": [67.36, 0.54], "fraction": [0.6847, 0.0073]},
    "quick": {"field": 1.0, "coupling": 1.0, "lengths": [8, 16, 24], "digits": 60,
              "hubble": [67.36, 0.54], "fraction": [0.6847, 0.0073]},
    "paper": {"field": 1.0, "coupling": 1.0, "lengths": [8, 16, 24, 32, 40], "digits": 60,
              "hubble": [67.36, 0.54], "fraction": [0.6847, 0.0073]},
}
TOLERANCE = {"rtol": 1e-9, "atol": 1e-12}
NOISE = ()
BOUND_BY = "the sixty-digit eigendecomposition of an interval's 2l Majorana correlations"

HBAR = 1.054571817e-34          # J s, exact in the SI since 2019 (h / 2 pi)
BOLTZMANN = 1.380649e-23        # J / K, exact
MEGAPARSEC = 3.0856775814913673e22  # m
LIGHT = 299792458.0             # m / s, exact


def speed(field: float, coupling: float, step: float = 1e-6) -> float:
    """The speed of the chain's long waves, from its dispersion 2 sqrt((J - g)^2 + 4 J g sin^2(q/2)) near q = 0,
    written with the sine so that the small difference does not cancel."""
    return float(2 * np.sqrt((coupling - field) ** 2 + 4 * coupling * field * np.sin(step / 2) ** 2) / step)


def correlation(p: int, q: int):
    """<i g_{2p} g_{2q+1}> in the infinite critical chain's ground state: -2 / (pi (2(q - p) + 1))."""
    return -mp.mpf(2) / (mp.pi * (2 * (q - p) + 1))


def generator(length: int, digits: int) -> np.ndarray:
    """kappa, real and antisymmetric, minus the coefficients k of the entanglement Hamiltonian K, with the interval's
    state proportional to exp(-K), of `length` cells of the infinite critical
    chain: the correlations <i g_{2p} g_{2q+1}> = -2 / (pi (2(q - p) + 1)), and kappa = 2 artanh of them."""
    mp.mp.dps = digits
    m = mp.matrix(2 * length, 2 * length)
    for p in range(length):
        for q in range(length):
            value = correlation(p, q)
            m[2 * p, 2 * q + 1] = mp.mpc(0, value)
            m[2 * q + 1, 2 * p] = mp.mpc(0, -value)
    levels, modes = mp.eighe(m)
    k = modes * mp.diag([2 * mp.atanh(e) for e in levels]) * modes.transpose_conj()
    return np.array([[float(mp.re(-1j * k[a, b])) for b in range(2 * length)] for a in range(2 * length)])


def nearest_weights(kappa: np.ndarray, coupling: float) -> np.ndarray:
    """The coefficient of each nearest term, in units of the chain's own, at the 2l - 1 places a term sits."""
    return np.array([kappa[a, a + 1] / (-2 * coupling) for a in range(len(kappa) - 1)])


@equation("eq:moment")
def moment_weights(kappa: np.ndarray, coupling: float) -> np.ndarray:
    """beta_a = sum over odd m of m k_{b, b+m} / (2J), computed as m kappa_{b, b+m} / (-2J) with kappa = -k, over the
    terms of every range centred where the nearest
    term (a, a+1) sits: each term weighted by what it adds to a long wave's energy (eq:moment)."""
    n = len(kappa)
    out = np.zeros(n - 1)
    for a in range(n - 1):
        for m in range(1, n, 2):
            b = a - (m - 1) // 2
            if b < 0 or b + m >= n:
                break
            out[a] += m * kappa[b, b + m] / (-2 * coupling)
    return out


@equation("eq:patch")
def profile(x: np.ndarray, length: int, wave_speed: float) -> np.ndarray:
    """(2 pi / v) x (l - x) / l: the weights of a part of a critical whole, vanishing at its edges (eq:patch)."""
    return 2 * np.pi / wave_speed * x * (length - x) / length


@equation("eq:desitter")
def desitter_temperature(hubble: float, fraction: float) -> float:
    """T = hbar H0 sqrt(Omega_Lambda) / (2 pi k_B), in kelvin: the temperature kappa / 2 pi of a horizon whose surface
    gravity is sqrt(Lambda / 3) c, with Lambda = 3 Omega_Lambda H0^2 / c^2 (eq:desitter). H0 in km/s/Mpc."""
    rate = hubble * 1e3 / MEGAPARSEC
    return HBAR * rate * np.sqrt(fraction) / (2 * np.pi * BOLTZMANN)


def _fit(weights: np.ndarray, x: np.ndarray, length: int) -> tuple[float, float]:
    shape = x * (length - x) / length
    coefficient = float(weights @ shape / (shape @ shape))
    return coefficient, float(np.max(np.abs(weights - coefficient * shape)) / np.max(weights))


def run(config: dict) -> dict:
    if config["field"] != config["coupling"]:
        raise ValueError("the correlations are the critical chain's: field equal to coupling")
    wave_speed = speed(config["field"], config["coupling"])
    intervals = []
    for length in config["lengths"]:
        kappa = generator(length, config["digits"])
        x = (np.arange(2 * length - 1) + 1) / 2
        near, full = nearest_weights(kappa, config["coupling"]), moment_weights(kappa, config["coupling"])
        near_c, near_gap = _fit(near, x, length)
        full_c, full_gap = _fit(full, x, length)
        intervals.append({"length": length, "x": x.tolist(), "near": near.tolist(), "full": full.tolist(),
                          "near_coefficient": near_c, "near_gap": near_gap,
                          "coefficient": full_c, "gap": full_gap})
    (hubble, hubble_error), (fraction, fraction_error) = config["hubble"], config["fraction"]
    temperature = desitter_temperature(hubble, fraction)
    spread = float(np.hypot(hubble_error / hubble, 0.5 * fraction_error / fraction))
    return {"speed": wave_speed, "two_pi_over_speed": float(2 * np.pi / wave_speed), "intervals": intervals,
            "desitter": {"temperature": float(temperature), "relative_error": spread}}
