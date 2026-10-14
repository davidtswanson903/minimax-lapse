"""Maximum entropy under constraints on local energies, and its two uses here: the coldest
equilibrium a common bound on every local energy allows (Appendix I), and a fitted entanglement
Hamiltonian, which is the maximum-entropy state of its part's measured local energies."""
from __future__ import annotations

import numpy as np
from scipy.optimize import minimize

from repokit.equations import equation


def _log_partition(held: np.ndarray):
    energies, vectors = np.linalg.eigh(held)
    weights = np.exp(-(energies - energies[0]))
    return -energies[0] + np.log(weights.sum()), (vectors * (weights / weights.sum())) @ vectors.conj().T


def thermal_state(terms: list, multipliers) -> np.ndarray:
    """exp(-sum_k beta_k h_k) / Z."""
    return _log_partition(sum(float(b) * t for b, t in zip(multipliers, terms)))[1]


def _solve(objective, size: int, bounded: bool, start=None) -> np.ndarray:
    found = minimize(objective, np.zeros(size) if start is None else np.asarray(start, dtype=float), jac=True,
                     method="L-BFGS-B", bounds=[(0, None)] * size if bounded else None,
                     options={"maxiter": 50000, "ftol": 1e-15, "gtol": 1e-12})
    return found.x


def coldest_multipliers(dens: list, bound: float, start=None) -> np.ndarray:
    """The multipliers beta_x(c) of the state of greatest entropy holding every local energy at or below
    the bound c, rho_c = exp(-sum_x beta_x D_x) / Z, from the dual: minimize log Z(beta) + c sum_x beta_x over
    beta >= 0. Their proportions are coldest_shape's."""
    def dual(beta):
        logz, rho = _log_partition(sum(b * d for b, d in zip(beta, dens)))
        return logz + bound * beta.sum(), np.array([bound - np.real(np.trace(rho @ d)) for d in dens])
    return _solve(dual, len(dens), bounded=True, start=start)


@equation("eq:coldest")
def coldest_shape(multipliers) -> np.ndarray:
    """beta(u) / tau(u), with tau(u) = sum_x beta_x(u): the proportions of the coldest equilibrium's multipliers,
    which tend to the saddle weights N* as the bound u falls to the common value (eq:coldest)."""
    multipliers = np.asarray(multipliers, dtype=float)
    return multipliers / multipliers.sum()


@equation("eq:entanglement")
def entanglement_weights(terms: list, moments, start=None) -> np.ndarray:
    """The weights beta_k of rho_A = exp(-sum_k beta_k h_k) / Z whose terms have the measured expectations
    m_k: the maximum-entropy state of the part's measured local energies, from minimizing
    log Z(beta) + sum_k beta_k m_k (eq:entanglement)."""
    moments = np.asarray(moments, dtype=float)
    def dual(beta):
        logz, rho = _log_partition(sum(b * t for b, t in zip(beta, terms)))
        return logz + beta @ moments, moments - np.array([np.real(np.trace(rho @ t)) for t in terms])
    return _solve(dual, len(terms), bounded=False, start=start)
