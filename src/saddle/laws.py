"""The even XXZ rings' law: as the depth grows, a deepened cell's weight times its depth tends
to the field times the other cells' weights, averaged with how strongly each one's spin moves
their total (Appendix E)."""
from __future__ import annotations

import numpy as np

from repokit.equations import equation
from saddle import model


def response(cells: int, weights, J: float = 1.0, delta: float = 0.5) -> np.ndarray:
    """w_y = sum over x != 0 of chi_xy, for y != 0, where chi_xy = 2 sum_{m != 0} <0|Z_x|m><m|Z_y|0> /
    (E_m - E_0) are the static susceptibilities of the weighted links' ground state |0> in the
    zero-magnetization sector, with the field switched off (eq:law)."""
    X, Y, Z = model.PAULI_X, model.PAULI_Y, model.PAULI_Z
    links = model.lattice("ring", cells)
    bond = lambda a, b: (J * (model.two_site(X, X, a, b, cells) + model.two_site(Y, Y, a, b, cells))
                         + delta * model.two_site(Z, Z, a, b, cells))
    held = sum(0.5 * (weights[a] + weights[b]) * bond(a, b) for a, b in links)
    spins = [model.on_site(Z, x, cells) for x in range(cells)]
    sector = np.where(np.isclose(np.diag(sum(spins)), 0))[0]
    energies, vectors = np.linalg.eigh(held[np.ix_(sector, sector)])
    ground = vectors[:, 0]
    amplitudes = [vectors.T @ (s[np.ix_(sector, sector)] @ ground) for s in spins]
    chi = np.array([[2 * np.sum(amplitudes[x][1:] * amplitudes[y][1:] / (energies[1:] - energies[0]))
                     for y in range(1, cells)] for x in range(1, cells)])
    return chi.sum(axis=0)


@equation("eq:law")
def even_ring_limit(cells: int, field: float, weights, J: float = 1.0, delta: float = 0.5) -> float:
    """lim mu N_0 = b sum_y w_y N_y / sum_y w_y over the other cells y, to first order in the field b (eq:law)."""
    weights = np.asarray(weights, dtype=float)
    moves = response(cells, weights, J, delta)
    return float(field * (moves @ weights[1:]) / moves.sum())
