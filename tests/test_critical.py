"""The stored critical rings carry the common value's two parts and the central charge read two ways."""
import numpy as np

from repokit import sources
from saddle import model, programs
from saddle.experiments import critical


def _critical() -> dict:
    return sources.results()["critical"]


def test_the_fermions_are_the_ring():
    """The free-fermion solution is the paper's ring: its ground energy and an interval's entropy agree with dense
    diagonalization of the densities' sum, and its common value with the uniform ring's."""
    cells, length = 8, 3
    values, vectors = np.linalg.eigh(sum(model.ring(cells, 1.0, 1.0)))
    solved = critical.solve(cells, 1.0, 1.0)
    assert abs(solved["energy"] - values[0]) < 1e-9
    p = np.linalg.svd(vectors[:, 0].reshape(2 ** length, -1), compute_uv=False) ** 2
    p = p[p > 1e-15]
    assert abs(critical.interval_entropy(solved["correlations"], length) - float(-(p * np.log(p)).sum())) < 1e-9
    six = model.ring(6, 1.0, 1.0)
    assert abs(critical.solve(6, 1.0, 1.0)["energy"] / 6 - programs.uniform_value(six)) < 1e-9


def test_the_two_central_charges_agree():
    stored = _critical()
    assert abs(stored["charge_size"] - stored["entropy"]["charge"]) < 1e-3
    common = [r["common"] for r in stored["rings"]]
    assert all(a < b for a, b in zip(common, common[1:])), "the common value rises toward its per-place part"
    assert all(c < stored["per_place"] for c in common)
