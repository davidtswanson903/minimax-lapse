"""The stored records protocol: Bell's bound factored, the gluing program and its dual, the violation tracking the
outside copy; and the history state, which holds every step at once."""
import numpy as np

from repokit import sources
from saddle.experiments import records


def _records() -> dict:
    return sources.results()["records"]


def test_the_factorization_bounds():
    stored = _records()
    assert stored["local"] == 2 and abs(stored["quantum"] - 2 * np.sqrt(2)) < 1e-12
    assert not stored["own_commute"] and stored["across_commute"]


def test_the_gluing_program_and_its_dual_agree():
    glued = _records()["gluing"]
    assert abs(glued["share"] - glued["dual"]) < 1e-6
    assert abs(1 - glued["share"] - (np.sqrt(2) - 1)) < 1e-6
    certificate = np.array(glued["certificate"])
    chsh_pattern = np.array([[0, 1, 1, 0], [0, 1, 1, 0], [0, 1, 1, 0], [1, 0, 0, 1]])
    assert np.allclose(certificate, chsh_pattern, atol=1e-5)


def test_the_records_follow_the_copy():
    stored = _records()
    for row in stored["protocol"]:
        assert abs(row["chsh"] - row["closed_form"]) < 1e-9 and abs(row["energy"]) < 1e-12 and row["signal"] < 1e-12
        if abs(row["overlap"] - stored["threshold"]) > 1e-6:
            assert (row["chsh"] > 2) == (row["overlap"] > stored["threshold"])


def test_the_history_holds_every_step():
    """Random unitary steps on two qubits: the history Hamiltonian's null space is the history state, and each step's
    state is the rest conditioned on the register's reading."""
    rng = np.random.default_rng(3)
    steps = []
    for _ in range(3):
        q, r = np.linalg.qr(rng.normal(size=(4, 4)) + 1j * rng.normal(size=(4, 4)))
        steps.append(q * (np.diag(r) / np.abs(np.diag(r))))
    first = rng.normal(size=4) + 1j * rng.normal(size=4)
    first = first / np.linalg.norm(first)
    levels, states = np.linalg.eigh(records.history_hamiltonian(steps, first))
    assert abs(levels[0]) < 1e-12 and levels[1] > 1e-3
    whole = states[:, 0].reshape(len(steps) + 1, 4)
    state = first
    for t in range(len(steps) + 1):
        if t:
            state = steps[t - 1] @ state
        rest = whole[t] / np.linalg.norm(whole[t])
        assert abs(abs(np.vdot(rest, state)) - 1) < 1e-10
