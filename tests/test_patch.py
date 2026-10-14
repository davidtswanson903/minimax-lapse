"""The stored static patch: an interval's entanglement-Hamiltonian weights take the profile with their longer-range
terms, and the correlations they start from are the infinite chain's."""
import numpy as np

from repokit import sources
from saddle.experiments import critical, patch


def _patch() -> dict:
    return sources.results()["patch"]


def test_the_correlations_are_the_long_rings():
    """An interval's entropy from the infinite chain's correlations agrees with a long ring's."""
    length = 6
    m = np.zeros((2 * length, 2 * length))
    for p in range(length):
        for q in range(length):
            m[2 * p, 2 * q + 1] = float(patch.correlation(p, q))
            m[2 * q + 1, 2 * p] = -m[2 * p, 2 * q + 1]
    nu = np.clip(np.abs(np.linalg.eigvalsh(1j * m)), 0.0, 1.0 - 1e-15)
    p_, q_ = (1 + nu) / 2, (1 - nu) / 2
    infinite = float(0.5 * np.sum(-p_ * np.log(p_) - q_ * np.log(q_)))
    ring = critical.interval_entropy(critical.solve(400, 1.0, 1.0)["correlations"], length)
    assert abs(infinite - ring) < 1e-4


def test_the_weights_take_the_profile_with_their_longer_range_terms():
    stored = _patch()
    two_pi = stored["two_pi_over_speed"]
    last = stored["intervals"][-1]
    assert abs(last["coefficient"] / two_pi - 1) < 1e-4 and last["gap"] < 1e-4
    assert last["near_coefficient"] / two_pi > 1.04 and last["near_gap"] > 0.01
    gaps = [abs(i["coefficient"] / two_pi - 1) for i in stored["intervals"]]
    assert all(a > b for a, b in zip(gaps, gaps[1:])), "the coefficient closes on 2 pi / v as the interval grows"
