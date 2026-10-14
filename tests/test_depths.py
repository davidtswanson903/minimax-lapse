"""The stored Ising drop-out depths from the weighting side: they match the state side where both reach, and they
rise and level off."""
import numpy as np

from repokit import sources


def _depths() -> dict:
    return sources.results()["depths"]


def test_the_closed_form_matches_the_state_side():
    weighting = {r["cells"]: r["depth"] for r in _depths()["rings"]}
    for row in sources.results()["scaling"]["ising"]:
        assert abs(weighting[row["cells"]] - row["depth"]) < 1e-3


def test_the_depth_rises_and_levels_off():
    stored = _depths()
    depth = np.array([r["depth"] for r in stored["rings"]])
    steps = np.diff(depth)
    assert (steps > 0).all() and (np.diff(steps) < 0).all()
    early, late = stored["fits"]
    assert abs(early["limit"] - late["limit"]) < 2e-3 and early["limit"] > depth[-1]
