"""The saddle's own weights at size: the stored results carry the claims made about them."""
import numpy as np

from repokit import sources


def _rings() -> list:
    return sources.results()["bump"]["rings"]


def test_the_saddle_is_certified_by_the_bracket():
    for ring in _rings():
        assert 0 <= ring["gap"] < 1e-11


def test_the_weights_form_one_bump():
    for ring in _rings():
        weights = np.array(ring["weights"]) / max(ring["weights"])
        peak = int(np.argmax(weights))
        assert np.all(np.diff(weights[:peak + 1]) >= -1e-9) and np.all(np.diff(weights[peak:]) <= 1e-9)
        assert ring["antipode"] < 1e-3


def test_the_weights_are_the_floors_only_maximizer():
    """At the bump's weights, all positive, the floor curves down in every direction along the simplex and is smooth
    there, so the weights are its only maximizer; on the uniform ring of the same size the floor is far flatter around
    uniform weights, flattest along the ring's longest waves."""
    for ring in _rings():
        assert min(ring["weights"]) > 0
        bent, flat = ring["curvature"], ring["uniform_curvature"]
        assert bent["most"] <= bent["least"] < 0 and bent["single"] > 0
        assert flat["least"] < 0 and -flat["least"] < 0.01 * -bent["least"]
        assert flat["longest"] > 0.99


def test_the_width_grows_as_the_fourth_root_of_size_squared_over_variation():
    """Where the variation is several times the uniform ring's flattest curvature; where it is not, the weights are
    near uniform plus the ring's longest wave, and the bump spans about half the ring."""
    from saddle.experiments.bump import RESPONSE

    stored = sources.results()["bump"]
    assert all(0 <= r["gap"] < 1e-6 for r in stored["sweep"])
    law = stored["law"]
    assert law["fitted"] >= 3 and 0.4 < law["cells_power"] < 0.6 and -0.3 < law["eps_power"] < -0.2
    assert law["ratio_most"] / law["ratio_least"] < 1.2
    for row in stored["sweep"]:
        if row["response"] < RESPONSE:
            assert 0.35 < row["width"] / row["cells"] < 0.55


def test_long_waves_take_the_optical_time_on_the_bump():
    for w in _rings()[-1]["waves"]:
        assert abs(w["measured"] / w["optical"] - 1) < 3e-3


def test_the_saddle_state_is_the_vacuum_on_the_metric():
    for ring in _rings():
        rows = ring["entanglement"]["intervals"]
        assert max(abs(r["miss"]) for r in rows) < 1e-4
        assert max(abs(r["optical_miss"]) for r in rows) > 0.05
        assert abs(ring["entanglement"]["power"] - 1) < 1e-3
