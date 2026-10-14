"""The spatial part on weightings chosen by hand: the stored results carry the claims made about them."""
from repokit import sources


def _space() -> dict:
    return sources.results()["space"]


def _widest(rows: list) -> dict:
    return max(rows, key=lambda r: r["width"])


def test_long_waves_take_the_optical_time():
    for rows in _space()["waves"].values():
        row = _widest(rows)
        assert abs(row["measured"] / row["optical"] - 1) < 1e-3
        assert row["local_speed_gap"] < 5e-3


def test_toward_a_zero_the_time_grows_as_a_logarithm():
    zero = _space()["zero"]["rows"]
    widest = max(r["width"] for r in zero)
    assert all(abs(r["measured"] / r["law"] - 1) < 3e-3 for r in zero if "law" in r and r["width"] == widest)


def test_a_step_turns_back_a_share():
    shares = _space()["step"]
    assert all(0.02 < s < 0.04 for s in shares.values())
    assert abs(shares["down"] - shares["up"]) < 1e-3


def test_the_ground_state_is_the_vacuum_on_the_metric():
    for e in _space()["entanglement"].values():
        assert e["miss"] < 5e-4
        assert abs(e["power"] - 1) < 2e-3
        assert e["optical_miss"] > 0.05
