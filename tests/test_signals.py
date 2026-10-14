"""Signals under the weighted evolution, read from the stored results: the time to cross the open chain's
interior goes as the inverse of the interior's weight, and at the saddle, which leaves the interior no weight,
a disturbance never leaves the cell it starts in."""
from repokit import sources

from saddle.experiments.signals import KICKED, NO_WEIGHT


def test_crossing_time_goes_as_the_inverse_weight():
    signals = sources.results()["signals"]
    crossing = [p["far"] * p["interior_weight"] for p in signals["path"]
                if p["far"] is not None and p["interior_weight"] > NO_WEIGHT]
    assert abs(signals["exponent"] + 1.0) < 0.01, signals["exponent"]
    assert (max(crossing) - min(crossing)) / max(crossing) < 0.05, crossing


def test_the_saddle_stops_every_signal():
    signals = sources.results()["signals"]
    assert signals["interior_saddle_weight"] < NO_WEIGHT
    assert signals["reached_at_saddle"] == KICKED
    assert signals["wall"] < 1e-8
