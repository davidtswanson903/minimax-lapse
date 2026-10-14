"""Time read from a timeless whole, from the stored results: conditioned on its clock's reading, every probe
turns at its place's weight to within half the clock's resolution, the gap shrinking as the clock's levels grow,
and past the deep cell's drop-out its probe does not turn."""
from repokit import sources


def test_conditioned_probes_turn_at_their_weights():
    """By comparing each probe's rate, read from the whole's state conditioned on the clock, with the saddle's weight
    at its place: within half the clock's resolution at every clock size, and closer at the finest than the coarsest."""
    timeless = sources.results()["timeless"]
    for ring in timeless["rings"]:
        for clock in ring["clocks"]:
            assert clock["gap"] <= clock["resolution"] / 2 + 1e-9, (ring["depth"], clock["levels"], clock["gap"])
        assert ring["clocks"][-1]["gap"] < ring["clocks"][0]["gap"]


def test_the_stopped_probe_does_not_turn():
    """By reading the deep cell's probe past its drop-out, where its weight is zero: its rate is zero."""
    timeless = sources.results()["timeless"]
    past = [r for r in timeless["rings"] if r["weights"][0] == 0.0]
    assert past, "no depth past the drop-out"
    assert all(abs(c["rates"][0]) < 1e-12 for r in past for c in r["clocks"])
