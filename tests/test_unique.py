"""Uniqueness of the optimal weighting and the coldest limit on a set: the stored results carry the claims made about
them."""
from repokit import sources

from saddle.experiments.unique import WIDE


def _unique() -> dict:
    return sources.results()["unique"]


def test_the_certificate_and_the_spreads_sort_every_model():
    stored = _unique()
    models, summary = stored["models"], stored["summary"]
    assert summary["unique"] + summary["narrow"] + summary["wide"] == summary["models"] == len(models)
    for row in models:
        if row["unique"]:
            assert row["gap"] > 1e-6 and (row["curvature"] is None or row["curvature"] < -1e-6)
        else:
            assert row["fails"] in ("lowest level not single", "flat") and row["spread"] >= 0


def test_every_open_ising_chain_and_star_has_a_set_of_optima():
    """Before any cell is deepened their ends, or their leaves, hold the common value apart from one another, so the
    weight between them is free; deepening one end of a chain leaves the other to hold it alone."""
    rows = [r for r in _unique()["models"] if r["spec"]["family"] == "ising" and r["spec"]["lattice"] in ("chain", "star")
            and "depth" not in r["spec"]]
    assert rows and all(not r["unique"] and r["spread"] > WIDE for r in rows)


def test_the_readings_ring_is_certified_before_its_drop_out():
    rings = {r["spec"]["depth"]: r for r in _unique()["rings"]}
    assert rings[0.0]["unique"] and rings[0.2]["unique"]
    assert not rings[1.0]["unique"] and rings[1.0]["fails"] == "lowest level not single"
    assert rings[1.0]["spread"] < 5e-3


def test_the_open_chains_coldest_proportions_stay_symmetric_and_close_on_the_even_split():
    steps = sources.results()["bridges"]["coldest_chain"]["steps"]
    ends = [s["ends"] for s in steps]
    assert all(abs(a - b) < 1e-6 for a, b in ends), "the reflection keeps the two ends equal"
    assert all(x[0] < y[0] < 0.5 for x, y in zip(ends, ends[1:])), "each end rises toward one half"
    distances = [s["distance"] for s in steps]
    assert all(a > b for a, b in zip(distances, distances[1:]))
