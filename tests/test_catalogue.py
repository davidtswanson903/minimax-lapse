"""The stored catalogue carries the breadth claims made about it."""
from repokit import sources


def _catalogue() -> dict:
    return sources.results()["catalogue"]


def test_general_claims_hold_across_the_catalogue():
    summary = _catalogue()["summary"]
    assert summary["programs"] == summary["forbids"] == summary["rates"] == summary["models"]
    assert summary["below_places"] > 0
    assert summary["uniform"] == summary["uniform_checked"] > 0


def test_deep_cell_drops_out_in_every_ising_model():
    rows = _catalogue()["models"]
    ising = [r for r in rows if r["spec"]["family"] == "ising"]
    assert ising and all(r["dropped"] for r in ising)
