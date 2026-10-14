"""The stored exits carry the claims about where results hold."""
from repokit import sources


def _exits() -> dict:
    return sources.results()["exits"]


def test_uniform_value_needs_the_rotation():
    rotation = _exits()["rotation"]
    assert rotation["floor"] - rotation["ground_per_cell"] > 1e-2


def test_coldest_needs_a_nondegenerate_saddle():
    assert _exits()["degenerate"]["largest_multiplier"] < 1e-9


def test_even_ring_keeps_its_cell_without_conservation():
    assert all(b["drop_depth"] is None for b in _exits()["conservation"])
