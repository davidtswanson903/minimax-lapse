"""The ring's stored results carry the claims made about it."""
from repokit import sources

from saddle.experiments.ring import NO_WEIGHT


def _ring() -> dict:
    return sources.results()["ring"]


def test_uniform_ring_matches_its_ground_energy():
    ring = _ring()
    assert abs(ring["uniform"]["level"] - ring["ground_per_cell"]) < 1e-7
    assert max(abs(w - 1 / ring["cells"]) for w in ring["uniform"]["weights"]) < 1e-7


def test_deep_cell_loses_its_weight_past_a_finite_depth():
    ring = _ring()
    depth, coarse = ring["dropout"]["depth"], ring["dropout"]["coarse"]
    assert 0 < depth < ring["deep"]["depth"] and abs(depth - coarse) < 1e-3
    before = [s["weights"][0] for s in ring["sweep"] if 0 < s["depth"] < depth - 0.02]
    after = [s["weights"][0] for s in ring["sweep"] if s["depth"] > depth + 0.002]
    assert before and min(before) > 1e-4, "before the drop-out depth the cell keeps a weight"
    assert after and max(after) < 1e-5, "past it the reported weight is at the solver's noise floor"

def test_past_the_drop_out_the_deep_density_is_not_fixed():
    """By the others' optimal pair: past the drop-out the optimal states hold the deep cell anywhere from the common
    value down to the pair's least density, d - sqrt((g xi)^2 + mu^2); the state side's one state lies in that range."""
    from saddle.experiments.depths import deepest_density

    ring = _ring()
    same = next(r for r in sources.results()["depths"]["rings"] if r["cells"] == ring["cells"])
    assert abs(same["common"] - ring["deep"]["level"]) < 1e-6
    lowest = deepest_density(same["link"], same["flip"], ring["_provenance"]["config"]["field"], ring["deep"]["depth"])
    margin = ring["deep"]["level"] - ring["deep"]["densities"][0]
    assert 0 < margin < same["common"] - lowest


def test_weights_are_rates_of_the_common_value():
    for ring in _ring()["rates"].values():
        for cell in ring:
            assert abs(cell["rate"] - cell["weight"]) < 1e-4


def test_the_optimal_weighting_is_unique_to_its_spread():
    """By the rate programs: on each ring every cell's weight over the weightings near the common value lies in an
    interval of width below 5e-3 that holds both its computed weight and its rate (Proposition rates)."""
    for ring in _ring()["rates"].values():
        for cell in ring:
            assert cell["greatest"] - cell["least"] < 5e-3
            assert cell["least"] - 1e-6 <= cell["weight"] <= cell["greatest"] + 1e-6


def test_grouping_the_cells_in_pairs_changes_the_saddle():
    """By the programs on the grouped ring: past the drop-out, the pairs' densities, the cells' sums, have a common
    value other than twice the cells', and the pair holding the deep cell carries no weight while the cell beside it
    carries some."""
    from saddle import model, programs

    fine = model.ring(6, 1.0, 1.0, [1.0] + [0.0] * 5)
    cells = programs.saddle(fine)
    pairs = programs.saddle([fine[0] + fine[1], fine[2] + fine[3], fine[4] + fine[5]])
    assert abs(pairs["level"] - 2 * cells["level"]) > 1e-2
    assert pairs["weights"][0] < 1e-5 and cells["weights"][1] > 0.05
