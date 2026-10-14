"""The Ising drop-out at size by free fermions: the stored results carry the claims made about them."""
from repokit import sources


def _large() -> dict:
    return sources.results()["largedrop"]


def test_free_fermions_meet_the_sequential_quadratic_depths():
    """By the other route: the stored depths of the sequential quadratic program over the deep cell's spin-up sector."""
    stored = {r["cells"]: r for r in sources.results()["depths"]["rings"]}
    shared = [r for r in _large()["rings"] if r["cells"] in stored]
    assert shared
    for row in shared:
        other = stored[row["cells"]]
        for key in ("common", "link", "flip", "depth"):
            assert abs(row[key] - other[key]) < 1e-6, (row["cells"], key)


def test_the_closed_forms_hypotheses_hold_at_every_size():
    """The even sector's lowest level is single, the other cells' floor curves down at their optimum, so their optimal
    weighting is unique, and the bracket certifies it."""
    for row in _large()["rings"]:
        assert row["single"] > 0 and row["curvature"] < 0 and 0 <= row["gap"] < 1e-4, row["cells"]
        assert min(row["weights"]) > 0


def test_the_depth_rises_by_halving_steps():
    rings = sorted(_large()["rings"], key=lambda r: r["cells"])
    depths = [r["depth"] for r in rings]
    assert all(a < b for a, b in zip(depths, depths[1:]))
    by_size = {r["cells"]: r["depth"] for r in rings}
    steps = [by_size[2 * n] - by_size[n] for n in (40, 80, 160) if 2 * n in by_size]
    assert all(0.4 < b / a < 0.6 for a, b in zip(steps, steps[1:]))
    assert all(depths[-1] < fit["limit"] < depths[-1] + 0.01 for fit in _large()["fits"])
