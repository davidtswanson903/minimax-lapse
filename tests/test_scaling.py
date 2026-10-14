"""The stored scaling study carries the claims made about size."""
from repokit import sources
from repokit.equations import checks


def _scaling() -> dict:
    return sources.results()["scaling"]


def test_ising_threshold_exists_at_every_size():
    rows = _scaling()["ising"]
    assert rows and all(0 < r["depth"] < 1 and abs(r["depth"] - r["coarse"]) < 1e-3 for r in rows)
    assert all(a["depth"] < b["depth"] for a, b in zip(rows, rows[1:]))


def test_xxz_outcome_follows_parity():
    for row in _scaling()["xxz"]:
        flags = [p["dropped"] for p in row["curve"]]
        if row["cells"] % 2:
            fall = flags.index(True)
            assert all(flags[fall:]) and not any(flags[:fall])
            assert row["curve"][fall - 1]["depth"] <= row["exact"] <= row["curve"][fall]["depth"], "the exact depth lies in the sampled fall"
            weights = [p["weight"] for p in row["curve"][:fall]]
            assert weights[-1] > weights[0], "on odd rings the weight grows before it falls"
        else:
            assert row["exact"] is None and not any(flags)
            products = [p["weight"] * p["depth"] for p in row["curve"] if p["depth"] >= 2]
            assert max(products) / min(products) < 1.2, "on even rings the weight falls as the inverse of the depth"


def test_odd_rings_follow_the_others_optimal_states():
    """On odd XXZ rings the outcome is the drop-out criterion's: where some optimal state of the others gives the
    deepened spin a positive expectation the cell drops out, by the deepest depth; where none does it keeps its
    weight, which falls as the inverse of the depth, with the common value above the others'. At the curves'
    anisotropy the pair's depth meets the drop-out program's."""
    scaling = _scaling()
    exact = {r["cells"]: r["exact"] for r in scaling["xxz"]}
    deepest = scaling["_provenance"]["config"]["xxz_depths"][-1]
    for entry in scaling["odd"]:
        pair = entry["pair"]
        assert pair is not None, (entry["delta"], entry["cells"], "the others' ground level is not one pair")
        if pair["depth"] is None:
            assert max(pair["spin"]) < 0 and entry["lift"] > 1e-5 and entry["product"] > 1e-4
        else:
            assert max(pair["spin"]) > 0 and pair["depth"] < deepest and entry["product"] < 1e-5
        if entry["delta"] == 0.5:
            assert abs(pair["depth"] - exact[entry["cells"]]) < 1e-3
    assert any(e["pair"]["depth"] is None for e in scaling["odd"]), "one odd ring keeps its deepened cell"


def test_even_xxz_product_follows_the_field():
    ratios: dict = {}
    for f in _scaling()["fields"]:
        ratios.setdefault(f["cells"], []).append(f["product"] / f["field"])
    assert ratios and all((max(r) - min(r)) / (sum(r) / len(r)) < 0.1 for r in ratios.values())


@checks("eq:law")
def test_even_ring_law_matches_its_limit():
    """By solving the saddle at the deepest depth and multiplying the weight by the depth: the product meets the derived limit, closer as the field weakens."""
    by_size: dict = {}
    for f in _scaling()["fields"]:
        by_size.setdefault(f["cells"], []).append((f["field"], abs(f["product"] - f["predicted"]) / f["predicted"]))
    for gaps in by_size.values():
        gaps.sort()
        assert gaps[0][1] < 0.02, "at the weakest field the product meets its derived limit"
        assert all(a[1] < b[1] for a, b in zip(gaps, gaps[1:])), "the gap shrinks with the field"


def test_the_laws_gap_at_the_deepest_depth_is_mostly_the_depths_factor():
    """By dividing the product at depth mu by mu / (mu + b): on the smallest even ring it then meets the derived limit
    at every field, and on the largest what is left is small and grows with the field."""
    scaling = _scaling()
    deepest = max(scaling["_provenance"]["config"]["xxz_depths"])
    left: dict = {}
    for f in scaling["fields"]:
        left.setdefault(f["cells"], []).append(
            (f["field"], abs(f["product"] * (deepest + f["field"]) / deepest / f["predicted"] - 1)))
    smallest, largest = sorted(left[min(left)]), sorted(left[max(left)])
    assert all(gap < 1e-4 for _, gap in smallest), "on the smallest even ring the gap is the depth's factor alone"
    assert all(gap < 0.01 for _, gap in largest), "on the largest little is left"
    assert all(a[1] < b[1] for a, b in zip(largest, largest[1:])), "and what is left grows with the field"


def test_parity_hypotheses_hold_on_the_even_rings():
    """On every even ring of the study of size, at each field: the other cells' weighted Hamiltonian has a unique
    ground state, with M = 0; the weighted links alone have a unique ground state in that sector; and the
    response-weighted sum of the other cells' weights is positive (Proposition parity's hypotheses)."""
    import numpy as np

    from saddle import model, programs
    from saddle.experiments import scaling

    config = scaling.CONFIGS["survey"]
    for cells in [n for n in config["xxz_sizes"] if n % 2 == 0]:
        spins = [model.on_site(model.PAULI_Z, x, cells) for x in range(cells)]
        sector = np.flatnonzero(np.isclose(np.diag(sum(spins)).real, 0.0))
        for field in config["xxz_fields"]:
            others = scaling._xxz(cells, 0.0, field)[1:]
            weights = programs.saddle(others)["weights"]
            values, vectors = np.linalg.eigh(programs.weighted(others, weights))
            assert values[1] - values[0] > 1e-6, (cells, field, "the others' optimal state is not unique")
            assert abs(np.real(vectors[:, 0].conj() @ sum(spins) @ vectors[:, 0])) < 1e-8, (cells, field)
            links = programs.weighted(scaling._xxz(cells, 0.0, 0.0)[1:], weights)[np.ix_(sector, sector)]
            levels, states = np.linalg.eigh(links)
            assert levels[1] - levels[0] > 1e-6, (cells, field, "the links' ground state is not unique")
            amplitudes = np.array([states.conj().T @ z[np.ix_(sector, sector)] @ states[:, 0] for z in spins])
            chi = 2 * np.real(np.einsum("xm,ym,m->xy", amplitudes[:, 1:].conj(), amplitudes[:, 1:],
                                        1 / (levels[1:] - levels[0])))
            response = chi[1:, 1:].sum(axis=0)
            assert float(response @ np.asarray(weights)) > 0, (cells, field)


def test_the_four_cell_others_floor_is_flat_without_the_field():
    """By sampling weightings: with the field off, the three other cells of the even XXZ ring of four share a ground
    state at one level, so every weighting of them certifies the same floor; the limits in Proposition parity's law
    need not be their own optimal weighting."""
    import numpy as np

    from saddle import programs
    from saddle.experiments import scaling

    others = scaling._xxz(4, 0.0, 0.0)[1:]
    rng = np.random.default_rng(3)
    floors = [np.linalg.eigvalsh(programs.weighted(others, rng.dirichlet(np.ones(3))))[0] for _ in range(20)]
    assert max(floors) - min(floors) < 1e-9
