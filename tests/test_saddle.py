"""The saddle on random instances: both programs agree, the four conditions hold,
and a place below the common value carries no weight."""
import numpy as np

from saddle import model, programs

SEEDS = range(5)
TERMS = (model.PAULI_X, model.PAULI_Z)


def random_instance(seed: int, cells: int = 3):
    """Random real pieces on a ring: the densities, and the Hamiltonian built apart."""
    rng = np.random.default_rng(seed)
    links = model.ring_links(cells)
    places = [sum(rng.normal() * model.on_site(t, x, cells) for t in TERMS) for x in range(cells)]
    bonds = [sum(rng.normal() * model.on_site(s, a, cells) @ model.on_site(t, b, cells) for s in TERMS for t in TERMS)
             for a, b in links]
    return model.densities(places, bonds, links), sum(places) + sum(bonds)


def test_densities_sum_to_the_hamiltonian():
    for seed in SEEDS:
        dens, hamiltonian = random_instance(seed)
        assert np.allclose(sum(dens), hamiltonian)


def test_programs_agree_on_random_instances():
    for seed in SEEDS:
        assert programs.saddle(random_instance(seed)[0])["gap"] < 1e-7


def test_saddle_conditions_hold_on_random_instances():
    for seed in SEEDS:
        result = programs.saddle(random_instance(seed)[0])
        assert result["constraint"] < 1e-4
        assert max(result["above"], result["slackness"], result["normalization"]) < 1e-7


def test_a_place_below_the_common_value_carries_no_weight():
    instances = [random_instance(seed)[0] for seed in SEEDS] + [model.ring(3, 1.0, 1.0, [1.0, 0.0, 0.0])]
    found_below = 0
    for dens in instances:
        result = programs.saddle(dens)
        below = result["densities"] < result["level"] - 1e-5
        found_below += int(below.sum())
        assert np.all(result["weights"][below] < 1e-6)
    assert found_below > 0


def test_every_weighting_is_some_divisions():
    """Shifting each place's zero against the values a ground state of the target's weighted sum gives them, by
    constants summing to zero, makes the target weighting optimal: its floor reaches the common value."""
    dens = model.ring(5, 1.0, 1.0)
    for target in (np.array([0.5, 0.0, 0.0, 0.25, 0.25]), np.array([0.1, 0.3, 0.2, 0.15, 0.25])):
        ground = np.linalg.eigh(programs.weighted(dens, target))[1][:, 0]
        values = np.array([np.real(ground.conj() @ d @ ground) for d in dens])
        shifts = values.mean() - values
        shifted = [d + c * np.eye(d.shape[0]) for d, c in zip(dens, shifts)]
        assert abs(shifts.sum()) < 1e-12 and np.allclose(sum(shifted), sum(dens))
        assert abs(programs.certificate(shifted, target) - programs.saddle(shifted)["level"]) < 1e-6


def test_the_models_divide_their_energy_canonically():
    """Every model's densities are the canonical division of its Hamiltonian, deepened cells included, and the
    shifts that can make any weighting optimal break it."""
    specs = [{"family": "ising", "lattice": "ring", "size": 4, "g": 1.0, "J": 1.0, "depth": 0.7},
             {"family": "ising", "lattice": "chain", "size": 4, "g": 0.5, "J": 1.0},
             {"family": "xxz", "lattice": "star", "size": 4, "J": 1.0, "delta": 0.5, "h": 0.2},
             {"family": "random", "lattice": "random", "size": 4, "seed": 1}]
    for spec in specs:
        dens = model.build(spec)
        canonical = model.canonical_division(sum(dens), spec["size"])
        assert all(np.allclose(a, b) for a, b in zip(dens, canonical)), spec
    dens = model.ring(4, 1.0, 1.0)
    shifted = [d + c * np.eye(d.shape[0]) for d, c in zip(dens, (0.3, -0.1, -0.1, -0.1))]
    assert not all(np.allclose(a, b) for a, b in zip(shifted, model.canonical_division(sum(shifted), 4)))


def test_the_canonical_division_meets_the_four_asks():
    """The canonical division is linear and unchanged by a change of basis at each place; the lowest-zero rule is
    not linear."""
    rng = np.random.default_rng(11)
    cells = 3

    def hermitian():
        raw = rng.normal(size=(8, 8)) + 1j * rng.normal(size=(8, 8))
        return raw + raw.conj().T

    def unitary():
        q, r = np.linalg.qr(rng.normal(size=(2, 2)) + 1j * rng.normal(size=(2, 2)))
        return q * (np.diag(r) / np.abs(np.diag(r)))

    def lowest_zero(h):
        parts = model.components(h, cells)
        return [sum((part - np.linalg.eigvalsh(part)[0] * np.eye(8)) / len(s) for s, part in parts.items() if x in s)
                + parts[frozenset()] / cells for x in range(cells)]

    h, g = hermitian(), hermitian()
    divide = lambda m: model.canonical_division(m, cells)
    assert all(np.allclose(a, b + c) for a, b, c in zip(divide(h + g), divide(h), divide(g)))
    local = np.kron(np.kron(unitary(), unitary()), unitary())
    turned = divide(local @ h @ local.conj().T)
    assert all(np.allclose(a, local @ b @ local.conj().T) for a, b in zip(turned, divide(h)))
    assert not all(np.allclose(a, b + c) for a, b, c in zip(lowest_zero(h + g), lowest_zero(h), lowest_zero(g)))


def test_any_state_and_weighting_bracket_the_common_value():
    """By sampling: random states have peaks at or above the computed common value and random weightings floors at
    or below it, and the computed pair's peak exceeds its floor by no more than the programs' gap."""
    rng = np.random.default_rng(7)
    for seed in SEEDS:
        dens = random_instance(seed)[0]
        result = programs.saddle(dens)
        level, dim = result["level"], dens[0].shape[0]
        for _ in range(20):
            weights = rng.dirichlet(np.ones(len(dens)))
            assert np.linalg.eigvalsh(programs.weighted(dens, weights))[0] <= level + 1e-7
            vectors = rng.normal(size=(dim, 3)) + 1j * rng.normal(size=(dim, 3))
            rho = vectors @ vectors.conj().T / np.real(np.trace(vectors @ vectors.conj().T))
            assert max(float(np.real(np.trace(rho @ d))) for d in dens) >= level - 1e-7
        floor = np.linalg.eigvalsh(programs.weighted(dens, result["weights"]))[0]
        assert max(result["densities"]) - floor <= result["gap"] + 1e-6


def test_the_conditions_certify_only_at_the_lowest_eigenvalue():
    """By example: with D_1 = diag(-1, 0, 5) and D_2 = diag(5, 0, -1) the common value is 0, and the middle basis state
    with the weighting (1, 0) meets the four conditions at 0, an eigenvalue of the weighted sum but not its lowest,
    while that weighting's floor is -1."""
    dens = [np.diag([-1.0, 0.0, 5.0]), np.diag([5.0, 0.0, -1.0])]
    rho = np.diag([0.0, 1.0, 0.0])
    found = programs.residuals(dens, rho, [1.0, 0.0], 0.0)
    assert max(found.values()) < 1e-12
    assert abs(programs.saddle(dens)["level"]) < 1e-7
    assert np.linalg.eigvalsh(dens[0])[0] == -1.0
