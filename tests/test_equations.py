"""Each equation function, checked by another route than the code it checks. The first line
of each docstring says which route."""
import itertools
from functools import reduce

import numpy as np

from repokit import equations, sources
from repokit.equations import checks
from saddle import entropy, laws, model, probes, programs


def pauli_string(ops: dict, n: int) -> np.ndarray:
    """The operator of a Pauli string, built by acting on basis states bit by bit, qubit 0 the
    most significant bit: Z multiplies by (-1)^bit, X flips the bit, Y flips it and multiplies by
    i (-1)^bit."""
    dim = 2 ** n
    out = np.zeros((dim, dim), dtype=complex)
    for i in range(dim):
        j, amplitude = i, 1 + 0j
        for x, op in ops.items():
            bit = (i >> (n - 1 - x)) & 1
            if op == "Z":
                amplitude *= -1 if bit else 1
            elif op == "X":
                j ^= 1 << (n - 1 - x)
            elif op == "Y":
                j ^= 1 << (n - 1 - x)
                amplitude *= -1j if bit else 1j
        out[j, i] += amplitude
    return out


def real(ops: dict, n: int) -> np.ndarray:
    matrix = pauli_string(ops, n)
    assert np.allclose(matrix.imag, 0)
    return matrix.real


@checks("eq:sites")
def test_site_operators_act_bit_by_bit():
    """By acting on basis states bit by bit: every one- and two-site operator is the stated tensor product."""
    n = 3
    named = {"X": model.PAULI_X, "Y": model.PAULI_Y, "Z": model.PAULI_Z}
    for x in range(n):
        for name in "XYZ":
            assert np.allclose(model.on_site(named[name], x, n), pauli_string({x: name}, n))
    for x, y in itertools.permutations(range(n), 2):
        for a, b in itertools.product("XYZ", repeat=2):
            expected = pauli_string({x: a, y: b}, n)
            if np.allclose(expected.imag, 0):
                assert np.allclose(model.two_site(named[a], named[b], x, y, n), expected.real)


@checks("eq:density")
def test_densities_split_the_hamiltonian():
    """By summing every piece once and splitting each link by hand: the densities add to the Hamiltonian, and each holds its cell's piece and half of each of its links."""
    n, rng = 4, np.random.default_rng(7)
    links = [(0, 1), (1, 2), (2, 3), (3, 0), (0, 2)]
    places = [rng.normal() * real({x: "X"}, n) for x in range(n)]
    bonds = [rng.normal() * real({a: "Z", b: "Z"}, n) for a, b in links]
    dens = model.densities(places, bonds, links)
    assert np.allclose(sum(dens), sum(places) + sum(bonds))
    for x in range(n):
        assert np.allclose(dens[x], places[x] + 0.5 * sum(b for b, link in zip(bonds, links) if x in link))


@checks("eq:ising")
def test_ising_pieces_are_the_stated_ones():
    """By building the stated Ising pieces from bitwise Pauli strings: they match, cell by cell and link by link."""
    n, links = 4, [(0, 1), (1, 2), (2, 3), (3, 0)]
    places, bonds = model.ising_pieces(n, links, g=0.7, J=-1.3, h=0.4)
    for x in range(n):
        assert np.allclose(places[x], -0.7 * real({x: "X"}, n) - 0.4 * real({x: "Z"}, n))
    for bond, (a, b) in zip(bonds, links):
        assert np.allclose(bond, 1.3 * real({a: "Z", b: "Z"}, n))


@checks("eq:xxz")
def test_xxz_pieces_are_the_stated_ones():
    """By building the stated XXZ pieces from bitwise Pauli strings: they match, cell by cell and link by link."""
    n, links = 4, [(0, 1), (1, 2), (2, 3)]
    places, bonds = model.xxz_pieces(n, links, J=0.9, delta=1.5, h=0.2)
    for x in range(n):
        assert np.allclose(places[x], -0.2 * real({x: "Z"}, n))
    for bond, (a, b) in zip(bonds, links):
        expected = 0.9 * (real({a: "X", b: "X"}, n) + real({a: "Y", b: "Y"}, n)) + 1.5 * real({a: "Z", b: "Z"}, n)
        assert np.allclose(bond, expected)


@checks("eq:random")
def test_random_pieces_follow_the_stated_draws():
    """By drawing the stated coefficients in the stated order and building the pieces from bitwise Pauli strings: they match."""
    n, links, seed = 3, [(0, 1), (1, 2)], 5
    places, bonds = model.random_pieces(n, links, seed)
    rng = np.random.default_rng(seed)
    for x in range(n):
        a, b = rng.normal(), rng.normal()
        assert np.allclose(places[x], a * real({x: "X"}, n) + b * real({x: "Z"}, n))
    for bond, (x, y) in zip(bonds, links):
        expected = sum(rng.normal() * real({x: p[0], y: p[1]}, n) for p in ("XX", "XZ", "ZX", "ZZ", "YY"))
        assert np.allclose(bond, expected)


@checks("eq:lattices")
def test_lattices_have_the_stated_links():
    """By counting links and walking each graph: every lattice has the stated number of links, and the random graphs are connected."""
    for n in (3, 4, 5, 6):
        assert len(set(map(frozenset, model.lattice("ring", n)))) == n
        assert len(model.lattice("chain", n)) == n - 1
        assert len(model.lattice("complete", n)) == n * (n - 1) // 2
        assert sorted(model.lattice("star", n)) == [(0, x) for x in range(1, n)]
        for seed in range(5):
            links = model.lattice("random", n, seed)
            assert len(links) >= n - 1 and len(set(map(frozenset, links))) == len(links)
            reached, frontier = {0}, [0]
            while frontier:
                here = frontier.pop()
                for a, b in links:
                    for other in ((b,) if a == here else (a,) if b == here else ()):
                        if other not in reached:
                            reached.add(other)
                            frontier.append(other)
            assert reached == set(range(n))
    for n in (4, 6):
        assert len(model.lattice("ladder", n)) == 3 * n // 2 - 2


@checks("eq:deepen")
def test_deepening_moves_only_the_deep_cell():
    """By differencing the densities with and without depth: only the deepened cell changes, by minus the depth times its spin."""
    for spec, y in (({"family": "ising", "lattice": "chain", "size": 4, "g": 1.0, "J": 1.0, "deep_at": "middle"}, 2),
                    ({"family": "xxz", "lattice": "star", "size": 4, "J": 1.0, "delta": 0.5, "h": 0.2, "deep_at": "last"}, 3)):
        plain, deep = model.build(spec), model.build({**spec, "depth": 0.7})
        for x in range(4):
            assert np.allclose(deep[x] - plain[x], -0.7 * real({y: "Z"}, 4) if x == y else 0.0)


@checks("eq:saddle")
def test_both_sides_meet_a_search_over_weightings():
    """By searching the simplex of weightings on a grid, refined near its best point, without either program: the best weighted ground energy found matches the right side, and the left side is not below it."""
    dens = model.ring(3, 1.0, 1.0, [0.2, 0.0, 0.0])
    least = lambda n: np.linalg.eigvalsh(sum(w * d for w, d in zip(n, dens)))[0]
    coarse = [(least((a, b, 1 - a - b)), a, b) for a, b in
              ((i / 100, j / 100) for i in range(101) for j in range(101 - i))]
    best, a0, b0 = max(coarse)
    for i, j in itertools.product(range(-20, 21), repeat=2):
        a, b = a0 + i / 2000, b0 + j / 2000
        if a >= 0 and b >= 0 and a + b <= 1:
            best = max(best, least((a, b, 1 - a - b)))
    right, left = programs.greatest_least(dens)["value"], programs.least_greatest(dens)["value"]
    assert abs(best - right) < 1e-3 and left >= best - 1e-7


@checks("eq:certificate")
def test_every_weighting_bounds_every_state():
    """By drawing random states and weightings: no state's largest density falls below any weighting's weighted ground energy."""
    rng = np.random.default_rng(3)
    dens = model.build({"family": "random", "lattice": "random", "size": 3, "seed": 1})
    for _ in range(200):
        g = rng.normal(size=(8, 8)) + 1j * rng.normal(size=(8, 8))
        rho = g @ g.conj().T
        rho /= np.trace(rho).real
        peak = max(float(np.real(np.trace(rho @ d))) for d in dens)
        assert peak >= programs.certificate(dens, rng.dirichlet(np.ones(3))) - 1e-12


@checks("eq:conditions")
def test_conditions_hold_at_a_saddle_built_by_hand():
    """By building the uniform ring's saddle by hand, its ground state with equal weights, and moving off it: every residual vanishes there, and the first does not off it."""
    dens = model.ring(3, 1.0, 1.0)
    energies, vectors = np.linalg.eigh(sum(dens))
    rho = np.outer(vectors[:, 0], vectors[:, 0])
    at = programs.residuals(dens, rho, np.full(3, 1 / 3), energies[0] / 3)
    off = programs.residuals(dens, rho, np.array([0.5, 0.25, 0.25]), energies[0] / 3)
    assert max(at.values()) < 1e-10 and off["constraint"] > 1e-3


@checks("eq:rates")
def test_rates_are_the_one_sided_differences():
    """By shifting one place's density and solving again: the forward and backward differences of the common value meet the greatest and least weight the place takes, even where those differ."""
    dens = model.ring(3, 1.0, 1.0, [0.1, 0.0, 0.0])
    split = [dens[0], dens[0], dens[1], dens[2]]  # one place split in two: its weight may sit in either copy
    low, high = programs.rate_bounds(split, 0)
    common = programs.greatest_least(split)["value"]
    def moved(step):
        out = list(split)
        out[0] = out[0] + step * np.eye(8)
        return programs.greatest_least(out)["value"]
    forward, backward = (moved(1e-3) - common) / 1e-3, (common - moved(-1e-3)) / 1e-3
    assert high - low > 0.1 and abs(forward - high) < 5e-3 and abs(backward - low) < 5e-3


@checks("eq:evolution")
def test_weighted_hamiltonian_two_ways_and_its_stationary_state():
    """By assembling the weighted Hamiltonian from cells and links directly, and by solving for a saddle: the two forms agree, and the saddle state commutes with it."""
    n, links = 3, model.lattice("ring", 3)
    places, bonds = model.ising_pieces(n, links, 1.0, 1.0)
    places = model.deepen(places, n, 0, 0.2)
    dens = model.densities(places, bonds, links)
    weights = np.array([0.2, 0.5, 0.3])
    assert np.allclose(programs.weighted(dens, weights), programs.weighted_from_pieces(places, bonds, links, weights))
    held = programs.weighted(dens, programs.saddle(dens)["weights"])
    rho = programs.least_greatest(dens)["state"]
    assert np.linalg.norm(held @ rho - rho @ held) < 1e-4, "the saddle state commutes with it to the solver's accuracy"


@checks("eq:uniform")
def test_uniform_value_is_what_the_programs_find():
    """By solving both programs on uniform rings and a complete graph: their common value is the ground energy per cell."""
    for spec in ({"family": "ising", "lattice": "ring", "size": 4, "g": 1.0, "J": 1.0},
                 {"family": "xxz", "lattice": "ring", "size": 4, "J": 1.0, "delta": 0.5, "h": 0.2},
                 {"family": "ising", "lattice": "complete", "size": 4, "g": 1.0, "J": -1.0}):
        dens = model.build(spec)
        assert abs(programs.saddle(dens)["level"] - programs.uniform_value(dens)) < 1e-6


@checks("eq:dropout")
def test_drop_depth_is_where_the_common_value_stops_moving():
    """By solving the whole problem at depths around it: just past the depth the common value is the other cells' own, and just before it is higher."""
    spec = {"family": "ising", "lattice": "ring", "size": 3, "g": 1.0, "J": 1.0}
    depth = programs.drop_depth(model.build(spec), 0, model.on_site(model.PAULI_Z, 0, 3))
    others = programs.least_greatest(model.build(spec)[1:])["value"]
    after = programs.saddle(model.build({**spec, "depth": depth + 0.01}))["level"]
    before = programs.saddle(model.build({**spec, "depth": depth - 0.02}))["level"]
    assert abs(after - others) < 1e-6 and before - others > 1e-5


@checks("eq:law")
def test_susceptibilities_match_finite_differences():
    """By switching on a small field at each other cell and differencing the spins: the summed responses match the susceptibility formula."""
    n, weights, step = 4, np.array([0.0, 1 / 3, 1 / 3, 1 / 3]), 1e-5
    links = [(x, (x + 1) % n) for x in range(n)]
    held = sum(0.5 * (weights[a] + weights[b]) * (real({a: "X", b: "X"}, n) + real({a: "Y", b: "Y"}, n) + 0.5 * real({a: "Z", b: "Z"}, n))
               for a, b in links)
    spins = [real({x: "Z"}, n) for x in range(n)]
    sector = np.where(np.isclose(np.diag(sum(spins)), 0))[0]
    def others_total(field_at, field):
        matrix = (held + field * spins[field_at])[np.ix_(sector, sector)]
        ground = np.linalg.eigh(matrix)[1][:, 0]
        return sum(ground @ s[np.ix_(sector, sector)] @ ground for s in spins[1:])
    differenced = [-(others_total(y, step) - others_total(y, -step)) / (2 * step) for y in range(1, n)]
    assert np.allclose(laws.response(n, weights), differenced, rtol=1e-4, atol=1e-8)


def test_every_model_in_the_catalogue_is_stated():
    """Every family and lattice the catalogue uses has its equation in the paper and its function in the code."""
    labels = equations.paper_labels(sources.paper())
    computed = {e.label for e in equations.functions()}
    used = {row["spec"]["family"] for row in sources.results()["catalogue"]["models"]}
    for family in used:
        label = model.FAMILY_EQUATIONS[family]
        assert label in labels and label in computed, family
    assert "eq:lattices" in labels and "eq:lattices" in computed


def _probed(dens: list, omega: float) -> list:
    """The densities with an idle probe qubit attached to every place, built directly: D'_x = D_x (x) I + K on probe x."""
    n, size = len(dens), dens[0].shape[0]
    piece = 0.5 * omega * np.diag([1.0, -1.0])
    out = []
    for x, d in enumerate(dens):
        on_probe = reduce(np.kron, [piece if k == x else np.eye(2) for k in range(n)])
        out.append(np.kron(d, np.eye(2 ** n)) + np.kron(np.eye(size), on_probe))
    return out


def _probe_operator(op: np.ndarray, x: int, n: int, size: int) -> np.ndarray:
    return np.kron(np.eye(size), reduce(np.kron, [op if k == x else np.eye(2) for k in range(n)]))


def _joint(depth: float = 0.1, omega: float = 1.0):
    dens = model.ring(3, 1.0, 1.0, [depth, 0.0, 0.0])
    found = programs.saddle(dens)
    weights = np.clip(found["weights"], 0.0, None)
    return dens, weights, found, _probed(dens, omega)


@checks("eq:probes")
def test_probes_leave_the_weights_where_they_were():
    """By solving the saddle of the system with its probes attached, in one space: the weights do not move, and the common value moves by the probes' ground energy."""
    dens, weights, found, probed = _joint()
    joint = programs.saddle(probed)
    piece = 0.5 * np.diag([1.0, -1.0])
    assert np.allclose(joint["weights"], weights, atol=1e-4)
    assert abs(joint["level"] - probes.probed_level(found["level"], piece)) < 1e-6


@checks("eq:clockrate")
def test_clocks_run_at_their_weights_in_the_joint_evolution():
    """By evolving the system and its probes together under the probed generator, without factoring it: each probe clock's phase advances at its place's weight."""
    dens, weights, _, probed = _joint()
    held = programs.weighted(probed, weights)
    energies, vectors = np.linalg.eigh(held)
    ground = np.linalg.eigh(programs.weighted(dens, weights))[1][:, 0]
    plus = np.array([1.0, 1.0]) / np.sqrt(2)
    state = np.kron(ground, reduce(np.kron, [plus] * 3))
    times = np.linspace(0.0, 10.0, 40)
    X, Y = np.array([[0, 1], [1, 0]], dtype=complex), np.array([[0, -1j], [1j, 0]])
    for x in range(3):
        sx, sy = _probe_operator(X, x, 3, 8), _probe_operator(Y, x, 3, 8)
        phases = []
        for t in times:
            evolved = vectors @ (np.exp(-1j * energies * t) * (vectors.conj().T @ state))
            phases.append(np.arctan2(np.real(np.vdot(evolved, sy @ evolved)), np.real(np.vdot(evolved, sx @ evolved))))
        rate = np.polyfit(times, np.unwrap(phases), 1)[0]
        assert abs(rate - probes.clock_rate(weights[x], 1.0)) < 1e-8


@checks("eq:tolman")
def test_thermometers_read_tolman_in_the_joint_equilibrium():
    """By forming the thermal state of the probed generator in one space and reading each probe's populations: temperature times weight is the global temperature at every place."""
    beta = 2.0
    dens, weights, _, probed = _joint()
    energies, vectors = np.linalg.eigh(programs.weighted(probed, weights))
    boltzmann = np.exp(-beta * (energies - energies.min()))
    thermal = (vectors * (boltzmann / boltzmann.sum())) @ vectors.conj().T
    up, down = np.diag([1.0, 0.0]), np.diag([0.0, 1.0])
    for x in range(3):
        excited = np.real(np.trace(thermal @ _probe_operator(up, x, 3, 8)))
        ground = np.real(np.trace(thermal @ _probe_operator(down, x, 3, 8)))
        reading = 1.0 / np.log(ground / excited)
        assert abs(reading - probes.local_temperature(weights[x], beta)) < 1e-8 * reading


@checks("eq:redshift")
def test_an_excitation_adds_its_weight_times_its_energy():
    """By exciting one probe in the joint ground state and differencing the probed generator's energy: the whole gains the place's weight times the probe's energy."""
    dens, weights, _, probed = _joint()
    held = programs.weighted(probed, weights)
    ground = np.linalg.eigh(programs.weighted(dens, weights))[1][:, 0]
    down, up = np.array([0.0, 1.0]), np.array([1.0, 0.0])
    resting = np.kron(ground, reduce(np.kron, [down] * 3))
    for x in range(3):
        excited = np.kron(ground, reduce(np.kron, [up if k == x else down for k in range(3)]))
        gain = np.vdot(excited, held @ excited).real - np.vdot(resting, held @ resting).real
        assert abs(gain - probes.global_energy(weights[x], 1.0)) < 1e-10


@checks("eq:phase")
def test_the_saddle_state_turns_at_the_common_value():
    """By evolving the state the state program finds, not an eigenvector, under the weighted generator: its overlap with itself turns at the common value."""
    dens = model.ring(3, 1.0, 1.0, [0.1, 0.0, 0.0])
    found = programs.saddle(dens)
    rho = programs.least_greatest(dens)["state"]
    energies, vectors = np.linalg.eigh(programs.weighted(dens, found["weights"]))
    times = np.linspace(0.0, 10.0, 40)
    overlap = [np.trace(rho @ (vectors * np.exp(-1j * energies * t)) @ vectors.conj().T) for t in times]
    rate = -np.polyfit(times, np.unwrap(np.angle(overlap)), 1)[0]
    assert abs(rate - probes.phase_rate(found["level"])) < 1e-6


@checks("eq:curvature")
def test_the_floor_curves_as_its_second_differences_say():
    """By second differences of the lowest eigenvalue along random pairs of directions, not perturbation theory: at a
    weighting with a single lowest state the floor's mixed second difference meets the Hessian's."""
    dens = model.build({"family": "ising", "lattice": "ring", "size": 5, "g": 1.0, "J": 1.0, "depth": 0.3})
    rng = np.random.default_rng(1)
    weights = rng.dirichlet(np.ones(5))
    hessian, gap = programs.floor_hessian(dens, weights)
    assert gap > 1e-3
    floor = lambda w: np.linalg.eigvalsh(programs.weighted(dens, w))[0]
    h = 1e-4
    for _ in range(3):
        v, u = rng.normal(size=5), rng.normal(size=5)
        mixed = (floor(weights + h * v + h * u) - floor(weights + h * v - h * u) - floor(weights - h * v + h * u)
                 + floor(weights - h * v - h * u)) / (4 * h * h)
        assert abs(mixed - v @ hessian @ u) < 1e-3 * max(1.0, abs(mixed))  # second differences at step 1e-4


@checks("eq:coldest")
def test_coldest_equilibria_approach_the_saddle_weights():
    """By checking each maximum-entropy state's own optimality conditions, and the saddle program separately: every state holds every local energy within its bound, carries multipliers only where an energy meets it, and as the bound tightens its normalized multipliers close on the program's weights while their sum grows."""
    dens = model.ring(3, 1.0, 1.0, [0.2, 0.0, 0.0])
    found = programs.saddle(dens)
    distances, totals, multipliers = [], [], None
    for gap in (0.3, 0.03, 0.003):
        bound = found["level"] + gap
        multipliers = entropy.coldest_multipliers(dens, bound, start=multipliers)
        state = entropy.thermal_state(dens, multipliers)
        energies = np.array([np.real(np.trace(state @ d)) for d in dens])
        assert np.all(energies <= bound + 1e-7)
        assert np.all(np.abs(multipliers * (energies - bound)) < 1e-5)
        distances.append(np.abs(entropy.coldest_shape(multipliers) - found["weights"]).max())
        totals.append(multipliers.sum())
    assert distances[0] > distances[1] > distances[2] and distances[2] < distances[0] / 3, "the proportions close on the weights"
    assert totals[0] < totals[1] < totals[2], "while the multipliers grow without bound"


@checks("eq:entanglement")
def test_fitting_recovers_the_weights_a_state_was_built_from():
    """By building a state as a weighted thermal state of the terms and fitting its measured terms: the fit returns the weights it was built from."""
    n = 3
    terms = [-model.on_site(model.PAULI_X, x, n) for x in range(n)] + [-model.two_site(model.PAULI_Z, model.PAULI_Z, x, x + 1, n) for x in range(n - 1)]
    built = np.array([0.4, 1.1, 1.9, 0.7, 1.6])
    state = entropy.thermal_state(terms, built)
    fitted = entropy.entanglement_weights(terms, [np.real(np.trace(state @ t)) for t in terms], start=np.ones(len(terms)))
    assert np.allclose(fitted, built, atol=1e-5)


@checks("eq:gain")
def test_gain_measures_what_reweighting_adds():
    """By searching the weightings on a grid and diagonalizing the Hamiltonian directly: the best grid weighting's lowest eigenvalue, less the ground energy per place, matches the gain, which vanishes on a uniform ring and not on an open chain."""
    def grid_gain(dens):
        best = max(np.linalg.eigvalsh(a * dens[0] + b * dens[1] + (1 - a - b) * dens[2])[0]
                   for a in np.linspace(0, 1, 81) for b in np.linspace(0, 1 - a, max(2, int(round((1 - a) * 80)) + 1)))
        return best - np.linalg.eigvalsh(dens[0] + dens[1] + dens[2])[0] / 3
    chain = model.build({"family": "ising", "lattice": "chain", "size": 3, "g": 1.0, "J": 1.0})
    ring = model.ring(3, 1.0, 1.0, [0.0, 0.0, 0.0])
    assert abs(programs.reweighting_gain(chain) - grid_gain(chain)) < 3e-3 and programs.reweighting_gain(chain) > 1e-2
    assert abs(programs.reweighting_gain(ring)) < 1e-6 and abs(grid_gain(ring)) < 3e-3


@checks("eq:conditioned")
def test_the_conditioned_rest_evolves_under_the_weighted_generator():
    """By the matrix exponential in the rest's own basis: a whole stationary under a generator plus a clock's energy, built
    from the generator's eigenvectors and clock levels that cancel them, conditioned on the clock's reading t, equals the
    rest at reading zero evolved by exp(-i A t), up to the phase exp(i lambda t)."""
    from scipy.linalg import expm

    from saddle.experiments.timeless import conditioned
    rng = np.random.default_rng(7)
    raw = rng.normal(size=(4, 4)) + 1j * rng.normal(size=(4, 4))
    generator = (raw + raw.conj().T) / 2
    energies, vectors = np.linalg.eigh(generator)
    level = 0.3                                        # lambda, the whole's eigenvalue
    clock_energy = level - energies                    # each clock level cancels one of the rest's
    amplitudes = rng.normal(size=4) + 1j * rng.normal(size=4)
    whole = np.diag(amplitudes)                        # rest level a paired with clock level a
    at = lambda t: vectors @ conditioned(whole, clock_energy, t)   # back to the rest's own basis
    for t in (0.0, 0.7, 2.3):
        expected = np.exp(1j * level * t) * expm(-1j * generator * t) @ at(0.0)
        assert np.allclose(at(t), expected, atol=1e-12)


@checks("eq:casimir")
def test_the_size_set_part_is_the_conformal_term():
    """By the closed form of the critical ring's ground energy, -2/sin(pi/2n): its per-cell excess over -4/pi is
    the conformal term at c = 1/2, v = 2, up to the next order in 1/n^2."""
    from saddle.experiments import critical

    for cells in (12, 40, 80):
        excess = critical.solve(cells, 1.0, 1.0)["energy"] / cells + 4 / np.pi
        assert abs(excess + 2 / (cells * np.sin(np.pi / (2 * cells))) - 4 / np.pi) < 1e-12
        assert abs(excess - critical.size_set_part(cells, 0.5, 2.0)) < 1 / cells ** 4


@checks("eq:chord")
def test_the_chord_carries_the_entropy():
    """By the slope of the computed entropies against the chord at c = 1/2 on a ring of forty cells."""
    from saddle.experiments import critical

    cells = 40
    correlations = critical.solve(cells, 1.0, 1.0)["correlations"]
    lengths = range(4, cells // 2 + 1)
    slope = np.polyfit([critical.chord(cells, n) for n in lengths],
                       [critical.interval_entropy(correlations, n) for n in lengths], 1)[0]
    assert abs(slope - 0.5) < 2e-3


@checks("eq:canonical")
def test_the_components_are_the_pauli_strings_on_each_set():
    """By the Pauli expansion: H_S is the sum of H's Pauli strings whose support is exactly S."""
    from saddle import model

    cells = 3
    rng = np.random.default_rng(5)
    raw = rng.normal(size=(8, 8)) + 1j * rng.normal(size=(8, 8))
    h = raw + raw.conj().T
    paulis = {"I": np.eye(2), "X": model.PAULI_X, "Y": model.PAULI_Y, "Z": model.PAULI_Z}
    expected = {}
    for word in itertools.product("IXYZ", repeat=cells):
        string = reduce(np.kron, [paulis[w] for w in word])
        support = frozenset(i for i, w in enumerate(word) if w != "I")
        expected[support] = expected.get(support, 0) + np.trace(string @ h) / 2 ** cells * string
    parts = model.components(h, cells)
    assert set(parts) == set(expected) and all(np.allclose(parts[s], expected[s]) for s in parts)
@checks("eq:moment")
def test_the_first_moment_is_the_long_wave_slope():
    """By the slope of a long wave's energy: for a translation-invariant chain with terms of ranges one and three, the
    first moment equals the ratio of its dispersion's slope at q = 0 to the chain's own."""
    from saddle.experiments import patch

    n, near, far = 40, -2 * 0.7, -2 * 0.1
    kappa = np.zeros((n, n))
    for a in range(n - 1):
        kappa[a, a + 1] = near
    for a in range(n - 3):
        kappa[a, a + 3] = far
    kappa = kappa - kappa.T
    q = 1e-4
    slope = (near * np.sin(q) + far * np.sin(3 * q)) / (-2 * np.sin(q))
    assert abs(patch.moment_weights(kappa, 1.0)[n // 2] - slope) < 1e-6


@checks("eq:patch")
def test_the_profile_has_slope_two_pi_over_v_at_its_edge():
    """By its derivative at the edge and its value at the centre: 2 pi / v and (2 pi / v) l / 4."""
    from saddle.experiments import patch

    length, wave_speed, h = 40, 2.0, 1e-6
    edge = (patch.profile(np.array([h]), length, wave_speed)[0] - 0.0) / h
    assert abs(edge - 2 * np.pi / wave_speed) < 1e-4
    assert abs(patch.profile(np.array([length / 2]), length, wave_speed)[0] - np.pi / wave_speed * length / 2) < 1e-12


@checks("eq:desitter")
def test_the_desitter_temperature_by_the_cosmological_constant():
    """By Lambda: Lambda = 3 Omega H0^2 / c^2, then T = hbar c sqrt(Lambda / 3) / (2 pi k_B)."""
    from saddle.experiments import patch

    hubble, fraction = 67.36, 0.6847
    rate = hubble * 1e3 / patch.MEGAPARSEC
    cosmological = 3 * fraction * rate ** 2 / patch.LIGHT ** 2
    temperature = patch.HBAR * patch.LIGHT * np.sqrt(cosmological / 3) / (2 * np.pi * patch.BOLTZMANN)
    assert abs(patch.desitter_temperature(hubble, fraction) / temperature - 1) < 1e-12
@checks("eq:chsh")
def test_the_chsh_value_is_the_bell_operators():
    """By the Bell operator's expectation: frequencies from the readings' projectors on its top state give S equal
    to the top eigenvalue."""
    from saddle.experiments import records

    found = records.bounds()
    project = lambda reading, outcome: (np.eye(2) + (-1) ** outcome * reading) / 2
    p = np.zeros((2, 2, 2, 2))
    for x, y, a, b in itertools.product((0, 1), repeat=4):
        p[x, y, a, b] = np.real(found["state"].conj() @ np.kron(project(records.ALICE[x], a), project(records.BOB[y], b)) @ found["state"])
    assert abs(records.chsh(p) - found["quantum"]) < 1e-12


@checks("eq:history")
def test_the_history_hamiltonian_is_a_path_laplacian_in_the_steps_frame():
    """By conjugation: with W the steps' products on each reading, W^dagger H W is the path Laplacian on the register
    tensored with the identity, plus the start's penalty at the first reading."""
    from saddle.experiments import records

    rng = np.random.default_rng(9)
    steps = [np.linalg.qr(rng.normal(size=(2, 2)) + 1j * rng.normal(size=(2, 2)))[0] for _ in range(3)]
    first = np.array([1.0, 0.0])
    count = len(steps)
    blocks, product = [np.eye(2)], np.eye(2)
    for unitary in steps:
        product = unitary @ product
        blocks.append(product)
    w = np.zeros((2 * (count + 1), 2 * (count + 1)), dtype=complex)
    for t, block in enumerate(blocks):
        w[2 * t:2 * t + 2, 2 * t:2 * t + 2] = block
    laplacian = np.zeros((count + 1, count + 1))
    for t in range(1, count + 1):
        edge = np.eye(count + 1)[t] - np.eye(count + 1)[t - 1]
        laplacian += 0.5 * np.outer(edge, edge)
    penalty = np.kron(np.diag(np.eye(count + 1)[0]), np.eye(2) - np.outer(first, first))
    assert np.allclose(w.conj().T @ records.history_hamiltonian(steps, first) @ w, np.kron(laplacian, np.eye(2)) + penalty)
@checks("eq:depth")
def test_the_closed_form_depth_is_where_the_full_saddle_drops_the_cell():
    """By the full saddle on the weighting side: on the ring of four cells the deep cell's optimal weight is positive
    just below the closed-form depth and zero just above it."""
    from scipy.optimize import minimize

    from saddle import model
    from saddle.experiments import depths

    found = depths.others_optimum(4, 1.0, 1.0)
    critical = depths.depth_from_others(found["link"], found["common"], found["flip"], 1.0)

    def deep_weight(mu):
        dens = model.ring(4, 1.0, 1.0, [mu, 0.0, 0.0, 0.0])

        def negative(w):
            levels, states = np.linalg.eigh(sum(wi * d for wi, d in zip(w, dens)))
            psi = states[:, 0]
            return -levels[0], -np.array([np.real(psi.conj() @ d @ psi) for d in dens])

        best = min((minimize(negative, start, jac=True, method="SLSQP", bounds=[(0, 1)] * 4,
                             constraints=[{"type": "eq", "fun": lambda w: w.sum() - 1, "jac": lambda w: np.ones(4)}],
                             options={"ftol": 1e-15, "maxiter": 1000})
                    for start in (np.ones(4) / 4, np.array([0.0, 1 / 3, 1 / 3, 1 / 3]))), key=lambda r: r.fun)
        return best.x[0]

    assert deep_weight(critical - 1e-3) > 1e-4 and deep_weight(critical + 1e-3) < 1e-6


@checks("eq:optical")
def test_the_optical_time_on_a_linear_lapse():
    """By the integral in closed form: on weights N = a + b x the time from x1 to x2 is ln((a + b x2)/(a + b x1)) / (v b),
    and on uniform weights it is the distance over v N."""
    from saddle.experiments import space

    x = np.arange(400, dtype=float)
    a, b = 1.0, -0.002
    weights = a + b * x
    exact = np.log((a + b * 350) / (a + b * 50)) / (space.SPEED * b)
    assert abs(space.optical_time(weights, 50.0, 350.0) - exact) < 1e-6 * exact
    assert abs(space.optical_time(np.full(400, 0.5), 10.0, 210.0) - 200.0) < 1e-9


@checks("eq:vacuum")
def test_the_vacuum_entropy_follows_the_free_fermions_on_the_uniform_chain():
    """By the free-fermion correlations of the uniform open critical chain: the formula's change from cut to cut is the
    computed entropy's, away from the ends, to within the chain's own corrections."""
    from saddle.experiments import space

    weights = np.ones(200)
    computed = space.cut_entropies(weights)
    formula = space.vacuum_entropy(weights)
    bulk = slice(40, 160)
    change = np.diff(computed[bulk]) - np.diff(formula[bulk])
    assert np.max(np.abs(change)) < 2e-3


@checks("eq:freefloor")
def test_the_free_floor_is_the_spin_rings_lowest_level():
    """By exact diagonalization of the spin ring: the floor is the lowest eigenvalue of sum_x N_x D_x, and the densities
    are its ground state's, with D_x the cell's field and half of each of its links."""
    from saddle.experiments import bump

    n = 5
    rng = np.random.default_rng(5)
    c, weights = rng.uniform(0.5, 1.5, n), rng.dirichlet(np.ones(n))
    X = [pauli_string({x: "X"}, n) for x in range(n)]
    ZZ = [pauli_string({x: "Z", (x + 1) % n: "Z"}, n) for x in range(n)]
    D = [-c[x] * X[x] - 0.5 * c[x] * ZZ[x] - 0.5 * c[x - 1] * ZZ[x - 1] for x in range(n)]
    levels, states = np.linalg.eigh(sum(w * d for w, d in zip(weights, D)))
    floor, densities = bump.free_floor(c, weights)
    assert abs(floor - levels[0]) < 1e-10
    assert np.max(np.abs(densities - [np.real(states[:, 0].conj() @ d @ states[:, 0]) for d in D])) < 1e-10


@checks("eq:interval")
def test_the_interval_entropy_follows_the_free_fermions_on_the_uniform_ring():
    """By the free-fermion correlations of a large uniform ring: for intervals short against the ring, the formula's
    change between lengths is the computed entropy's, (1/6) ln of their ratio."""
    from saddle.experiments import bump

    n = 1000
    G = bump._correlations(np.ones(n), np.ones(n) / n)
    lapse = np.ones(n)
    for l1, l2 in ((10, 20), (15, 30)):
        computed = bump._entropy(G, 100, 100 + l2) - bump._entropy(G, 100, 100 + l1)
        formula = bump.interval_entropy(lapse, 100, 100 + l2) - bump.interval_entropy(lapse, 100, 100 + l1)
        assert abs(computed - formula) < 2e-3
