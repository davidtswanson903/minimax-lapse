"""The setting: places, links between them, a Hermitian piece on each place and
on each link, and each place's local density, its own piece plus half of each
link piece it touches. The densities sum to the Hamiltonian."""
from __future__ import annotations

from functools import reduce

import numpy as np

from repokit.equations import equation

PAULI_I = np.eye(2)
PAULI_X = np.array([[0.0, 1.0], [1.0, 0.0]])
PAULI_Z = np.diag([1.0, -1.0])


@equation("eq:sites")
def on_site(matrix: np.ndarray, site: int, cells: int) -> np.ndarray:
    """O_x = I^(x) (x) O (x) I^(n-1-x), qubit 0 the leftmost factor (eq:sites)."""
    return reduce(np.kron, [matrix if k == site else PAULI_I for k in range(cells)])


@equation("eq:density")
def densities(place_pieces, link_pieces, links) -> list[np.ndarray]:
    """D_x = h_x + 1/2 sum of the pieces of the links at x, so that sum_x D_x = H (eq:density)."""
    out = [np.array(piece, dtype=float) for piece in place_pieces]
    for piece, (a, b) in zip(link_pieces, links):
        out[a] = out[a] + 0.5 * piece
        out[b] = out[b] + 0.5 * piece
    return out


def ring_links(cells: int) -> list[tuple[int, int]]:
    if cells < 3:
        raise ValueError("a ring needs at least three cells")
    return [(x, (x + 1) % cells) for x in range(cells)]


@equation("eq:ising", "eq:deepen")
def ring(cells: int, field: float, coupling: float, depths=None) -> list[np.ndarray]:
    """The ring of the introduction: h_x = -g X_x - mu_x Z_x and h_{x,x+1} = -J Z_x Z_{x+1},
    the Ising pieces (eq:ising) with cell x deepened by mu_x (eq:deepen)."""
    depths = list(depths) if depths is not None else [0.0] * cells
    links = ring_links(cells)
    places = [-field * on_site(PAULI_X, x, cells) - depths[x] * on_site(PAULI_Z, x, cells) for x in range(cells)]
    bonds = [-coupling * on_site(PAULI_Z, a, cells) @ on_site(PAULI_Z, b, cells) for a, b in links]
    return densities(places, bonds, links)


PAULI_Y = np.array([[0.0, -1.0j], [1.0j, 0.0]])
RANDOM_GRAPH_SEED_OFFSET = 1000
RANDOM_GRAPH_EXTRA_LINK = 0.3
FAMILY_EQUATIONS = {"ising": "eq:ising", "xxz": "eq:xxz", "random": "eq:random"}


@equation("eq:sites")
def two_site(a: np.ndarray, b: np.ndarray, x: int, y: int, cells: int) -> np.ndarray:
    """O_x O'_y, the product of two site operators (eq:sites); the products used here are real."""
    product = on_site(a, x, cells) @ on_site(b, y, cells)
    return np.real(product)


@equation("eq:lattices")
def lattice(kind: str, cells: int, seed: int = 0) -> list[tuple[int, int]]:
    """The links of each lattice (eq:lattices): ring {x, x+1 mod n}; open chain {x, x+1}; two-leg
    ladder legs {rm+i, rm+i+1} and rungs {i, i+m}, m = n/2; complete graph every pair; star
    {0, x}; random graph: each x >= 1 linked to an earlier cell drawn uniformly, then each
    remaining pair with probability RANDOM_GRAPH_EXTRA_LINK, seeded RANDOM_GRAPH_SEED_OFFSET + seed."""
    if kind == "ring":
        return ring_links(cells)
    if kind == "chain":
        return [(x, x + 1) for x in range(cells - 1)]
    if kind == "ladder":
        if cells % 2 or cells < 4:
            raise ValueError("a ladder needs an even number of cells, at least four")
        legs = cells // 2
        rows = [(r * legs + i, r * legs + i + 1) for r in (0, 1) for i in range(legs - 1)]
        return rows + [(i, i + legs) for i in range(legs)]
    if kind == "complete":
        return [(a, b) for a in range(cells) for b in range(a + 1, cells)]
    if kind == "star":
        return [(0, x) for x in range(1, cells)]
    if kind == "random":
        rng = np.random.default_rng(RANDOM_GRAPH_SEED_OFFSET + seed)
        links = [(int(rng.integers(0, x)), x) for x in range(1, cells)]  # a random spanning tree
        extra = [(a, b) for a in range(cells) for b in range(a + 1, cells) if (a, b) not in links]
        for pair in extra:
            if rng.random() < RANDOM_GRAPH_EXTRA_LINK:
                links.append(pair)
        return links
    raise ValueError(f"no lattice {kind!r}")


@equation("eq:ising")
def ising_pieces(cells: int, links, g: float, J: float, h: float = 0.0):
    """h_x = -g X_x - b Z_x and h_xy = -J Z_x Z_y, the longitudinal field b given as h (eq:ising)."""
    places = [-g * on_site(PAULI_X, x, cells) - h * on_site(PAULI_Z, x, cells) for x in range(cells)]
    bonds = [-J * two_site(PAULI_Z, PAULI_Z, a, b, cells) for a, b in links]
    return places, bonds


@equation("eq:xxz")
def xxz_pieces(cells: int, links, J: float, delta: float, h: float):
    """h_x = -b Z_x and h_xy = J (X_x X_y + Y_x Y_y) + Delta Z_x Z_y, the field b given as h (eq:xxz)."""
    places = [-h * on_site(PAULI_Z, x, cells) for x in range(cells)]
    bonds = [J * (two_site(PAULI_X, PAULI_X, a, b, cells) + two_site(PAULI_Y, PAULI_Y, a, b, cells))
             + delta * two_site(PAULI_Z, PAULI_Z, a, b, cells) for a, b in links]
    return places, bonds


RANDOM_PAIRS = ("XX", "XZ", "ZX", "ZZ", "YY")


@equation("eq:random")
def random_pieces(cells: int, links, seed: int):
    """h_x = a_x X_x + b_x Z_x and h_xy = sum over OO' in XX, XZ, ZX, ZZ, YY of c_xy^OO' O_x O'_y,
    every coefficient a standard normal draw from NumPy's default generator seeded with the
    seed, drawn in that order: a_x then b_x cell by cell, then link by link (eq:random)."""
    rng = np.random.default_rng(seed)
    named = {"X": PAULI_X, "Y": PAULI_Y, "Z": PAULI_Z}
    places = [rng.normal() * on_site(PAULI_X, x, cells) + rng.normal() * on_site(PAULI_Z, x, cells)
              for x in range(cells)]
    bonds = [sum(rng.normal() * two_site(named[p[0]], named[p[1]], a, b, cells) for p in RANDOM_PAIRS)
             for a, b in links]
    return places, bonds


@equation("eq:deepen")
def deepen(places: list, cells: int, at: int, depth: float) -> list:
    """h_y -> h_y - mu Z_y for the deepened place y (eq:deepen)."""
    out = list(places)
    out[at] = out[at] - depth * on_site(PAULI_Z, at, cells)
    return out


def build(spec: dict) -> list[np.ndarray]:
    """The local densities of one catalogue model: its family's pieces on its lattice, one cell
    deepened, split into densities."""
    cells, family = spec["size"], spec["family"]
    links = lattice(spec["lattice"], cells, spec.get("seed", 0))
    if family == "ising":
        places, bonds = ising_pieces(cells, links, spec["g"], spec["J"], spec.get("h", 0.0))
    elif family == "xxz":
        places, bonds = xxz_pieces(cells, links, spec["J"], spec["delta"], spec["h"])
    elif family == "random":
        places, bonds = random_pieces(cells, links, spec["seed"])
    else:
        raise ValueError(f"no family {family!r}")
    places = deepen(places, cells, deep_index(spec), spec.get("depth", 0.0))
    return densities(places, bonds, links)


def deep_index(spec: dict) -> int:
    """The cell a depth deepens: the first, unless the spec names the middle or the last."""
    return {"middle": spec["size"] // 2, "last": spec["size"] - 1}.get(spec.get("deep_at"), 0)


def translation_invariant(spec: dict) -> bool:
    """A rotation carries every density to the next: uniform pieces on a ring or a complete graph."""
    return spec["lattice"] in ("ring", "complete") and spec["family"] in ("ising", "xxz") and not spec.get("depth")


def _outside(h: np.ndarray, keep, cells: int) -> np.ndarray:
    """E_T(H): H averaged over every place outside T and replaced there by the identity, one place at a time by
    the twirl (H + X H X + Y H Y + Z H Z)/4 = tr_x(H)/2 (x) I_x."""
    for place in range(cells):
        if place not in keep:
            twirl = [on_site(p, place, cells) for p in (PAULI_X, PAULI_Y, PAULI_Z)]
            h = (h + sum(o @ h @ o for o in twirl)) / 4
    return h


@equation("eq:canonical")
def components(h: np.ndarray, cells: int) -> dict:
    """H_S = sum over T in S of (-1)^{|S - T|} E_T(H): the part of H acting on the places S and averaging to zero
    over each of them, for every set S of places (eq:canonical)."""
    from itertools import combinations
    sets = [frozenset(c) for size in range(cells + 1) for c in combinations(range(cells), size)]
    averaged = {s: _outside(h, s, cells) for s in sets}
    return {s: sum((-1) ** (len(s) - len(t)) * averaged[t] for t in sets if t <= s) for s in sets}


def canonical_division(h: np.ndarray, cells: int) -> list[np.ndarray]:
    """The canonical densities: each part H_S shared equally among the places of S, the constant part among all."""
    parts = components(h, cells)
    return [sum(part / len(s) for s, part in parts.items() if place in s) + parts[frozenset()] / cells
            for place in range(cells)]
