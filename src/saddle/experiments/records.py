"""Chance from records inside: Bell's assumptions factored, a gluing certificate, and a protocol solved once.

Three computations. The CHSH combination of two readings on each of two parts is bounded by 2 over every
assignment of outcomes, one description holding all four readings at once, and reaches 2 sqrt 2 on a maximally
entangled pair when only readings of different parts commute. The largest share of those frequencies one classical
description can carry is a linear program, and its dual is a Bell inequality. And two friends measure the pair and
record their outcomes, each record copied outside at an overlap gamma (1 for no copy, 0 for a full one); then Alice
either reads Charlie's record or undoes his measurement and measures her qubit in another basis, and Bob likewise
with Debbie, the settings themselves held in qubits. The whole, settings, pair, records, copies and a step register,
is solved once as the ground state of its history Hamiltonian (Eq. history), with no time parameter, and every
frequency is read from it conditioned on the last step and on the records."""
from __future__ import annotations

import itertools

import cvxpy as cp
import numpy as np

from repokit.equations import equation

CONFIGS = {
    "golden": {"overlaps": [1.0, 0.0]},
    "quick": {"overlaps": [1.0, 0.8, 0.0]},
    "paper": {"overlaps": [1.0, 0.95, 0.9, 0.85, 0.8, 0.75, 0.7, 0.65, 0.6, 0.5, 0.4, 0.2, 0.0]},
}
TOLERANCE = {"rtol": 1e-6, "atol": 1e-9}
NOISE = ()
BOUND_BY = "dense diagonalization of the history Hamiltonian: seven steps times eight qubits"

PAULI_X = np.array([[0.0, 1.0], [1.0, 0.0]])
PAULI_Z = np.diag([1.0, -1.0])
ALICE = (PAULI_Z, PAULI_X)
BOB = ((PAULI_Z + PAULI_X) / np.sqrt(2), (PAULI_Z - PAULI_X) / np.sqrt(2))
SETTING_X, SETTING_Y, CHARLIE, CHARLIE_RECORD, CHARLIE_COPY, DEBBIE, DEBBIE_RECORD, DEBBIE_COPY = range(8)
QUBITS = 8
DIM = 2 ** QUBITS


def rotation(angle: float) -> np.ndarray:
    return np.array([[np.cos(angle / 2), -np.sin(angle / 2)], [np.sin(angle / 2), np.cos(angle / 2)]])


@equation("eq:chsh")
def chsh(frequencies: np.ndarray) -> float:
    """S = E(0,0) + E(0,1) + E(1,0) - E(1,1), with E(x, y) = sum over a, b of (-1)^(a+b) P(a, b | x, y), the
    frequencies indexed [x, y, a, b] (eq:chsh)."""
    e = {(x, y): sum((-1) ** (a + b) * frequencies[x, y, a, b] for a in (0, 1) for b in (0, 1))
         for x in (0, 1) for y in (0, 1)}
    return float(e[0, 0] + e[0, 1] + e[1, 0] - e[1, 1])


def bounds() -> dict:
    """The bound over every assignment of outcomes, the top of the Bell operator, its top state, and which readings
    commute: a part's own, or readings of different parts."""
    local = max((-1) ** (a0 + b0) + (-1) ** (a0 + b1) + (-1) ** (a1 + b0) - (-1) ** (a1 + b1)
                for a0, a1, b0, b1 in itertools.product((0, 1), repeat=4))
    operator = sum(sign * np.kron(ALICE[x], BOB[y]) for (x, y), sign in {(0, 0): 1, (0, 1): 1, (1, 0): 1, (1, 1): -1}.items())
    levels, states = np.linalg.eigh(operator)
    across = all(np.allclose(np.kron(a, np.eye(2)) @ np.kron(np.eye(2), b), np.kron(np.eye(2), b) @ np.kron(a, np.eye(2)))
                 for a in ALICE for b in BOB)
    return {"local": int(local), "quantum": float(levels[-1]), "state": states[:, -1],
            "own_commute": bool(np.allclose(ALICE[0] @ ALICE[1], ALICE[1] @ ALICE[0])), "across_commute": bool(across)}


def gluing(state: np.ndarray) -> dict:
    """The largest share of the frequencies one classical description carries: weights on the sixteen assignments
    of outcomes, below every frequency, summed as high as they go; the constraints' multipliers are the dual."""
    project = lambda reading, outcome: (np.eye(2) + (-1) ** outcome * reading) / 2
    assignments = list(itertools.product((0, 1), repeat=4))
    rows, frequencies = [], []
    for x, y, a, b in itertools.product((0, 1), repeat=4):
        frequencies.append(float(np.real(state.conj() @ np.kron(project(ALICE[x], a), project(BOB[y], b)) @ state)))
        rows.append([1.0 if (g[x] == a and g[2 + y] == b) else 0.0 for g in assignments])
    rows, frequencies = np.array(rows), np.array(frequencies)
    weights = cp.Variable(len(assignments), nonneg=True)
    below = rows @ weights <= frequencies
    share = cp.Problem(cp.Maximize(cp.sum(weights)), [below]).solve()
    multipliers = np.asarray(below.dual_value, dtype=float)
    return {"share": float(share), "dual": float(frequencies @ multipliers),
            "certificate": np.round(multipliers / multipliers.max(), 6).reshape(4, 4).tolist()}


def _on(qubit: int, single: np.ndarray) -> np.ndarray:
    out = np.array([[1.0]])
    for q in range(QUBITS):
        out = np.kron(out, single if q == qubit else np.eye(2))
    return out


def _controlled(control: int, unitary: np.ndarray) -> np.ndarray:
    one = _on(control, np.diag([0.0, 1.0]))
    return np.eye(DIM) - one + one @ unitary


def protocol(overlap: float) -> list[np.ndarray]:
    """The steps in order. Charlie measures Z and Debbie (Z + X)/sqrt 2, each writing the outcome in a record; each
    record is copied outside with the copy's two states overlapping by `overlap`; with its setting 1, Alice undoes
    Charlie's measurement and measures X, and Bob undoes Debbie's and measures (Z - X)/sqrt 2."""
    def measures(system, record, angle):
        return _on(system, rotation(angle)) @ _controlled(system, _on(record, PAULI_X)) @ _on(system, rotation(-angle))

    charlie, debbie = measures(CHARLIE, CHARLIE_RECORD, 0.0), measures(DEBBIE, DEBBIE_RECORD, np.pi / 4)
    copy = lambda record, target: _controlled(record, _on(target, rotation(2 * np.arccos(overlap))))
    alice = _controlled(SETTING_X, _on(CHARLIE, rotation(-np.pi / 2)) @ charlie)
    bob = _controlled(SETTING_Y, _on(DEBBIE, rotation(np.pi / 4)) @ debbie)
    return [charlie, debbie, copy(CHARLIE_RECORD, CHARLIE_COPY), copy(DEBBIE_RECORD, DEBBIE_COPY), alice, bob]


def start() -> np.ndarray:
    """The settings each in (|0> + |1>)/sqrt 2, the pair (|00> + |11>)/sqrt 2, every record and copy at 0."""
    out = np.zeros(DIM)
    for bits in itertools.product((0, 1), repeat=QUBITS):
        if any(bits[q] for q in (CHARLIE_RECORD, CHARLIE_COPY, DEBBIE_RECORD, DEBBIE_COPY)):
            continue
        if bits[CHARLIE] == bits[DEBBIE]:
            out[int("".join(map(str, bits)), 2)] = 0.5 / np.sqrt(2)
    return out


@equation("eq:history")
def history_hamiltonian(steps: list[np.ndarray], first: np.ndarray) -> np.ndarray:
    """H = |0><0| (x) (I - |psi0><psi0|) + sum over t of (|t><t| + |t-1><t-1|)/2 (x) I - (|t><t-1| (x) U_t + h.c.)/2,
    on the step register and the rest: zero exactly on the history of the steps from psi0 (eq:history)."""
    count, dim = len(steps), len(first)
    step = lambda i, j: np.outer(np.eye(count + 1)[i], np.eye(count + 1)[j])
    h = np.kron(step(0, 0), np.eye(dim) - np.outer(first, first.conj()))
    for t, unitary in enumerate(steps, start=1):
        h = h + 0.5 * np.kron(step(t, t) + step(t - 1, t - 1), np.eye(dim))
        h = h - 0.5 * (np.kron(step(t, t - 1), unitary) + np.kron(step(t - 1, t), unitary.conj().T))
    return h


def frequencies(whole: np.ndarray) -> np.ndarray:
    """P(a, b | x, y), conditioned on the last step and on the settings: Alice's outcome is Charlie's record if
    x = 0 and her qubit if x = 1, Bob's likewise."""
    last = whole[-1] / np.linalg.norm(whole[-1])
    p = np.zeros((2, 2, 2, 2))
    for index, amplitude in enumerate(last):
        bits = [(index >> (QUBITS - 1 - q)) & 1 for q in range(QUBITS)]
        x, y = bits[SETTING_X], bits[SETTING_Y]
        a = bits[CHARLIE_RECORD] if x == 0 else bits[CHARLIE]
        b = bits[DEBBIE_RECORD] if y == 0 else bits[DEBBIE]
        p[x, y, a, b] += abs(amplitude) ** 2
    return p / p.sum(axis=(2, 3), keepdims=True)


def run(config: dict) -> dict:
    found = bounds()
    glued = gluing(found["state"])
    threshold = 2 ** 0.75 - 1
    first = start()
    rows = []
    for overlap in [*config["overlaps"], threshold]:
        steps = protocol(overlap)
        levels, states = np.linalg.eigh(history_hamiltonian(steps, first))
        whole = states[:, 0].reshape(len(steps) + 1, DIM)
        p = frequencies(whole)
        signal = max(max(np.abs(p[x, 0].sum(axis=1) - p[x, 1].sum(axis=1)).max() for x in (0, 1)),
                     max(np.abs(p[0, y].sum(axis=0) - p[1, y].sum(axis=0)).max() for y in (0, 1)))
        rows.append({"overlap": float(overlap), "chsh": chsh(p), "closed_form": float((1 + overlap) ** 2 / np.sqrt(2)),
                     "energy": float(levels[0]), "gap": float(levels[1] - levels[0]), "signal": float(signal)})
    return {"local": found["local"], "quantum": found["quantum"], "own_commute": found["own_commute"],
            "across_commute": found["across_commute"], "gluing": glued, "threshold": threshold,
            "steps": len(protocol(1.0)), "qubits": QUBITS, "protocol": rows}
