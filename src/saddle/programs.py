"""The two programs whose values the saddle equates.

    least over states of the greatest local density      (a semidefinite program in the state)
    greatest over weights on the simplex of the lowest
    eigenvalue of the weighted sum of densities           (a semidefinite program in the weights)

The pieces used here are real, so real symmetric states suffice: the real part of
an optimal state is a state with the same densities. Both programs are solved by
Clarabel through CVXPY."""
from __future__ import annotations

import cvxpy as cp
import numpy as np

from repokit.equations import equation

SOLVER = cp.CLARABEL


def _solve(problem: cp.Problem, what: str) -> None:
    """Solve with Clarabel. A run that stops short of its tolerances is run again with five times the
    iterations and longer iterative refinement; one still short is kept as inaccurate, since every
    caller judges the pair by its gap and the residuals of its conditions."""
    problem.solve(solver=SOLVER)
    if problem.status == cp.OPTIMAL_INACCURATE:
        problem.solve(solver=SOLVER, max_iter=1000, iterative_refinement_max_iter=50)
    if problem.status not in (cp.OPTIMAL, cp.OPTIMAL_INACCURATE):
        raise RuntimeError(f"the {what} program ended {problem.status}")


@equation("eq:saddle")
def least_greatest(dens: list[np.ndarray]) -> dict:
    """The left side of eq:saddle: min over states rho of max_x Tr(rho D_x), as a semidefinite program."""
    size = dens[0].shape[0]
    state = cp.Variable((size, size), symmetric=True)
    top = cp.Variable()
    bounds = [cp.trace(d @ state) <= top for d in dens]
    problem = cp.Problem(cp.Minimize(top), [state >> 0, cp.trace(state) == 1, *bounds])
    _solve(problem, "state")
    rho = (state.value + state.value.T) / 2
    return {
        "value": float(top.value),
        "state": rho,
        "densities": np.array([float(np.trace(d @ rho)) for d in dens]),
    }


@equation("eq:saddle")
def greatest_least(dens: list[np.ndarray]) -> dict:
    """The right side of eq:saddle: max over weightings N in the simplex of lambda_min(sum_x N_x D_x)."""
    count, size = len(dens), dens[0].shape[0]
    weights = cp.Variable(count, nonneg=True)
    floor = cp.Variable()
    weighted = sum(weights[i] * dens[i] for i in range(count))
    problem = cp.Problem(cp.Maximize(floor), [cp.sum(weights) == 1, weighted - floor * np.eye(size) >> 0])
    _solve(problem, "weight")
    return {"value": float(floor.value), "weights": np.clip(np.asarray(weights.value, dtype=float), 0.0, None)}


def saddle(dens: list[np.ndarray]) -> dict:
    """Both programs, and the residuals of the conditions that hold at a saddle pair:
    the weighted sum holds the state at the common value; no density exceeds it;
    weight sits only where a place reaches it; the weights sum to one."""
    lower, upper = least_greatest(dens), greatest_least(dens)
    level, weights, rho = upper["value"], upper["weights"], lower["state"]
    return {
        "level": level,
        "gap": abs(lower["value"] - level),
        "weights": weights,
        "densities": lower["densities"],
        **residuals(dens, rho, weights, level, lower["densities"]),
    }


@equation("eq:conditions")
def residuals(dens: list[np.ndarray], rho: np.ndarray, weights, level: float, densities=None) -> dict:
    """How far a pair is from the four conditions of eq:conditions: (sum_x N_x D_x) rho = e rho;
    Tr(rho D_x) <= e; N_x (Tr(rho D_x) - e) = 0; sum_x N_x = 1."""
    weights = np.asarray(weights, dtype=float)
    if densities is None:
        densities = np.array([float(np.real(np.trace(rho @ d))) for d in dens])
    held = weighted(dens, weights)
    return {
        "constraint": float(np.linalg.norm(held @ rho - level * rho)),
        "above": float(max(0.0, np.max(densities - level))),
        "slackness": float(np.max(np.abs(weights * (densities - level)))),
        "normalization": float(abs(np.sum(weights) - 1.0)),
    }


@equation("eq:certificate")
def certificate(dens: list[np.ndarray], weights) -> float:
    """lambda_min(sum_x N_x D_x), a lower bound on max_x Tr(rho D_x) for every state rho (eq:certificate)."""
    return float(np.linalg.eigvalsh(weighted(dens, weights))[0])


@equation("eq:evolution")
def weighted(dens: list[np.ndarray], weights) -> np.ndarray:
    """A_N = sum_x N_x D_x, the weighted Hamiltonian (eq:evolution)."""
    return sum(float(w) * d for w, d in zip(weights, dens))


@equation("eq:evolution")
def weighted_from_pieces(place_pieces, link_pieces, links, weights) -> np.ndarray:
    """A_N = sum_x N_x h_x + sum over links {x, y} of (N_x + N_y)/2 h_xy: each place's own piece at its
    weight and each link at the mean weight of its ends (eq:evolution)."""
    out = sum(float(weights[x]) * piece for x, piece in enumerate(place_pieces))
    return out + sum(0.5 * (float(weights[a]) + float(weights[b])) * piece for piece, (a, b) in zip(link_pieces, links))


@equation("eq:rates")
def rate_bounds(dens: list[np.ndarray], place: int, tol: float = 1e-7) -> tuple[float, float]:
    """The one-sided rates of eq:rates, as (left, right): the least and the greatest weight N_y over
    the optimal weightings Delta*, which are the rates at 0- and at 0+, found as two semidefinite
    programs over weightings whose weighted ground energy is within tol of the common value."""
    level = greatest_least(dens)["value"]
    dim = dens[0].shape[0]
    # A common value solved only to reduced accuracy can sit just above the true one, leaving
    # no weighting within tol of it; the set is then widened, as drop_depth widens its own.
    for widen in (1.0, 10.0, 100.0):
        out = []
        for sense in (cp.Minimize, cp.Maximize):
            weights = cp.Variable(len(dens), nonneg=True)
            held = sum(weights[i] * dens[i] for i in range(len(dens)))
            problem = cp.Problem(sense(weights[place]),
                                 [cp.sum(weights) == 1, held - (level - widen * tol) * np.eye(dim) >> 0])
            try:
                problem.solve(solver=SOLVER)
            except cp.error.SolverError:  # a set that is nearly one point can stop the solver; it is widened
                continue
            if problem.status in ("optimal", "optimal_inaccurate"):
                out.append(float(problem.value))
        if len(out) == 2:
            return out[0], out[1]
    raise RuntimeError("the rate programs failed at every widening")


@equation("eq:gain")
def reweighting_gain(dens: list[np.ndarray]) -> float:
    """Gamma = e* - lambda_min(H)/n with H = sum_x D_x: what the optimal weighting adds to the uniform one,
    which is H/n. Never negative, and zero exactly when uniform weights are optimal (eq:gain)."""
    return float(saddle(dens)["level"] - np.linalg.eigvalsh(sum(dens))[0] / len(dens))


@equation("eq:curvature")
def floor_hessian(dens: list[np.ndarray], weights) -> tuple[np.ndarray, float]:
    """H_xy = -2 sum_{j != 0} Re <0|D_x|j><j|D_y|0> / (E_j - E_0): the floor's second derivatives at a weighting whose
    weighted Hamiltonian has a single lowest state |0> (eq:curvature), with the gap E_1 - E_0 that makes it hold."""
    energies, states = np.linalg.eigh(weighted(dens, weights))
    amplitudes = np.array([states.conj().T @ d @ states[:, 0] for d in dens])[:, 1:]
    with np.errstate(divide="ignore", invalid="ignore"):  # a level not single has no Hessian; the gap says so
        hessian = -2 * np.real((amplitudes / (energies[1:] - energies[0])) @ amplitudes.conj().T)
    return hessian, float(energies[1] - energies[0])


UNIQUE_GAP = 1e-6  # a lowest level split from the next by less than this is not counted single
UNIQUE_AT = 1e-6  # a place within this of the common value counts as at it, which only widens what is checked
UNIQUE_CURVE = 1e-6  # a curvature nearer zero than this is not counted negative


def unique_weighting(dens: list[np.ndarray], weights, level: float) -> dict:
    """Proposition unique's certificate at a computed optimal weighting: its weighted Hamiltonian's lowest level is a
    single state; the places that state holds at the common value are the only ones any optimal weighting can
    weight; and on the weightings that move only those places, summing to zero, the floor's Hessian is negative
    definite. Where all three hold, the optimal weighting is unique. Reports which fails, if one does."""
    hessian, gap = floor_hessian(dens, weights)
    ground = np.linalg.eigh(weighted(dens, weights))[1][:, 0]
    at = [x for x, d in enumerate(dens) if level - float(np.real(ground.conj() @ d @ ground)) <= UNIQUE_AT]
    curvature = None  # one place at the value leaves nothing to move, and a level not single has no Hessian
    if len(at) > 1 and gap > UNIQUE_GAP:
        basis = np.linalg.qr(np.vstack([np.ones(len(at)), np.eye(len(at))[:-1]]).T)[0][:, 1:]
        curvature = float(np.linalg.eigvalsh(basis.T @ hessian[np.ix_(at, at)] @ basis).max())
    if gap <= UNIQUE_GAP:
        reason = "lowest level not single"
    elif curvature is not None and curvature >= -UNIQUE_CURVE:
        reason = "flat"
    else:
        reason = None
    return {"gap": gap, "at": len(at), "curvature": curvature, "unique": reason is None, "fails": reason}


@equation("eq:uniform")
def uniform_value(dens: list[np.ndarray]) -> float:
    """lambda_min(H)/n, the common value of a ring a rotation maps into itself (eq:uniform)."""
    return float(np.linalg.eigvalsh(sum(dens))[0]) / len(dens)


@equation("eq:dropout")
def drop_depth(dens: list[np.ndarray], deep: int, spin: np.ndarray, tol: float = 1e-9) -> float | None:
    """mu_c of eq:dropout, the depth from which a place deepened as D_deep - depth * spin carries no weight: the least depth at which some state optimal for the other places holds
    it at or below their common value. That is the least ratio (Tr(rho D_deep) - e') /
    Tr(rho spin) over the other places' optimal states with positive spin, found as one
    semidefinite program after the substitution sigma = t rho. None when no optimal state
    of the others gives the place a positive spin, so the program is infeasible: then it never
    drops out.

    Eq. dropout assumes the place has weight at depth zero, which makes mu_c positive. Without
    that, the least ratio is at most zero, or the program is unbounded below when the place
    already sits below e' in an optimal state of the others; either way it carries no weight from
    depth zero on, and the function returns 0.

    e' is the others' primal value, an upper bound on their common value, so the set of
    their optimal states is never empty; tol widens it slightly, and the depth's
    sensitivity to tol is recorded by the callers."""
    others = [d for i, d in enumerate(dens) if i != deep]
    level = least_greatest(others)["value"]
    dim = dens[0].shape[0]
    # A face that is nearly one point can make an interior-point step fail; then the set is
    # widened a little and solved again. The callers' fine-against-coarse bound covers this.
    for widen in (1.0, 3.0, 10.0):
        sigma, t = cp.Variable((dim, dim), symmetric=True), cp.Variable(nonneg=True)
        constraints = [sigma >> 0, cp.trace(sigma) == t, cp.trace(sigma @ spin) == 1]
        constraints += [cp.trace(sigma @ d) <= (level + widen * tol) * t for d in others]
        problem = cp.Problem(cp.Minimize(cp.trace(sigma @ dens[deep]) - level * t), constraints)
        try:
            problem.solve(solver=SOLVER)
        except BaseException as error:  # the solver's internal failures arrive as a panic
            if isinstance(error, (KeyboardInterrupt, SystemExit)):
                raise
            continue
        if problem.status in ("infeasible", "infeasible_inaccurate"):
            return None
        if problem.status in ("unbounded", "unbounded_inaccurate"):
            return 0.0
        if problem.value is not None:
            return max(0.0, float(problem.value))
    raise RuntimeError("the drop-out program failed at every widening")
