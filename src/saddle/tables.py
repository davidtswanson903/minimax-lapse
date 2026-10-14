"""Tables for the paper, built from stored results only."""
from __future__ import annotations

from repokit.tables import Sym, Table
from repokit.values import exact, inp, measured

from saddle.experiments.catalogue import distinct_below

PIECES = {"ising": "Ising", "xxz": "XXZ", "random": "random"}
LATTICES = {"ring": "ring", "chain": "open chain", "ladder": "two-leg ladder",
            "complete": "complete graph", "star": "star", "random": "random graph"}
SYMBOLS = {"g": Sym("g", "$g$"), "J": Sym("J", "$J$"), "h": Sym("b", "$b$"), "delta": Sym("Δ", r"$\Delta$")}
DEFAULTS = {"h": 0.0}
DASH = Sym("–", "--")
NONE = Sym("—", "---")


def _digits(x: float) -> int:
    """The significant digits a parameter was given with: 0.25 has two, 2.0 one."""
    body = repr(abs(float(x))).rstrip("0").rstrip(".").replace(".", "").lstrip("0")
    return max(1, len(body))


def _parameter(key: str, specs: list) -> list:
    xs = sorted({float(s.get(key, DEFAULTS.get(key, 0.0))) for s in specs})
    parts = [SYMBOLS[key], Sym(" = ", " $=$ ")]
    for i, x in enumerate(xs):
        parts += ([", "] if i else []) + [inp(f"table.{key}", x, sig=_digits(x))]
    return parts


def _parameters(family: str, lattice: str, specs: list) -> list:
    if family == "random":
        parts = [exact("table.seeds", len({s["seed"] for s in specs})), " draws"]
    else:
        keys = [k for k in ("g", "J", "h", "delta") if any(k in s for s in specs)]
        parts = []
        # the parameter with the most values, run once for each setting of the others that vary, then those
        # held fixed throughout: the cell names the combinations that ran, not a product of values
        value = lambda s, k: float(s.get(k, DEFAULTS.get(k, 0.0)))
        seen = {k: sorted({value(s, k) for s in specs}) for k in keys}
        varying = [k for k in keys if len(seen[k]) > 1]
        scan = max(varying, key=lambda k: len(seen[k])) if varying else keys[0]
        held = [k for k in varying if k != scan]
        groups: dict = {}
        for s in specs:
            groups.setdefault(tuple(value(s, k) for k in held), []).append(s)
        for i, (values, group) in enumerate(groups.items()):
            parts += (["; "] if i else []) + _parameter(scan, group)
            for key, x in zip(held, values):
                if key in DEFAULTS and x == DEFAULTS[key]:
                    continue
                parts += [", ", SYMBOLS[key], Sym(" = ", " $=$ "), inp(f"table.{key}", x, sig=_digits(x))]
        for key in keys:
            if key != scan and key not in varying:
                parts += ["; "] + _parameter(key, specs)
    if any(s.get("deep_at") for s in specs):
        parts += ["; deepening the hub or a leaf" if lattice == "star" else "; deepening an end or the middle"]
    return parts


def catalogue_table(catalogue: dict) -> Table:
    kinds: dict = {}
    for row in catalogue["models"]:
        kinds.setdefault((row["spec"]["family"], row["spec"]["lattice"]), []).append(row)
    body = []
    for (family, lattice), rows in kinds.items():
        body.append(_row(PIECES[family], LATTICES[lattice], rows, _parameters(family, lattice, [r["spec"] for r in rows])))
    every = catalogue["models"]
    body += [None, _row("All", "", every, "")]
    return Table(
        name="catalogue",
        columns=("Pieces", "Lattice", "Cells", "Parameters", "Models", "Places below", "Drops out", "Weight kept"),
        align=("l", "l", "c", "p{2.45in}", "r", "r", "r", "r"),
        rows=tuple(body),
    )


def _row(pieces: str, lattice: str, rows: list, parameters) -> tuple:
    sizes = [r["spec"]["size"] for r in rows]
    cells = exact("table.cells", min(sizes)) if min(sizes) == max(sizes) else \
        [exact("table.cells", min(sizes)), DASH, exact("table.cells", max(sizes))]
    weighted = [r for r in rows if r["weighted"]]
    kept = [r["deep_weight"] for r in weighted if not r["dropped"]]
    return (
        pieces, lattice, cells, parameters,
        exact("table.models", len(rows)),
        exact("table.below", distinct_below(rows)),
        [exact("table.dropped", sum(r["dropped"] for r in weighted)), " of ", exact("table.models", len(weighted))],
        measured("table.kept", max(kept), sig=2) if kept else NONE,
    )


MARK = Sym("yes", r"$\bullet$")


def bump_sweep_table(bump: dict) -> Table:
    """The bump's sweep, ring by ring: its size and variation, the uniform ring's flattest curvature, their ratio, and
    the bump's width at half height in cells, as a share of the ring and against (n^2/eps)^(1/4)."""
    from saddle.experiments.bump import RESPONSE
    rows = tuple(
        (exact("table.cells", r["cells"]), inp("table.eps", r["eps"], sig=1), measured("table.flat", r["flat"], sig=2),
         measured("table.response", r["response"], sig=2), exact("table.width", r["width"]),
         measured("table.share", r["width"] / r["cells"], sig=2),
         measured("table.scale", r["width"] / (r["cells"] ** 2 / r["eps"]) ** 0.25, sig=2),
         MARK if r["response"] >= RESPONSE else NONE)
        for r in sorted(bump["sweep"], key=lambda r: (r["cells"], r["eps"])))
    return Table(
        name="bump-sweep",
        columns=(Sym("n", "$n$"), Sym("ε", r"$\epsilon$"), "Flattest curvature", "Ratio", "Width",
                 Sym("width/n", "Width$/n$"), Sym("width/(n²/ε)^¼", r"Width$/(n^2/\epsilon)^{1/4}$"), "Fitted"),
        align=("r", "r", "r", "r", "r", "r", "r", "c"),
        rows=rows,
    )


def large_rings_table(large: dict) -> Table:
    """The Ising drop-out by free fermions, ring by ring: the depth, the bracket certifying the other cells' optimum,
    and the fits each size enters."""
    starts = large["_provenance"]["config"]["fit_from"]
    rows = tuple(
        (exact("table.cells", r["cells"]), measured("table.depth", r["depth"], sig=4),
         measured("table.bracket", r["gap"], sig=1), *[MARK if r["cells"] >= s else NONE for s in starts])
        for r in large["rings"])
    return Table(
        name="large-rings",
        columns=(Sym("n", "$n$"), Sym("μ_c", r"Depth $\mu_c$"), "Bracket",
                 *[["Fit from ", exact("table.from", s)] for s in starts]),
        align=("r", "r", "r", *["c"] * len(starts)),
        rows=rows,
    )


def tables(results: dict, values: dict) -> list[Table]:
    out = [catalogue_table(results["catalogue"])] if "catalogue" in results else []
    out += [bump_sweep_table(results["bump"])] if "bump" in results else []
    out += [large_rings_table(results["largedrop"])] if "largedrop" in results else []
    return out
