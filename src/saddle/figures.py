"""Figures, drawn from stored results only."""
from __future__ import annotations

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.ticker import LogLocator, MaxNLocator, NullFormatter, ScalarFormatter  # noqa: E402

from repokit import paths, sources  # noqa: E402


def draw() -> list:
    return [_deep_cell(), _size(), _simulation(), _timeless(), _critical(), _patch(), _records(), _bump()]


def _save(fig, name: str):
    out = paths.BUILD / "figures" / name
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, metadata={"CreationDate": None})
    plt.close(fig)
    return out


def _size():
    """(a) The depth at which a deepened cell drops out of an Ising ring, bracketed, by
    size; (b) the deepened cell's weight on XXZ rings against its depth, by size."""
    scaling = sources.results()["scaling"]
    fig, (top, bottom, product) = plt.subplots(1, 3, figsize=(7.0, 2.3))
    cells = [e["cells"] for e in scaling["ising"]]
    top.plot(cells, [e["depth"] for e in scaling["ising"]], marker="o", markersize=3, linewidth=1, color="C0",
             label="state side")
    larger = sources.results()["depths"]
    top.plot([r["cells"] for r in larger["rings"]], [r["depth"] for r in larger["rings"]], marker="o", markersize=3,
             linestyle="none", color="C1", markerfacecolor="white", label="weighting side")
    free = sources.results()["largedrop"]
    top.plot([r["cells"] for r in free["rings"]], [r["depth"] for r in free["rings"]], marker="s", markersize=2.5,
             linestyle="none", color="C2", label="free fermions")
    top.axhline(free["fits"][-1]["limit"], color="k", linestyle="--", linewidth=0.8, label=r"fitted limit $\mu_\infty$")
    top.set_xscale("log")
    top.legend(fontsize=6)
    top.set_xlabel("cells in the Ising ring")
    top.set_ylabel(r"drop-out depth $\mu_c$")
    top.set_title("(a)", loc="left", fontsize=8)
    markers = iter(("o", "s", "^", "D", "v", "P"))
    for entry in scaling["xxz"]:
        depths = [p["depth"] for p in entry["curve"]]
        weights = [p["weight"] for p in entry["curve"]]
        odd = entry["cells"] % 2 == 1
        style = {"marker": next(markers), "markersize": 3, "linewidth": 1, "markerfacecolor": "white" if odd else None,
                 "linestyle": "-" if odd else "--", "color": f"C{entry['cells'] - 3}"}
        flags = [p["dropped"] for p in entry["curve"]]
        cut = flags.index(True) if True in flags else len(flags)  # no line across the fall to zero
        bottom.plot(depths[:cut], weights[:cut], label=f"{entry['cells']} cells", **style)
        if cut < len(flags):  # the weightless points drawn on top, where the other curves would hide them
            bottom.plot(depths[cut:], weights[cut:], zorder=5, **style)
        if not odd:
            kept = [(d, w * d) for d, w in zip(depths, weights) if d > 0]
            product.plot([d for d, _ in kept], [x for _, x in kept], label=f"{entry['cells']} cells", **style)
    bottom.set_xlabel(r"depth $\mu$ of the deepened cell")
    bottom.set_ylabel(r"deep cell's weight $N^*_0$, XXZ ring")
    bottom.set_ylim(bottom=-0.01)
    bottom.legend(frameon=False, fontsize=7, loc="upper right")
    bottom.set_title("(b)", loc="left", fontsize=8)
    product.set_xlabel(r"depth $\mu$ of the deepened cell")
    product.set_ylabel(r"$\mu\,N^*_0$, even XXZ ring")
    product.set_ylim(bottom=0)
    product.legend(frameon=False, fontsize=7)
    product.set_title("(c)", loc="left", fontsize=8)
    fig.tight_layout()
    return _save(fig, "size.pdf")


def _deep_cell():
    """The deep cell's weight against its depth, the closed-form depth of Eq. depth dashed and the depth the readings
    use marked."""
    ring = sources.results()["ring"]
    depths = [s["depth"] for s in ring["sweep"]]
    weights = [s["weights"][0] for s in ring["sweep"]]
    edge = next(r["depth"] for r in sources.results()["depths"]["rings"] if r["cells"] == ring["cells"])
    shallow = ring["_provenance"]["config"]["shallow"]
    fig, ax = plt.subplots(figsize=(3.4, 2.2))
    ax.plot(depths, weights, marker="o", markersize=3, linewidth=1, color="C0")
    ax.axvline(edge, linestyle="--", linewidth=0.8, color="k", label="closed-form depth")
    read = next(s for s in ring["sweep"] if abs(s["depth"] - shallow) < 1e-9)
    ax.plot([shallow], [read["weights"][0]], marker="o", markersize=7, markerfacecolor="none", color="C3",
            linestyle="none", label="the readings' depth")
    ax.legend(frameon=False, fontsize=7)
    ax.set_xlabel(r"depth $\mu$ of the deepened cell")
    ax.set_ylabel(r"deep cell's weight $N^*_0$")
    ax.set_xlim(0, max(depths))
    ax.set_ylim(bottom=-0.005)
    fig.tight_layout()
    return _save(fig, "deep-cell-weight.pdf")


def _simulation():
    """(a) Probe clocks on the ring with its first cell deepened: each runs at its cell's weight,
    and past the drop-out the deep cell's clock stops. (b) Probe thermometers at every depth
    simulated, against the inverse weight: all on Tolman's line."""
    import numpy as np
    sim = sources.results()["simulation"]
    depths = sim["_provenance"]["config"]["depths"]
    rings = {r["depth"]: r for r in sim["rings"]}
    times = sim["times"]
    fig, (clocks, tolman) = plt.subplots(2, 1, figsize=(3.4, 3.6))
    shallow, deep = rings[depths[1]], rings[depths[2]]
    for x, phase in enumerate(shallow["phases"]):
        clocks.plot(times, phase, color="C3" if x == 0 else "C0", linewidth=1.4 if x == 0 else 0.8,
                    label="deepened cell" if x == 0 else ("the other cells" if x == 1 else None))
    clocks.plot(times, deep["phases"][0], color="C3", linestyle="--", linewidth=1.2, label="deepened cell, past its drop-out")
    clocks.set_xlabel("the clock's reading $t$")
    clocks.set_ylabel(r"probe clock phase $\phi_x$")
    clocks.legend(frameon=False, fontsize=7)
    clocks.set_title("(a)", loc="left", fontsize=8)
    markers = iter(("o", "s", "^", "D"))
    reach = []
    for depth in depths:
        ring = rings[depth]
        points = [(1.0 / w, T) for w, T in zip(ring["weights"], ring["temperatures"]) if T is not None and w > 0]
        reach += [p[0] for p in points]
        tolman.loglog([p[0] for p in points], [p[1] for p in points], marker=next(markers), markersize=4,
                      linestyle="none", color=f"C{depths.index(depth)}", markerfacecolor="white", label=f"depth {depth:g}")
    line = np.logspace(np.log10(min(reach)) - 0.02, np.log10(max(reach)) + 0.02, 20)
    tolman.loglog(line, line, color="k", linestyle="--", linewidth=0.8, label=r"$\beta T_x = 1/N^*_x$")
    # within a decade a log axis labels its minor ticks too, and they collide: label 1 and 3 times each
    # power of ten, as plain numbers, and leave the other minor ticks bare
    for axis in (tolman.xaxis, tolman.yaxis):
        axis.set_major_locator(LogLocator(base=10, subs=(1.0, 3.0)))
        axis.set_major_formatter(ScalarFormatter())
        axis.set_minor_formatter(NullFormatter())
    tolman.set_xlabel(r"inverse weight $1/N^*_x$")
    tolman.set_ylabel(r"thermometer reading $\beta T_x$")
    tolman.legend(frameon=False, fontsize=7)
    tolman.set_title("(b)", loc="left", fontsize=8)
    fig.tight_layout()
    return _save(fig, "simulation.pdf")


def _timeless():
    """(a) The largest gap between a probe's rate, read against the clock inside, and its weight, against the
    clock's resolution: within half the resolution, and falling with it. (b) Each probe's rate at the finest clock
    against its weight: on the diagonal, the dropped cell's probe at zero."""
    import numpy as np
    timeless = sources.results()["timeless"]
    fig, (gaps, rates) = plt.subplots(1, 2, figsize=(7.0, 2.5))
    markers = iter(("o", "s"))
    for ring in timeless["rings"]:
        marker = next(markers)
        clocks = ring["clocks"]
        gaps.loglog([c["resolution"] for c in clocks], [max(c["gap"], 1e-16) for c in clocks], marker=marker,
                    markersize=4, linewidth=1, color=f"C{timeless['rings'].index(ring)}", markerfacecolor="white", label=f"depth {ring['depth']:g}")
        rates.plot(ring["weights"], clocks[-1]["rates"], marker=marker, markersize=4, linestyle="none", color=f"C{timeless['rings'].index(ring)}",
                   markerfacecolor="white", label=f"depth {ring['depth']:g}")
    span = np.array([min(c["resolution"] for r in timeless["rings"] for c in r["clocks"]),
                     max(c["resolution"] for r in timeless["rings"] for c in r["clocks"])])
    gaps.loglog(span, span / 2, color="k", linestyle="--", linewidth=0.8, label="half the resolution")
    gaps.set_xlabel(r"the clock's resolution, in units of $\omega$")
    gaps.set_ylabel("largest gap, rate to weight")
    gaps.legend(frameon=False, fontsize=7)
    gaps.set_title("(a)", loc="left", fontsize=8)
    top = max(max(r["weights"]) for r in timeless["rings"]) * 1.1
    rates.plot([0, top], [0, top], color="k", linestyle="--", linewidth=0.8, label="rate = weight")
    rates.set_xlabel(r"weight $N^*_x$")
    rates.set_ylabel("probe's rate per unit reading")
    rates.legend(frameon=False, fontsize=7)
    rates.set_title("(b)", loc="left", fontsize=8)
    fig.tight_layout()
    return _save(fig, "timeless.pdf")


def _critical():
    """(a) The common value's size-set part on uniform critical rings against one over the size squared: the
    conformal line. (b) An interval's entanglement entropy on the largest ring against the chord: slope c."""
    import numpy as np
    critical = sources.results()["critical"]
    fig, (size, entropy) = plt.subplots(1, 2, figsize=(7.0, 2.5))
    inverse = np.array([r["cells"] ** -2.0 for r in critical["rings"]])
    part = np.array([r["common"] - critical["per_place"] for r in critical["rings"]])
    size.plot(inverse, part, marker="o", markersize=4, linestyle="none", color="C0", markerfacecolor="white",
              label="computed")
    line = np.linspace(0, inverse.max() * 1.05, 20)
    size.plot(line, critical["size_coefficient"] * line, color="k", linestyle="--", linewidth=0.8, label=r"$-\pi c v/(6n^2)$")
    size.xaxis.set_major_locator(MaxNLocator(4))
    size.set_xlabel(r"$1/n^2$")
    size.set_ylabel(r"size-set part $e^*(n) - e_\infty$")
    size.legend(frameon=False, fontsize=7)
    size.set_title("(a)", loc="left", fontsize=8)
    held = critical["entropy"]
    cells = held["cells"]
    chord = np.array([np.log(cells / np.pi * np.sin(np.pi * n / cells)) / 3 for n in held["lengths"]])
    entropies = np.array(held["entropies"])
    entropy.plot(chord, entropies, marker="o", markersize=3, linestyle="none", color="C0", markerfacecolor="white",
                 label=f"ring of {cells} cells")
    offset = float(np.mean(entropies - held["charge"] * chord))
    entropy.plot(chord, held["charge"] * chord + offset, color="k", linestyle="--", linewidth=0.8, label="slope $c$")
    entropy.set_xlabel(r"chord $\frac{1}{3}\ln[\frac{n}{\pi}\sin\frac{\pi\ell}{n}]$")
    entropy.set_ylabel(r"entanglement entropy $S(\ell)$")
    entropy.legend(frameon=False, fontsize=7)
    entropy.set_title("(b)", loc="left", fontsize=8)
    fig.tight_layout()
    return _save(fig, "critical.pdf")
def _patch():
    """(a) The longest interval's weights, scaled by the profile's peak, against the position across it: with their
    longer-range terms on the parabola, from the nearest terms alone above it. (b) The coefficient against length."""
    import numpy as np
    stored = sources.results()["patch"]
    two_pi = stored["two_pi_over_speed"]
    fig, (weights, coefficient) = plt.subplots(1, 2, figsize=(7.0, 2.5))
    last = stored["intervals"][-1]
    length, x = last["length"], np.array(last["x"])
    peak = two_pi * length / 4
    weights.plot(x / length, np.array(last["near"]) / peak, marker="o", markersize=2.5, linestyle="none", color="C1",
                 markerfacecolor="white", label="nearest terms alone")
    weights.plot(x / length, np.array(last["full"]) / peak, marker="s", markersize=2.5, linestyle="none", color="C0",
                 markerfacecolor="white", label="with the longer-range terms")
    u = np.linspace(0, 1, 101)
    weights.plot(u, 4 * u * (1 - u), color="k", linestyle="--", linewidth=0.8, label=r"$\frac{2\pi}{v}\,x(\ell - x)/\ell$")
    weights.set_xlabel(r"position across the interval $x/\ell$")
    weights.set_ylabel(r"weight $\beta$, scaled by the peak")
    weights.legend(frameon=False, fontsize=7)
    weights.set_title("(a)", loc="left", fontsize=8)
    lengths = [i["length"] for i in stored["intervals"]]
    coefficient.plot(lengths, [i["near_coefficient"] for i in stored["intervals"]], marker="o", markersize=4,
                     linewidth=1, color="C1", markerfacecolor="white", label="nearest terms alone")
    coefficient.plot(lengths, [i["coefficient"] for i in stored["intervals"]], marker="s", markersize=4, linewidth=1,
                     color="C0", markerfacecolor="white", label="with the longer-range terms")
    coefficient.axhline(two_pi, color="k", linewidth=0.8, linestyle="--", label=r"$2\pi/v$")
    coefficient.set_xlabel(r"interval length $\ell$")
    coefficient.set_ylabel("coefficient of the profile")
    coefficient.legend(frameon=False, fontsize=7)
    coefficient.set_title("(b)", loc="left", fontsize=8)
    fig.tight_layout()
    return _save(fig, "patch.pdf")
def _records():
    """The CHSH value of Alice's and Bob's outcomes against the outside copy's overlap: on (1 + gamma)^2 / sqrt 2, above
    the bound of one description of absolute outcomes only while the copy is weak."""
    import numpy as np
    stored = sources.results()["records"]
    rows = sorted(stored["protocol"], key=lambda r: r["overlap"])
    # the legend sits below the axes: inside, every corner is crossed by the bound or the threshold
    fig, ax = plt.subplots(figsize=(3.4, 3.2), layout="constrained")
    overlap = np.linspace(0, 1, 101)
    ax.plot(overlap, (1 + overlap) ** 2 / np.sqrt(2), color="C0", linewidth=1, label=r"$(1+\gamma)^2/\sqrt{2}$")
    ax.plot([r["overlap"] for r in rows], [r["chsh"] for r in rows], marker="o", markersize=4, linestyle="none",
            color="C0", markerfacecolor="white", label="the whole, solved once")
    ax.axhline(stored["local"], color="k", linestyle="--", linewidth=0.8, label="absolute outcomes")
    ax.axhline(stored["quantum"], color="k", linestyle=":", linewidth=0.8, label="Tsirelson's bound")
    ax.axvline(stored["threshold"], color="C3", linestyle="-.", linewidth=0.8, label="threshold")
    ax.set_xlabel(r"overlap of the outside copy $\gamma$")
    ax.set_ylabel(r"CHSH value $S$")
    fig.legend(loc="outside lower center", ncol=2, fontsize=7, frameon=False, columnspacing=1.0, handlelength=2.2)
    return _save(fig, "records.pdf")


def _bump():
    """(a) The saddle's own weights on each ring, scaled to their peak, against the distance from the peak. (b) The
    saddle state's entanglement across each interval, less the vacuum's on the metric, with one cell as the cutoff
    and with one unit of optical length."""
    import numpy as np

    stored = sources.results()["bump"]
    fig, (weights, misses) = plt.subplots(2, 1, figsize=(3.4, 4.4))
    for ring in stored["rings"]:
        w = np.array(ring["weights"]) / max(ring["weights"])
        offset = np.arange(len(w)) - int(np.argmax(w))
        weights.plot(offset, w, label=f"{ring['cells']} cells")
        rows = ring["entanglement"]["intervals"]
        length = [r["stop"] - r["start"] for r in rows]
        misses.plot(length, [r["miss"] for r in rows], marker="o", linestyle="none", label=f"one cell, {ring['cells']}")
        misses.plot(length, [r["optical_miss"] for r in rows], marker="x", linestyle="none", label=f"optical unit, {ring['cells']}")
    weights.set_xlabel("cells from the peak")
    weights.set_ylabel("weight / peak")
    weights.legend(fontsize=7)
    weights.set_title("(a)", loc="left", fontsize=8)
    misses.axhline(0.0, color="k", linewidth=0.6)
    misses.set_xlabel("interval's length, cells")
    misses.set_ylabel("entropy less the vacuum's")
    misses.legend(fontsize=6)
    misses.set_title("(b)", loc="left", fontsize=8)
    fig.tight_layout()
    return _save(fig, "bump.pdf")
