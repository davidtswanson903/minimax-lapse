"""Figures: drawn from stored results by the work's own `figures` module, if it has one.

    python -m repokit.figures
"""
from __future__ import annotations

import importlib
import sys

from . import paths


STYLE = {"font.size": 9, "axes.grid": True, "grid.alpha": 0.3, "legend.frameon": False}


def style() -> None:
    """The figures' look: matplotlib's ten default colours for the series, a 9-point font and a light grid. A work's
    figures keep a marker or a line style on every series too, so that a figure still reads in grey."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update(STYLE)
    plt.rcParams["axes.prop_cycle"] = matplotlib.rcParamsDefault["axes.prop_cycle"]


def main() -> int:
    try:
        module = importlib.import_module(f"{paths.package()}.figures")
    except ModuleNotFoundError:
        print("the work has no figures module; nothing to draw")
        return 0
    style()
    for path in module.draw():
        print(f"drew {paths.rel(path)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
