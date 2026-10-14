"""Uniqueness of the optimal weighting, model by model, by Proposition unique's certificate.

At a computed optimal weighting the certificate asks three things: that the weighted Hamiltonian's lowest level is a
single state; which places that state holds at the common value, the only places any optimal weighting can weight;
and that the floor's Hessian, Eq. curvature, is negative definite on the weightings that move only those places and
sum to zero. Where all three hold the optimal weighting is unique. The models are the catalogue's, each before
deepening and with its deepened cell at the catalogue's depth, and the six-cell ring the readings use, at its three
depths."""
from __future__ import annotations

from saddle import model, programs
from saddle.experiments import catalogue

CONFIGS = {
    "golden": {"models": catalogue.CONFIGS["golden"]["models"], "deep": 2.0,
               "ring": {"cells": 3, "depths": [0.0, 0.2]}},
    "quick": {"models": catalogue.CONFIGS["quick"]["models"], "deep": 2.0,
              "ring": {"cells": 4, "depths": [0.0, 0.2, 1.0]}},
    "survey": {"models": catalogue.CONFIGS["survey"]["models"], "deep": 2.0,
               "ring": {"cells": 6, "depths": [0.0, 0.2, 1.0]}},
    "certification": {"models": catalogue.CONFIGS["certification"]["models"], "deep": 2.0,
                      "ring": {"cells": 6, "depths": [0.0, 0.2, 1.0]}},
}
TOLERANCE = {"rtol": 1e-4, "atol": 1e-6}
NOISE = ("gap", "curvature", "spread")
WIDE = 0.05  # a spread wider than this marks a set of optimal weightings, not a weighting fixed to a width
BOUND_BY = "double-precision dense linear algebra"


NEAR = 1e-7  # the weightings within this of the common value, over which an uncertified model's weights spread


def _certify(spec: dict) -> dict:
    """The certificate; where it fails, the widest spread of any weight over the weightings near the common value,
    found by the rate programs, which tells a weighting fixed to a width from a set."""
    dens = model.build(spec)
    upper = programs.greatest_least(dens)
    row = {"spec": spec, **programs.unique_weighting(dens, upper["weights"], upper["value"])}
    if not row["unique"]:
        bounds = [programs.rate_bounds(dens, x, tol=NEAR) for x in range(len(dens))]
        row["spread"] = max(high - low for low, high in bounds)
    return row


def run(config: dict) -> dict:
    models = []
    for spec in catalogue.specs(config):
        models.append(_certify(spec))
        models.append(_certify({**spec, "depth": config["deep"]}))
    ring = config["ring"]
    rings = [_certify({"family": "ising", "lattice": "ring", "size": ring["cells"], "g": 1.0, "J": 1.0, "depth": depth})
             for depth in ring["depths"]]
    fails: dict = {}
    for row in models:
        if row["fails"]:
            fails[row["fails"]] = fails.get(row["fails"], 0) + 1
    open_rows = [r for r in models if not r["unique"]]
    return {"models": models, "rings": rings,
            "summary": {"models": len(models), "unique": sum(r["unique"] for r in models), "fails": fails,
                        "narrow": sum(r["spread"] <= WIDE for r in open_rows),
                        "wide": sum(r["spread"] > WIDE for r in open_rows)}}
