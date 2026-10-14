"""The values the surfaces print, each read from the stored ring results."""
from __future__ import annotations

import math

from repokit.values import percent, bound, exact, inp, measured

from saddle import model
from saddle.experiments.catalogue import BELOW as BELOW_CATALOGUE, NO_WEIGHT as NO_WEIGHT_CATALOGUE
from saddle.experiments import space as space_experiment
from saddle.experiments.depths import deepest_density
from saddle.experiments.ring import NEAR as RING_NEAR
from saddle.experiments.scaling import ANISOTROPY, CURVE_FIELD
from saddle.tables import _digits
from saddle.experiments.ring import NO_WEIGHT
from saddle.experiments.signals import NO_WEIGHT as NO_WEIGHT_SIGNALS, SADDLE_WINDOW


def values(results: dict) -> list:
    ring = results["ring"]
    every = [ring["uniform"], ring["deep"], *ring["sweep"]]
    cells = ring["cells"]
    strictly_below = [s for s in ring["sweep"] if s["level"] - s["densities"][0] > 1e-6]
    return [
        exact("ringCells", cells),
        measured("groundEnergyPerCell", ring["ground_per_cell"], sig=6),
        bound("uniformMatchesGroundBound", 1e-7, worst=abs(ring["uniform"]["level"] - ring["ground_per_cell"])),
        bound("programsAgreeBound", 1e-7, worst=max(s["gap"] for s in every)),
        bound("constraintAtLevelBound", 1e-4, worst=max(s["constraint"] for s in every)),
        bound("slacknessBound", 1e-7, worst=max(max(s["slackness"], s["above"], s["normalization"]) for s in every)),
        bound("uniformWeightsBound", 1e-7, worst=max(abs(w - 1 / cells) for w in ring["uniform"]["weights"])),
        bound("noWeightBound", NO_WEIGHT, worst=max(s["weights"][0] for s in strictly_below)),
        measured("dropDepth", ring["dropout"]["depth"], sig=3),
        bound("dropToleranceBound", 1e-3, worst=abs(ring["dropout"]["depth"] - ring["dropout"]["coarse"])),
        inp("deepDepth", ring["deep"]["depth"], sig=2),
        measured("deepCellMargin", ring["deep"]["level"] - ring["deep"]["densities"][0], sig=3),
        measured("deepCellMarginMost", _deepest_margin(ring, results["depths"]), sig=3),
        measured("ringShallowDeepLeast", ring["rates"]["shallow"][0]["least"], sig=3),
        measured("ringShallowDeepGreatest", ring["rates"]["shallow"][0]["greatest"], sig=3),
        measured("deepLevel", ring["deep"]["level"], sig=6),
        bound("rateMatchesWeightBound", 1e-4, worst=max(abs(c["rate"] - c["weight"]) for r in ring["rates"].values() for c in r)),
        bound("ringSpreadBound", 5e-3, worst=max(c["greatest"] - c["least"] for r in ring["rates"].values() for c in r)),
        inp("ringShallowDepth", ring["_provenance"]["config"]["shallow"], sig=1),
        inp("ringShift", ring["_provenance"]["config"]["shift"], sig=1, note="the constant each cell's density is shifted by to read its rate"),
        inp("ringNearLevel", RING_NEAR, sig=1, note="how close to the common value a weighting's floor must be for the spread of the weights"),
        inp("randomGraphExtraLink", model.RANDOM_GRAPH_EXTRA_LINK, sig=1),
        exact("randomGraphSeedOffset", model.RANDOM_GRAPH_SEED_OFFSET),
        *_catalogue(results["catalogue"]),
        *_scaling(results["scaling"]),
        *_space(results["space"]),
        *_bump(results["bump"]),
        *_unique(results["unique"]),
        *_coldest_chain(results["bridges"]["coldest_chain"]),
        *_largedrop(results["largedrop"], results["depths"]),
        *_bump_sweep(results["bump"]),
        *_simulation(results["simulation"]),
        *_bridges(results["bridges"]),
        *_exits(results["exits"]),
        *_signals(results["signals"]),
        *_timeless(results["timeless"]),
        *_critical(results["critical"]),
        *_patch(results["patch"]),
        *_records(results["records"]),
        *_depths(results["depths"], results["scaling"]),
    ]


def _signals(signals: dict) -> list:
    crossing = [p["far"] * p["interior_weight"] for p in signals["path"]
                if p["far"] is not None and p["interior_weight"] > NO_WEIGHT_SIGNALS]
    return [
        exact("signalCells", signals["cells"]),
        measured("signalExponent", signals["exponent"], sig=3),
        measured("signalCrossing", crossing[-1], sig=3),
        bound("signalCrossingSpreadBound", 0.05, worst=(max(crossing) - min(crossing)) / max(crossing)),
        inp("signalWindow", SADDLE_WINDOW, sig=1),
        bound("signalWallBound", 1e-8, worst=signals["wall"]),
    ]


def _timeless(timeless: dict) -> list:
    shallow = timeless["rings"][0]
    first, last = shallow["clocks"][0], shallow["clocks"][-1]
    stopped = [abs(c["rates"][0]) for r in timeless["rings"] if r["weights"][0] == 0.0 for c in r["clocks"]]
    return [
        exact("timelessCells", timeless["cells"]),
        inp("timelessDepth", shallow["depth"], sig=1),
        exact("timelessLevelsFirst", first["levels"]),
        exact("timelessLevels", last["levels"]),
        exact("timelessReadings", timeless["_provenance"]["config"]["readings"]),
        measured("timelessGapFirst", first["gap"], sig=2),
        measured("timelessGapLast", last["gap"], sig=2),
        bound("timelessStoppedBound", 1e-12, worst=max(stopped)),
    ]


def _critical(critical: dict) -> list:
    sizes = [r["cells"] for r in critical["rings"]]
    gap = abs(critical["charge_size"] - critical["entropy"]["charge"])
    return [
        exact("criticalSizeFirst", sizes[0]),
        exact("criticalSizeLast", sizes[-1]),
        measured("criticalPerPlace", critical["per_place"], sig=6),
        measured("criticalSpeed", critical["speed"], sig=4),
        measured("criticalChargeSize", critical["charge_size"], sig=4),
        measured("criticalChargeEntropy", critical["entropy"]["charge"], sig=4),
        exact("criticalShortest", min(critical["entropy"]["lengths"])),
        exact("criticalLongest", max(critical["entropy"]["lengths"])),
        bound("criticalChargeGapBound", 1e-3, worst=gap),
    ]


def _patch(patch: dict) -> list:
    first, last = patch["intervals"][0], patch["intervals"][-1]
    config = patch["_provenance"]["config"]
    ratio = lambda interval, key="coefficient": interval[key] / patch["two_pi_over_speed"]
    return [
        exact("patchLengthFirst", first["length"]),
        exact("patchLengthLast", last["length"]),
        measured("patchCoefficient", last["coefficient"], sig=6),
        measured("patchTwoPi", patch["two_pi_over_speed"], sig=6),
        measured("patchNearCoefficient", last["near_coefficient"], sig=4),
        bound("patchShapeBound", 1e-4, worst=last["gap"]),
        percent("patchNearShape", last["near_gap"], sig=2),
        measured("patchRatioGapFirst", abs(ratio(first) - 1), sig=2),
        measured("patchRatioGapLast", abs(ratio(last) - 1), sig=2),
        measured("patchNearRatio", ratio(last, "near_coefficient"), sig=3),
        inp("planckHubble", config["hubble"][0], sig=4),
        inp("planckFraction", config["fraction"][0], sig=4),
        measured("desitterTemperature", patch["desitter"]["temperature"], sig=3, unit="K"),
        percent("desitterSpread", patch["desitter"]["relative_error"], sig=1),
    ]


def _records(records: dict) -> list:
    rows = sorted(records["protocol"], key=lambda r: -r["overlap"])
    glued = records["gluing"]
    return [
        exact("bellLocal", records["local"]),
        measured("bellQuantum", records["quantum"], sig=6),
        measured("gluingShare", glued["share"], sig=4),
        measured("gluingContextual", 1 - glued["share"], sig=4),
        bound("gluingDualBound", 1e-6, worst=abs(glued["share"] - glued["dual"])),
        exact("recordsQubits", records["qubits"]),
        exact("recordsSteps", records["steps"]),
        measured("recordsChshFirst", rows[0]["chsh"], sig=4),
        measured("recordsChshLast", rows[-1]["chsh"], sig=4),
        measured("recordsThreshold", records["threshold"], sig=4),
        measured("recordsHistoryGap", min(r["gap"] for r in rows), sig=2),
        bound("recordsGapBound", 1e-9, worst=max(abs(r["chsh"] - r["closed_form"]) for r in rows)),
        bound("recordsEnergyBound", 1e-12, worst=max(abs(r["energy"]) for r in rows)),
        bound("recordsSignalBound", 1e-12, worst=max(r["signal"] for r in rows)),
    ]


def _depths(depths: dict, scaling: dict) -> list:
    rings = depths["rings"]
    by_size = {r["cells"]: r["depth"] for r in rings}
    agreement = max(abs(by_size[r["cells"]] - r["depth"]) for r in scaling["ising"] if r["cells"] in by_size)
    early, late = depths["fits"]
    return [
        exact("depthsSmallest", rings[0]["cells"]),
        exact("depthsLargest", rings[-1]["cells"]),
        measured("depthsFirst", rings[0]["depth"], sig=4),
        measured("depthsLast", rings[-1]["depth"], sig=4),
        measured("depthsLastStep", rings[-1]["depth"] - rings[-2]["depth"], sig=2),
        exact("depthsFitFrom", early["from"]),
        exact("depthsFitFromLate", late["from"]),
        measured("depthsLimit", early["limit"], sig=4),
        measured("depthsLimitLate", late["limit"], sig=4),
        measured("depthsPower", early["power"], sig=2),
        bound("depthsAgreementBound", 1e-3, worst=agreement),
    ]


def _exits(exits: dict) -> list:
    return [
        measured("exitChainGroundPerCell", exits["rotation"]["ground_per_cell"], sig=4),
        measured("exitChainFloor", exits["rotation"]["floor"], sig=4),
        measured("exitGain", exits["rotation"]["floor"] - exits["rotation"]["ground_per_cell"], sig=2),
        exact("exitChainCells", exits["_provenance"]["config"]["chain_cells"]),
        bound("exitIdentityBound", 1e-9, worst=exits["degenerate"]["largest_multiplier"]),
        inp("exitTransverseField", max(b["field"] for b in exits["conservation"]), sig=1),
    ]


def _bridges(bridges: dict) -> list:
    profile, coldest = bridges["profile"], bridges["coldest"]
    first, last = coldest["steps"][0], coldest["steps"][-1]
    return [
        measured("profileSlope", profile["slope"], sig=2), measured("profileCutFraction", profile["cut_fraction"], sig=1),
        measured("profileNearRatio", profile["nearest_over_largest"], sig=2),
        exact("coldestCells", coldest["cells"]), inp("coldestDepth", coldest["depth"], sig=1),
        inp("coldestFirstGap", first["gap"], sig=1), inp("coldestLastGap", last["gap"], sig=1),
        measured("coldestFirstDistance", first["distance"], sig=2), measured("coldestLastDistance", last["distance"], sig=2),
        measured("coldestFirstTotal", first["total"], sig=2), measured("coldestLastTotal", last["total"], sig=2),
    ]


def _simulation(sim: dict) -> list:
    config = sim["_provenance"]["config"]
    rings = {r["depth"]: r for r in sim["rings"]}
    shallow, deep = rings[config["depths"][1]], rings[config["depths"][2]]
    uniform = rings[0.0]
    return [
        exact("simCells", config["cells"]),
        inp("simShallowDepth", config["depths"][1], sig=1), inp("simDeepDepth", config["depths"][2], sig=1),
        measured("deepClockShallow", shallow["clock_rates"][0], sig=2),
        measured("deepHotRatio", shallow["temperatures"][0] / uniform["temperatures"][0], sig=2),
        exact("xxzSimCells", sim["xxz"]["cells"]), inp("xxzSimDepth", sim["xxz"]["depth"], sig=1),
        measured("xxzDeepClock", sim["xxz"]["clock_rates"][0], sig=2),
        bound("clockRateGapBound", 1e-9, worst=sim["gaps"]["clock"]),
        bound("tolmanGapBound", 1e-9, worst=sim["gaps"]["tolman"]),
        bound("redshiftGapBound", 1e-12, worst=sim["gaps"]["redshift"]),
        bound("phaseGapBound", 1e-7, worst=sim["gaps"]["phase"]),
    ]


def _space(space: dict) -> list:
    config = space["_provenance"]["config"]
    widest = max(config["widths"])
    at = lambda rows: next(r for r in rows if r["width"] == widest)
    zero = [r for r in space["zero"]["rows"] if "law" in r and r["width"] == widest]
    entropy = space["entanglement"].values()
    return [
        exact("spaceCells", config["cells"]),
        inp("spaceFallMin", space_experiment.FALL_MIN, sig=1),
        inp("spaceWellMin", space_experiment.WELL_MIN, sig=1),
        inp("spaceWellWidth", space_experiment.WELL_WIDTH, sig=1),
        inp("spaceStepsDrop", space_experiment.STEPS_DROP, sig=2),
        inp("spaceStepsSmooth", space_experiment.STEPS_SMOOTH, sig=2),
        exact("spaceStart", config["start"]), exact("spaceStop", config["stop"]),
        exact("spaceZeroFrom", config["zero_from"]),
        exact("spaceZeroMid", max(config["zero_to"])), exact("spaceZeroNear", min(config["zero_to"])),
        inp("spaceWidth", widest, sig=2),
        bound("spaceWaveBound", 1e-3, worst=max(abs(at(rows)["measured"] / at(rows)["optical"] - 1) for rows in space["waves"].values())),
        bound("spaceSpeedBound", 5e-3, worst=max(at(rows)["local_speed_gap"] for rows in space["waves"].values())),
        bound("spaceZeroBound", 3e-3, worst=max(abs(r["measured"] / r["law"] - 1) for r in zero)),
        inp("spaceStepWeight", config["step"], sig=1),
        percent("spaceStepShare", sum(space["step"].values()) / 2, sig=2),
        exact("spaceEntropyCells", config["entropy_cells"]),
        bound("spaceEntropyBound", 5e-4, worst=max(e["miss"] for e in entropy)),
        measured("spacePowerLow", min(e["power"] for e in entropy), sig=4),
        measured("spacePowerHigh", max(e["power"] for e in entropy), sig=4),
    ]


def _bump(bump: dict) -> list:
    config = bump["_provenance"]["config"]
    first, last = bump["rings"][0], bump["rings"][-1]
    every = [row for ring in bump["rings"] for row in ring["entanglement"]["intervals"]]
    return [
        exact("bumpCellsFirst", first["cells"]),
        exact("bumpCellsLast", last["cells"]),
        inp("bumpEpsFirst", first["eps"], sig=1),
        inp("bumpEpsLast", last["eps"], sig=1),
        bound("bumpGapBound", 1e-11, worst=max(ring["gap"] for ring in bump["rings"])),
        exact("bumpWidthFirst", first["width"]),
        exact("bumpWidthLast", last["width"]),
        inp("bumpPacket", config["packet"], sig=2),
        inp("bumpLevelLast", config["levels"][-1], sig=2),
        bound("bumpWaveBound", 3e-3, worst=max(abs(w["measured"] / w["optical"] - 1) for w in last["waves"])),
        bound("bumpEntropyBound", 1e-4, worst=max(abs(row["miss"]) for row in every)),
        measured("bumpPowerLow", min(ring["entanglement"]["power"] for ring in bump["rings"]), sig=4),
        measured("bumpPowerHigh", max(ring["entanglement"]["power"] for ring in bump["rings"]), sig=4),
        measured("bumpCurvatureFirst", -first["curvature"]["least"], sig=2),
        measured("bumpCurvatureLast", -last["curvature"]["least"], sig=2),
        measured("bumpFlatFirst", -first["uniform_curvature"]["least"], sig=2),
        measured("bumpFlatLast", -last["uniform_curvature"]["least"], sig=2),
        measured("bumpFlatPower", math.log(first["uniform_curvature"]["least"] / last["uniform_curvature"]["least"])
                                        / math.log(last["cells"] / first["cells"]), sig=2),
    ]


def _unique(unique: dict) -> list:
    """Proposition single's certificate across the catalogue, and on the readings' ring."""
    from saddle.experiments.unique import NEAR as UNIQUE_NEAR, WIDE as UNIQUE_WIDE
    from saddle.programs import UNIQUE_CURVE, UNIQUE_GAP
    summary, open_rows = unique["summary"], [r for r in unique["models"] if not r["unique"]]
    certified = [r for r in unique["rings"] if r["unique"]]
    return [
        exact("uniqueModels", summary["models"]),
        exact("uniqueCertified", summary["unique"]),
        exact("uniqueNarrow", summary["narrow"]),
        exact("uniqueWide", summary["wide"]),
        inp("uniqueNear", UNIQUE_NEAR, sig=1, note="how close to the common value a weighting's floor must be for the spread"),
        bound("uniqueSpreadBound", UNIQUE_WIDE, worst=max(r["spread"] for r in open_rows if r["spread"] <= UNIQUE_WIDE)),
        measured("uniqueRingCurvature", min(-r["curvature"] for r in certified), sig=2),
        inp("uniqueCurve", UNIQUE_CURVE, sig=1, note="how far below zero a curvature must be to count as negative"),
        inp("uniqueGap", UNIQUE_GAP, sig=1, note="how far the next level must sit for the lowest to count as single"),
        inp("uniqueRingDepth", max(r["spec"]["depth"] for r in certified), sig=1),
    ]


def _largedrop(large: dict, depths: dict) -> list:
    """The Ising drop-out at size, by free fermions, against the sequential quadratic depths where both are computed."""
    rings = {r["cells"]: r for r in large["rings"]}
    stored = {r["cells"]: r["depth"] for r in depths["rings"]}
    shared = [n for n in rings if n in stored]
    doublings = [n for n in sorted(rings) if 2 * n in rings and 4 * n in rings and n >= 40]
    shares = [(rings[4 * n]["depth"] - rings[2 * n]["depth"]) / (rings[2 * n]["depth"] - rings[n]["depth"]) for n in doublings]
    early, late = large["fits"][0], large["fits"][-1]
    return [
        exact("largeSmallest", min(rings)), exact("largeLargest", max(rings)),
        bound("largeAgreementBound", 1e-6, worst=max(abs(rings[n]["depth"] - stored[n]) for n in shared)),
        exact("largeSharedLargest", max(shared)),
        bound("largeGapBound", 1e-4, worst=max(r["gap"] for r in rings.values())),
        measured("largeDepthLast", rings[max(rings)]["depth"], sig=4),
        measured("largeDoublingShare", sum(shares) / len(shares), sig=2),
        exact("largeFitFrom", early["from"]), exact("largeFitFromLate", late["from"]),
        measured("largeLimit", early["limit"], sig=4), measured("largeLimitLate", late["limit"], sig=4),
        measured("largePower", early["power"], sig=2),
        measured("largePowerLate", late["power"], sig=2),
        measured("largeCurvatureSix", -rings[6]["curvature"], sig=2),
    ]


def _bump_sweep(bump: dict) -> list:
    """The bump's width across sizes and variations: the law where the variation is several times the uniform ring's
    flattest curvature, and the share of the ring the bump spans where it is not."""
    from saddle.experiments.bump import RESPONSE
    sweep, law = bump["sweep"], bump["law"]
    gentle = [r for r in sweep if r["response"] < RESPONSE]
    return [
        exact("bumpSweepCellsLeast", min(r["cells"] for r in sweep)), exact("bumpSweepCellsMost", max(r["cells"] for r in sweep)),
        inp("bumpSweepEpsLeast", min(r["eps"] for r in sweep), sig=1), inp("bumpSweepEpsMost", max(r["eps"] for r in sweep), sig=1),
        exact("bumpSweepCount", len(sweep)), exact("bumpFitCount", law["fitted"]),
        inp("bumpResponse", RESPONSE, sig=1, note="the variation over the uniform ring's flattest curvature above which the law is fitted"),
        bound("bumpSweepGapBound", 1e-6, worst=max(r["gap"] for r in sweep)),
        measured("bumpCellsPower", law["cells_power"], sig=2), measured("bumpEpsPower", law["eps_power"], sig=2),
        measured("bumpScaleLeast", law["ratio_least"], sig=2), measured("bumpScaleMost", law["ratio_most"], sig=2),
        measured("bumpGentleResponse", max(r["response"] for r in gentle), sig=2),
        measured("bumpGentleShareLeast", min(r["width"] / r["cells"] for r in gentle), sig=2),
        measured("bumpGentleShareMost", max(r["width"] / r["cells"] for r in gentle), sig=2),
    ]


def _coldest_chain(chain: dict) -> list:
    first, last = chain["steps"][0], chain["steps"][-1]
    return [
        exact("coldestChainCells", chain["cells"]),
        measured("coldestChainFirstEnd", first["ends"][0], sig=3),
        measured("coldestChainLastEnd", last["ends"][0], sig=3),
    ]


def _deepest_margin(ring: dict, depths: dict) -> float:
    """How far below the common value the deep cell can sit, over the optimal states past the drop-out, at the
    ring's deep depth: e' less the least density the others' optimal pair gives it."""
    config = ring["_provenance"]["config"]
    same = next(r for r in depths["rings"] if r["cells"] == ring["cells"])
    return same["common"] - deepest_density(same["link"], same["flip"], config["field"], ring["deep"]["depth"])


def _scaling(scaling: dict) -> list:
    ising, xxz = scaling["ising"], scaling["xxz"]
    even = [e for e in xxz if e["cells"] % 2 == 0]
    return [
        exact("scalingSmallest", min(e["cells"] for e in ising + xxz)),
        exact("scalingLargest", max(e["cells"] for e in ising + xxz)),
        measured("scalingFirst", ising[0]["depth"], sig=3),
        measured("scalingLast", ising[-1]["depth"], sig=3),
        *_xxz(xxz, scaling["fields"]),
        *_odd(scaling["odd"], scaling["_provenance"]["config"]),
    ]


def _odd(odd: list, config: dict) -> list:
    """The odd XXZ rings across anisotropies: how many drop out, and the one that keeps its deepened cell."""
    kept = [e for e in odd if e["pair"]["depth"] is None]
    one = kept[0]
    return [
        inp("xxzOddDeltaLeast", min(config["anisotropies"]), sig=1),
        inp("xxzOddDeltaMost", max(config["anisotropies"]), sig=2),
        exact("xxzOddRings", len(odd)),
        exact("xxzOddDropped", len(odd) - len(kept)),
        inp("xxzKeptDelta", one["delta"], sig=2),
        exact("xxzKeptCells", one["cells"]),
        measured("xxzKeptShare", sum(one["pair"]["share"]) / 2, sig=3),
        measured("xxzKeptSpin", sum(one["pair"]["spin"]) / 2, sig=2),
        measured("xxzKeptProduct", one["product"], sig=2),
        measured("xxzKeptLift", one["lift"], sig=2),
    ]


def _fall(entry: dict) -> tuple[float, float]:
    """The sampled depths between which an odd ring's deepened cell falls to zero."""
    curve = entry["curve"]
    i = next(k for k, p in enumerate(curve) if p["dropped"])
    return curve[i - 1]["depth"], curve[i]["depth"]


def _xxz(xxz: list, fields: list) -> list:
    odd = [e for e in xxz if e["cells"] % 2 == 1]
    even = [e for e in xxz if e["cells"] % 2 == 0]
    product = lambda e: e["curve"][-1]["weight"] * e["curve"][-1]["depth"]
    ratios: dict = {}
    for f in fields:
        ratios.setdefault(f["cells"], []).append(f["product"] / f["field"])
    spread = max((max(r) - min(r)) / (sum(r) / len(r)) for r in ratios.values())
    # At a finite depth mu the law's gap carries the depth's own factor mu / (mu + b): what is left once it is divided
    # out, on the smallest even ring, where it is all of the gap, and on the largest.
    deepest = even[0]["curve"][-1]["depth"]
    left: dict = {}
    for f in fields:
        left.setdefault(f["cells"], []).append(abs(f["product"] * (deepest + f["field"]) / deepest / f["predicted"] - 1))
    return [
        bound("xxzDepthFactorBound", 1e-4, worst=max(left[min(left)])),
        percent("xxzDepthLeftLeast", min(left[max(left)]), sig=1), percent("xxzDepthLeftMost", max(left[max(left)]), sig=1),
        inp("xxzDeepest", even[0]["curve"][-1]["depth"], sig=_digits(even[0]["curve"][-1]["depth"])),
        inp("xxzAnisotropy", ANISOTROPY, sig=1, note="the XXZ rings' anisotropy, with coupling one"),
        inp("xxzFieldLeast", min(f["field"] for f in fields), sig=1),
        inp("xxzCurveField", CURVE_FIELD, sig=1, note="the field of the XXZ drop-out depths and curves"),
        inp("xxzFieldMost", max(f["field"] for f in fields), sig=1),
        exact("xxzOddFirstCells", odd[0]["cells"]), exact("xxzOddLastCells", odd[-1]["cells"]),
        measured("xxzOddFirstDepth", odd[0]["exact"], sig=3), measured("xxzOddLastDepth", odd[-1]["exact"], sig=3),
        percent("xxzLawGap", max(abs(f["product"] - f["predicted"]) / f["predicted"]
                                 for f in fields if f["field"] == min(g["field"] for g in fields)), sig=2),
        exact("xxzEvenFirstCells", even[0]["cells"]), exact("xxzEvenLastCells", even[-1]["cells"]),
        measured("xxzKappaFirst", product(even[0]), sig=2), measured("xxzKappaLast", product(even[-1]), sig=2),
        percent("xxzFieldSpread", spread, sig=1),
    ]


def _catalogue(catalogue: dict) -> list:
    summary, rows = catalogue["summary"], catalogue["models"]
    ising = [r for r in rows if r["spec"]["family"] == "ising" and r["weighted"]]
    others = [r for r in rows if r["spec"]["family"] != "ising" and r["weighted"]]
    return [
        exact("catalogueModels", summary["models"]),
        exact("catalogueFamilies", summary["families"]),
        exact("catalogueLattices", summary["lattices"]),
        exact("catalogueLargest", summary["largest"]),
        exact("catalogueSeedFirst", min(r["spec"]["seed"] for r in rows if "seed" in r["spec"])),
        exact("catalogueSeedLast", max(r["spec"]["seed"] for r in rows if "seed" in r["spec"])),
        exact("catalogueBelowPlaces", summary["below_places"]),
        inp("catalogueBelowMargin", BELOW_CATALOGUE, sig=1, note="how far below the common value a place must sit to count as below it"),
        exact("catalogueUniformRings", summary["uniform_checked"]),
        bound("catalogueGapBound", 1e-7, worst=summary["worst"]["gap"]),
        bound("catalogueNoWeightBound", NO_WEIGHT_CATALOGUE, worst=summary["worst"]["below_weight_max"]),
        inp("catalogueDeepDepth", catalogue["_provenance"]["config"]["deep"], sig=2),
        exact("catalogueIsingModels", len(ising)),
        exact("catalogueIsingDropped", sum(r["dropped"] for r in ising)),
        exact("catalogueOtherModels", len(others)),
        exact("catalogueOtherDropped", sum(r["dropped"] for r in others)),
        measured("catalogueKeptLargest", max(r["deep_weight"] for r in ising + others if not r["dropped"]), sig=2),
    ]
