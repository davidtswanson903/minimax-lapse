"""LaTeX the reader's tree needs: the values file, and the keys and claims the paper marks. The composer's writers
of the paper's other generated parts extend this module from the private root."""
from __future__ import annotations

import pathlib
import re
import unicodedata

from . import paths

from .render import banner
from .values import Value

COMMENT = re.compile(r"(?<!\\)%.*")
ABSTRACT = re.compile(r"\\begin\{abstract\}(.*?)\\end\{abstract\}", re.S)
CLAIM = re.compile(r"%\s*claim:\s*([^\n]+)")
CITE = re.compile(r"\\cite[a-zA-Z]*\*?(?:\[[^\]]*\]){0,2}\{([^}]*)\}")
VAL = re.compile(r"\\val\{([^}]*)\}")
MATH = re.compile(r"\$([^$]*)\$")

GREEK = {
    "alpha": "α", "beta": "β", "gamma": "γ", "delta": "δ", "epsilon": "ε", "varepsilon": "ε",
    "zeta": "ζ", "eta": "η", "theta": "θ", "kappa": "κ", "lambda": "λ", "mu": "μ", "nu": "ν",
    "xi": "ξ", "pi": "π", "rho": "ρ", "sigma": "σ", "tau": "τ", "phi": "φ", "varphi": "φ",
    "chi": "χ", "psi": "ψ", "omega": "ω", "Gamma": "Γ", "Delta": "Δ", "Theta": "Θ",
    "Lambda": "Λ", "Xi": "Ξ", "Pi": "Π", "Sigma": "Σ", "Phi": "Φ", "Psi": "Ψ", "Omega": "Ω",
}
SYMBOLS = {
    "times": "×", "cdot": "·", "approx": "≈", "sim": "∼", "simeq": "≃", "le": "≤", "leq": "≤",
    "ge": "≥", "geq": "≥", "neq": "≠", "ne": "≠", "pm": "±", "infty": "∞", "partial": "∂",
    "to": "→", "rightarrow": "→", "propto": "∝", "ell": "ℓ", "hbar": "ħ", "langle": "⟨",
    "rangle": "⟩", "sum": "Σ", "int": "∫", "nabla": "∇", "ll": "≪", "gg": "≫", "in": "∈",
    "min": "min", "max": "max", "inf": "inf", "sup": "sup", "lim": "lim", "log": "log", "exp": "exp",
    "det": "det", "equiv": "≡",
}
ACCENTS = {
    '"': "\u0308", "'": "\u0301", "`": "\u0300", "^": "\u0302", "~": "\u0303", "=": "\u0304",
    ".": "\u0307", "c": "\u0327", "v": "\u030c", "u": "\u0306", "H": "\u030b", "r": "\u030a",
    "k": "\u0328",
}
ACCENT = re.compile(r"\\([\"'`^~=.])\{?([A-Za-z])\}?|\\([cvuHrk])\{([A-Za-z])\}")
UNWRAP_TEXT = re.compile(r"\\(?:emph|textit|textbf|textrm|textsc|textsf|mbox|text)\{([^{}]*)\}")
UNWRAP_MATH = re.compile(r"\\(?:mathrm|mathit|mathbf|mathsf|operatorname|text|textrm)\{([^{}]*)\}")
FRAC = re.compile(r"\\[dt]?frac\{([^{}]*)\}\{([^{}]*)\}")
SQRT = re.compile(r"\\sqrt\{([^{}]*)\}")
SUPERSCRIPT = str.maketrans("0123456789+-−=()n", "⁰¹²³⁴⁵⁶⁷⁸⁹⁺⁻⁻⁼⁽⁾ⁿ")
SUBSCRIPT = str.maketrans("0123456789+-−=()", "₀₁₂₃₄₅₆₇₈₉₊₋₋₌₍₎")


def strip_comments(tex: str) -> str:
    return COMMENT.sub("", tex)


def abstract(paper: str) -> str:
    match = ABSTRACT.search(paper)
    if not match:
        raise ValueError("the paper has no abstract environment")
    return match.group(1)


LABEL = re.compile(r"\\label\{([^}]*)\}")


def appendix_labels(paper: str) -> set[str]:
    """Labels defined after \\appendix: where a reader finds a proof in full."""
    body = strip_comments(paper)
    at = body.find("\\appendix")
    return set(LABEL.findall(body[at:])) if at >= 0 else set()


def claim_keys(paper: str) -> list[str]:
    keys: list[str] = []
    for match in CLAIM.finditer(paper):
        keys += [key.strip() for key in match.group(1).split(",") if key.strip()]
    return keys


def cite_keys(paper: str) -> list[str]:
    keys: list[str] = []
    for match in CITE.finditer(strip_comments(paper)):
        keys += [key.strip() for key in match.group(1).split(",") if key.strip()]
    return keys


def val_keys(tex: str) -> list[str]:
    return VAL.findall(strip_comments(tex))


def values_tex(vals: dict[str, Value]) -> str:
    lines = [
        banner("results/ through the work's values module", "tex").rstrip("\n"),
        "\\makeatletter",
        "\\providecommand{\\val}[1]{\\@ifundefined{repokit@val@#1}"
        "{\\textbf{??#1??}\\PackageWarning{repokit}{undefined value #1}}"
        "{\\csname repokit@val@#1\\endcsname}}",
        "\\makeatother",
    ]
    for key in sorted(vals):
        lines.append(
            f"\\expandafter\\def\\csname repokit@val@{key}\\endcsname{{\\ensuremath{{{vals[key].tex()}}}}}"
        )
    return "\n".join(lines) + "\n"
