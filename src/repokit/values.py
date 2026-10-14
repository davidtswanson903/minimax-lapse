"""Named values and the one formatter that prints them.

Every number a reader meets is a Value: computed from stored results, given a
kind, and printed here for text and for TeX alike, so the paper and the
Markdown can never print one value two ways.

Kinds:
  measured  a computed or measured number, printed to its significant figures
  input     a parameter taken from a source, printed the same way
  percent   a fraction, printed as a percentage to its significant figures
  exact     an integer or a Fraction, printed exactly, optionally signed
  bound     a chosen threshold with the worst value it covers; printed as chosen
"""
from __future__ import annotations

import math
from dataclasses import dataclass
from fractions import Fraction

KINDS = ("measured", "input", "percent", "exact", "bound")
BOUND_SUFFIX = "Bound"
_SUP = str.maketrans("0123456789-+", "⁰¹²³⁴⁵⁶⁷⁸⁹⁻⁺")
_MINUS = "\u2212"


@dataclass(frozen=True)
class Value:
    key: str
    kind: str
    x: float | int | Fraction
    sig: int | None = None
    worst: float | None = None
    unit: str = ""
    note: str = ""
    signed: bool = False

    def __post_init__(self) -> None:
        if self.kind not in KINDS:
            raise ValueError(f"{self.key}: unknown kind {self.kind!r}")

    def text(self) -> str:
        return _with_unit(_number(self, tex=False), self.unit, tex=False)

    def tex(self) -> str:
        return _with_unit(_number(self, tex=True), self.unit, tex=True)


def measured(key: str, x: float, sig: int = 3, unit: str = "", note: str = "") -> Value:
    return Value(key, "measured", float(x), sig=sig, unit=unit, note=note)


def inp(key: str, x: float, sig: int = 3, unit: str = "", note: str = "") -> Value:
    return Value(key, "input", float(x), sig=sig, unit=unit, note=note)


def percent(key: str, fraction: float, sig: int = 3, note: str = "") -> Value:
    return Value(key, "percent", float(fraction), sig=sig, note=note)


def exact(key: str, x: int | Fraction, signed: bool = False, unit: str = "", note: str = "") -> Value:
    if isinstance(x, float):
        if not x.is_integer():
            raise ValueError(f"{key}: an exact value is an integer or a Fraction")
        x = int(x)
    return Value(key, "exact", x, signed=signed, unit=unit, note=note)


def bound(key: str, x: float, worst: float, unit: str = "", note: str = "") -> Value:
    """A chosen bound and the worst measured value it covers. The build fails if
    the measured value reaches the bound: that is a regression, not noise."""
    x = float(x)
    if not x > 0:
        raise ValueError(f"{key}: a bound is positive")
    if not abs(float(worst)) < x:
        raise ValueError(
            f"{key}: the measured value {float(worst):.3e} is not below the bound {x:.0e}; "
            "investigate before loosening the bound"
        )
    return Value(key, "bound", x, worst=float(worst), unit=unit, note=note)


def _number(v: Value, tex: bool) -> str:
    if v.kind == "exact":
        return _exact(v.x, v.signed, tex)
    if v.kind == "bound":
        return _bound(float(v.x), tex)
    if v.kind == "percent":
        return _significant(100 * float(v.x), v.sig or 3, tex) + ("\\%" if tex else "%")
    return _significant(float(v.x), v.sig or 3, tex)


def _exponent(a: float, sig: int) -> int:
    rounded = float(f"{a:.{sig - 1}e}")
    return math.floor(math.log10(rounded) + 1e-12)


def _power(mantissa: str, e: int, tex: bool) -> str:
    if tex:
        return f"10^{{{e}}}" if mantissa == "1" else f"{mantissa}\\times10^{{{e}}}"
    sup = str(e).translate(_SUP)
    return f"10{sup}" if mantissa == "1" else f"{mantissa} × 10{sup}"


def _signed(negative: bool, body: str, tex: bool) -> str:
    return (("-" if tex else _MINUS) + body) if negative else body


def _significant(x: float, sig: int, tex: bool) -> str:
    if x == 0:
        return "0"
    a = abs(x)
    e = _exponent(a, sig)
    if -3 <= e < sig:
        body = f"{a:.{max(sig - 1 - e, 0)}f}"
    else:
        mantissa = f"{a / 10 ** e:.{sig - 1}f}"
        if float(mantissa) >= 10:
            e += 1
            mantissa = f"{a / 10 ** e:.{sig - 1}f}"
        body = _power(mantissa, e, tex)
    return _signed(x < 0, body, tex)


def _bound(x: float, tex: bool) -> str:
    e = math.floor(math.log10(x) + 1e-12)
    if -3 <= e < 4:
        return f"{x:.12g}"
    return _power(f"{x / 10 ** e:.12g}", e, tex)


def _exact(x: int | Fraction, signed: bool, tex: bool) -> str:
    if isinstance(x, Fraction) and x.denominator != 1:
        body, negative = f"{abs(x.numerator)}/{x.denominator}", x < 0
    else:
        body, negative = str(abs(int(x))), int(x) < 0
    if negative:
        return _signed(True, body, tex)
    return ("+" + body) if signed and x != 0 else body


def _with_unit(number: str, unit: str, tex: bool) -> str:
    if not unit:
        return number
    return f"{number}\\,\\mathrm{{{unit}}}" if tex else f"{number} {unit}"
