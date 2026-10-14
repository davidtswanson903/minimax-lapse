"""Every number a reader meets is printed one way, and every bound holds over what it covers."""
from fractions import Fraction

from repokit.values import BOUND_SUFFIX, bound, exact, measured, percent


def test_formatter_prints_each_kind():
    cases = [
        (measured("a", 6.453e-8, sig=3), "6.45 × 10⁻⁸", "6.45\\times10^{-8}"),
        (measured("b", 4.0012, sig=3), "4.00", "4.00"),
        (measured("c", 0.52741, sig=4), "0.5274", "0.5274"),
        (measured("d", -1234.5, sig=3), "−1.23 × 10³", "-1.23\\times10^{3}"),
        (measured("e", 85.2, sig=2, unit="s"), "85 s", "85\\,\\mathrm{s}"),
        (percent("f", 0.0159, sig=3), "1.59%", "1.59\\%"),
        (exact("g", 1, signed=True), "+1", "+1"),
        (exact("h", -11), "−11", "-11"),
        (exact("i", Fraction(124, 21)), "124/21", "124/21"),
        (bound("jBound", 1e-12, worst=2e-16), "10⁻¹²", "10^{-12}"),
        (bound("kBound", 5e-3, worst=0.0), "0.005", "0.005"),
    ]
    wrong = [(v.key, v.text(), v.tex()) for v, text, tex in cases if (v.text(), v.tex()) != (text, tex)]
    assert not wrong, wrong


def test_bound_names_match_kinds(values):
    wrong = [key for key, value in values.items() if (value.kind == "bound") != key.endswith(BOUND_SUFFIX)]
    assert not wrong, f"a bound's key ends in {BOUND_SUFFIX}, and only a bound's: {wrong}"


def test_bounds_hold(values):
    broken = [key for key, value in values.items() if value.kind == "bound" and not abs(value.worst) < value.x]
    assert not broken, broken
