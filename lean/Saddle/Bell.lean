import Mathlib.Algebra.Star.CHSH

/-!
# Bell's assumptions, factored

Proposition `prop:factorization` of the appendix on Bell's assumptions. For readings self-adjoint with square one, each
of one part's commuting with each of the other's, the CHSH combination is at most `2 √2`; if all four commute, at most
`2`. Both are Mathlib's (`tsirelson_inequality`, `CHSH_inequality_of_comm`), for ordered star-algebras over the reals;
stated here in the paper's terms. The gap between the bounds is a part's own readings commuting.
-/

namespace Saddle

variable {R : Type*}

/-- **Proposition `prop:factorization`, the bound of one description.** When all four readings commute, so that
one description holds both of a part's settings, the CHSH combination is at most `2`. -/
theorem chsh_of_one_description [CommRing R] [PartialOrder R] [StarRing R] [StarOrderedRing R] [Algebra ℝ R]
    [IsOrderedModule ℝ R] (A₀ A₁ B₀ B₁ : R) (T : IsCHSHTuple A₀ A₁ B₀ B₁) :
    A₀ * B₀ + A₀ * B₁ + A₁ * B₀ - A₁ * B₁ ≤ 2 :=
  CHSH_inequality_of_comm A₀ A₁ B₀ B₁ T

/-- **Proposition `prop:factorization`, the bound of the parts.** When only readings of different parts commute,
the CHSH combination is at most `2 √2`, written `√2 ^ 3`. -/
theorem chsh_of_parts [Ring R] [PartialOrder R] [StarRing R] [StarOrderedRing R] [Algebra ℝ R]
    [IsOrderedModule ℝ R] [StarModule ℝ R] (A₀ A₁ B₀ B₁ : R) (T : IsCHSHTuple A₀ A₁ B₀ B₁) :
    A₀ * B₀ + A₀ * B₁ + A₁ * B₀ - A₁ * B₁ ≤ √2 ^ 3 • (1 : R) :=
  tsirelson_inequality A₀ A₁ B₀ B₁ T

end Saddle
