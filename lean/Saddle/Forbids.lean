import Saddle.Conditions

/-!
# What the saddle forbids

Corollary `cor:forbids` of the appendix on the saddle: a place some optimal state leaves below the common value
carries no weight in any optimal weighting, and a place with weight sits at the common value in every optimal
state. It is the third condition of Eq. `eq:conditions`, read for any optimal pair.
-/

open Matrix

namespace Saddle

variable {n X : Type*} [Fintype n] [DecidableEq n] [Fintype X]

/-- **Corollary `cor:forbids` (what is forbidden).** For any state `ρ` attaining the left side of Eq. `eq:saddle`
and any weighting `N` attaining the right: if `ρ` leaves place `x` below the common value, `N` gives it no weight;
and if `N` gives `x` weight, `ρ` holds it exactly at the common value. -/
theorem saddle_forbids [Nonempty n] [Nonempty X] {D : X → Matrix n n ℂ} (hD : ∀ x, (D x).IsHermitian)
    {ρ : Matrix n n ℂ} {N : X → ℝ} (hρ : ρ ∈ states n) (hN : N ∈ simplex X)
    (hmin : IsMinOn (peak D) (states n) ρ) (hmax : IsMaxOn (floor hD) (simplex X) N) (x : X) :
    (expectation ρ (D x) < peak D ρ → N x = 0) ∧ (0 < N x → expectation ρ (D x) = peak D ρ) := by
  have hslack := (saddle_conditions hD hρ hN hmin hmax).2.2.1 x
  refine ⟨fun hlow => ?_, fun hpos => ?_⟩
  · exact (mul_eq_zero.mp hslack).resolve_right (sub_ne_zero.mpr hlow.ne)
  · exact sub_eq_zero.mp ((mul_eq_zero.mp hslack).resolve_left hpos.ne')

end Saddle
