import Saddle.Minimax

/-!
# The bracket and the certificate

Corollary `cor:gap` of the appendix on the saddle. Every weighting's floor lies at or below the common value and
every state's peak at or above it, so a computed state and weighting whose peak and floor agree within a margin are
each within that margin of the common value. A state that holds every density at or below a weighting's floor
attains the left side of Eq. `eq:saddle`, and the weighting the right: the conditions of Eq. `eq:conditions`, with
the common value the lowest eigenvalue of the weighted sum, certify a saddle pair.
-/

open Matrix
open scoped ComplexOrder

namespace Saddle

variable {n X : Type*} [Fintype n] [DecidableEq n] [Fintype X]

/-- **Corollary `cor:gap` (the bracket).** The common value `e*`, attained by a state's peak and a weighting's
floor, lies at or above every weighting's floor and at or below every state's peak. -/
theorem saddle_bracket [Nonempty n] [Nonempty X] {D : X → Matrix n n ℂ} (hD : ∀ x, (D x).IsHermitian) :
    ∃ e : ℝ, (∃ ρ ∈ states n, peak D ρ = e) ∧ (∃ N ∈ simplex X, floor hD N = e) ∧
      (∀ N ∈ simplex X, floor hD N ≤ e) ∧ (∀ ρ ∈ states n, e ≤ peak D ρ) := by
  obtain ⟨ρ₀, hρ₀, N₀, hN₀, hmin, hmax, h₀⟩ := saddle_value hD
  exact ⟨peak D ρ₀, ⟨ρ₀, hρ₀, rfl⟩, ⟨N₀, hN₀, h₀.symm⟩,
    fun N hN => (isMaxOn_iff.mp hmax N hN).trans h₀.ge, fun ρ hρ => isMinOn_iff.mp hmin ρ hρ⟩

/-- **Corollary `cor:gap` (the certificate).** A state that holds every density at or below a weighting's floor
attains the least peak over all states, and the weighting the greatest floor over all weightings. -/
theorem saddle_certificate [Nonempty n] [Nonempty X] {D : X → Matrix n n ℂ} (hD : ∀ x, (D x).IsHermitian)
    {ρ : Matrix n n ℂ} {N : X → ℝ} (hρ : ρ ∈ states n) (hN : N ∈ simplex X)
    (hle : ∀ x, expectation ρ (D x) ≤ floor hD N) :
    IsMinOn (peak D) (states n) ρ ∧ IsMaxOn (floor hD) (simplex X) N := by
  have hpeak : peak D ρ ≤ floor hD N := ciSup_le hle
  exact ⟨isMinOn_iff.mpr fun σ hσ => hpeak.trans (floor_le_peak hD hσ hN),
    isMaxOn_iff.mpr fun M hM => (floor_le_peak hD hρ hM).trans hpeak⟩

end Saddle
