import Saddle.Minimax

/-!
# The conditions at the saddle

Proposition `prop:conditions` of the appendix on the saddle: any state attaining the left side of Eq. `eq:saddle`
and any weighting attaining the right satisfy the four conditions of Eq. `eq:conditions`. The proof is the paper's
chain `e* = f(N*) ≤ Tr(ρ* A_{N*}) = ∑ₓ N*ₓ Tr(ρ* Dₓ) ≤ e*`, whose equalities give the conditions.
-/

open Matrix
open scoped ComplexOrder

namespace Saddle

variable {n X : Type*} [Fintype n] [DecidableEq n] [Fintype X]

/-- **Proposition `prop:conditions` (the conditions).** If the state `ρ` attains the left side of
Eq. `eq:saddle` and the weighting `N` the right, the four conditions of Eq. `eq:conditions` hold at the
common value `e* = peak D ρ`: the weighted sum holds the state at `e*`, no density exceeds `e*`, weight sits
only where a density reaches `e*`, and the weights sum to one. -/
theorem saddle_conditions [Nonempty n] [Nonempty X] {D : X → Matrix n n ℂ} (hD : ∀ x, (D x).IsHermitian)
    {ρ : Matrix n n ℂ} {N : X → ℝ} (hρ : ρ ∈ states n) (hN : N ∈ simplex X)
    (hmin : IsMinOn (peak D) (states n) ρ) (hmax : IsMaxOn (floor hD) (simplex X) N) :
    weighted D N * ρ = (peak D ρ : ℂ) • ρ ∧
      (∀ x, expectation ρ (D x) ≤ peak D ρ) ∧
      (∀ x, N x * (expectation ρ (D x) - peak D ρ) = 0) ∧
      ∑ x, N x = 1 := by
  -- the least peak is the greatest floor, so the pair's two values are one, `e*`
  obtain ⟨ρ₀, hρ₀, N₀, hN₀, -, -, h₀⟩ := saddle_value hD
  have hvalue : peak D ρ = floor hD N := by
    have := isMinOn_iff.mp hmin ρ₀ hρ₀
    have := isMaxOn_iff.mp hmax N₀ hN₀
    have := floor_le_peak hD hρ hN
    linarith
  have hle := expectation_le_peak D ρ
  -- `e* = f(N) ≤ Tr(ρ A_N) = ∑ₓ Nₓ Tr(ρ Dₓ) ≤ e*`, so both are equalities
  have hlow : floor hD N ≤ expectation ρ (weighted D N) := (lowest_isLeast _).2 ⟨ρ, hρ, rfl⟩
  have hhigh : ∑ x, N x * expectation ρ (D x) ≤ peak D ρ := by
    calc ∑ x, N x * expectation ρ (D x) ≤ ∑ x, N x * peak D ρ :=
          Finset.sum_le_sum fun x _ => mul_le_mul_of_nonneg_left (hle x) (hN.1 x)
      _ = peak D ρ := by rw [← Finset.sum_mul, hN.2, one_mul]
  rw [expectation_weighted] at hlow
  refine ⟨?_, hle, ?_, hN.2⟩
  · rw [hvalue]
    exact (expectation_eq_lowest_iff (weighted_isHermitian hD N) hρ).mp
      (by rw [expectation_weighted]; exact le_antisymm (hhigh.trans hvalue.le) hlow)
  · have hzero : ∑ x, N x * (peak D ρ - expectation ρ (D x)) = 0 := by
      simp only [mul_sub, Finset.sum_sub_distrib, ← Finset.sum_mul, hN.2, one_mul]
      linarith
    have hterm := (Finset.sum_eq_zero_iff_of_nonneg fun x _ =>
      mul_nonneg (hN.1 x) (sub_nonneg.mpr (hle x))).mp hzero
    intro x
    have := hterm x (Finset.mem_univ x)
    linarith [this, mul_sub (N x) (peak D ρ) (expectation ρ (D x)), mul_sub (N x) (expectation ρ (D x)) (peak D ρ)]

end Saddle
