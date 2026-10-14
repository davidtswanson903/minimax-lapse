import Saddle.Minimax

/-!
# Every weighting is some division's

Proposition `prop:division` of the appendix on rates. Shift each place's density by its own constant,
`D_x → D_x + c_x I`, with the constants summing to zero: the densities still sum to the same Hamiltonian, and every
floor moves by `∑ₓ Nₓ cₓ`. Given any weighting `N₀`, let `ψ` be a state at the floor of `A_{N₀}`, `gₓ` the value it
gives each density, and `cₓ = ḡ - gₓ` with `ḡ` their mean. Every shifted floor is at most `ḡ`, since `ψ` gives the
shifted weighted sum `∑ₓ Nₓ (gₓ + cₓ) = ḡ`; and the shifted floor at `N₀` is at least `ḡ`, since any state gives
`A_{N₀}` at least `∑ₓ N₀ₓ gₓ`. So `N₀` is an optimal weighting of the shifted densities.
-/

open Matrix

namespace Saddle

variable {n X : Type*} [Fintype n] [DecidableEq n] [Fintype X]

/-- The densities with each place's shifted by its own constant, `D_x + c_x I`. -/
noncomputable def divided (D : X → Matrix n n ℂ) (c : X → ℝ) : X → Matrix n n ℂ :=
  fun x => D x + (c x : ℂ) • 1

omit [Fintype n] [Fintype X] in
/-- The shifted densities are Hermitian. -/
theorem divided_isHermitian {D : X → Matrix n n ℂ} (hD : ∀ x, (D x).IsHermitian) (c : X → ℝ) :
    ∀ x, (divided D c x).IsHermitian := fun x => by
  unfold Matrix.IsHermitian divided
  rw [conjTranspose_add, (hD x).eq, conjTranspose_smul, conjTranspose_one, Complex.star_def, Complex.conj_ofReal]

omit [Fintype n] in
/-- Constants that sum to zero leave the Hamiltonian, the densities' sum, unchanged. -/
theorem sum_divided (D : X → Matrix n n ℂ) {c : X → ℝ} (hc : ∑ x, c x = 0) :
    ∑ x, divided D c x = ∑ x, D x := by
  simp only [divided, Finset.sum_add_distrib, ← Finset.sum_smul, ← Complex.ofReal_sum, hc, Complex.ofReal_zero,
    zero_smul, add_zero]

omit [Fintype X] in
/-- A state gives a shifted density its value of the density plus the constant. -/
theorem expectation_divided {ρ : Matrix n n ℂ} (hρ : ρ ∈ states n) (D : X → Matrix n n ℂ) (c : X → ℝ) (x : X) :
    expectation ρ (divided D c x) = expectation ρ (D x) + c x := by
  simp only [divided, expectation, Matrix.mul_add, trace_add, Complex.add_re, Matrix.mul_smul, trace_smul, hρ.2,
    smul_eq_mul, mul_one, Complex.ofReal_re]

/-- A state gives the shifted weighted sum its value of the weighted sum plus `∑ₓ Nₓ cₓ`. -/
theorem expectation_weighted_divided {ρ : Matrix n n ℂ} (hρ : ρ ∈ states n) (D : X → Matrix n n ℂ) (c : X → ℝ)
    (N : X → ℝ) : expectation ρ (weighted (divided D c) N) = expectation ρ (weighted D N) + ∑ x, N x * c x := by
  simp only [expectation_weighted, expectation_divided hρ, mul_add, Finset.sum_add_distrib]

/-- **Proposition `prop:division`.** For every weighting `N₀` there are constants `c` summing to zero such that the
densities `D_x + c_x I`, which sum to the same Hamiltonian, have `N₀` as an optimal weighting. -/
theorem every_weighting_some_division [Nonempty n] {D : X → Matrix n n ℂ} (hD : ∀ x, (D x).IsHermitian)
    {N₀ : X → ℝ} (hN₀ : N₀ ∈ simplex X) :
    ∃ c : X → ℝ, ∑ x, c x = 0 ∧ ∑ x, divided D c x = ∑ x, D x ∧
      IsMaxOn (floor (divided_isHermitian hD c)) (simplex X) N₀ := by
  have hX : Nonempty X := by
    by_contra h
    rw [not_nonempty_iff] at h
    have := hN₀.2
    simp at this
  have hcard : (Fintype.card X : ℝ) ≠ 0 := Nat.cast_ne_zero.mpr Fintype.card_ne_zero
  obtain ⟨⟨ψ, hψ, hψeq⟩, -⟩ := lowest_isLeast (weighted_isHermitian hD N₀)
  simp only at hψeq
  set g : X → ℝ := fun x => expectation ψ (D x) with hg
  set m : ℝ := (∑ x, g x) / Fintype.card X with hm
  have hc : ∑ x, (m - g x) = 0 := by
    rw [Finset.sum_sub_distrib, Finset.sum_const, Finset.card_univ, nsmul_eq_mul, hm, mul_div_cancel₀ _ hcard,
      sub_self]
  have hsum : ∀ N ∈ simplex X, ∑ x, N x * (g x + (m - g x)) = m := fun N hN => by
    simp only [add_sub_cancel, ← Finset.sum_mul, hN.2, one_mul]
  refine ⟨fun x => m - g x, hc, sum_divided D hc, isMaxOn_iff.mpr fun N hN => ?_⟩
  -- every shifted floor is at most `m`: the state `ψ` gives the shifted weighted sum exactly `m`
  have hupper : floor (divided_isHermitian hD fun x => m - g x) N ≤ m := by
    have h := (lowest_isLeast (weighted_isHermitian (divided_isHermitian hD fun x => m - g x) N)).2 ⟨ψ, hψ, rfl⟩
    simp only at h
    rw [floor]
    calc _ ≤ _ := h
      _ = m := by
        rw [expectation_weighted]
        simp only [expectation_divided hψ]
        exact hsum N hN
  -- the shifted floor at `N₀` is at least `m`: no state gives `A_{N₀}` less than `ψ` does
  have hlower : m ≤ floor (divided_isHermitian hD fun x => m - g x) N₀ := by
    obtain ⟨⟨ρ, hρ, hρeq⟩, -⟩ :=
      lowest_isLeast (weighted_isHermitian (divided_isHermitian hD fun x => m - g x) N₀)
    simp only at hρeq
    have hge := (lowest_isLeast (weighted_isHermitian hD N₀)).2 ⟨ρ, hρ, rfl⟩
    simp only at hge
    rw [floor, ← hρeq, expectation_weighted_divided hρ]
    have hψN : expectation ψ (weighted D N₀) = ∑ x, N₀ x * g x := expectation_weighted ψ D N₀
    have := hsum N₀ hN₀
    simp only [mul_add, Finset.sum_add_distrib] at this
    linarith
  linarith
