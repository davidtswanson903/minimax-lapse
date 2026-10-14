import Saddle.Minimax

/-!
# The gain of reweighting

Corollary `cor:gain` of the appendix on uniform rings. Uniform weights give the Hamiltonian divided by the number of
places, `A_U = H / n`, whose lowest eigenvalue is `λ_min(H) / n`; it is one of the floors whose greatest is the
common value. So the gain of reweighting, `Γ = e* - λ_min(H)/n` (Eq. `eq:gain`), is never negative, and it
vanishes exactly when uniform weights are optimal.
-/

open Matrix
open scoped ComplexOrder

namespace Saddle

variable {n X : Type*} [Fintype n] [DecidableEq n] [Fintype X]

/-- The Hamiltonian, `H = ∑ₓ Dₓ`: the weighted sum with every weight one. -/
noncomputable def hamiltonian (D : X → Matrix n n ℂ) : Matrix n n ℂ :=
  weighted D fun _ => 1

omit [Fintype n] [DecidableEq n] in
/-- The Hamiltonian is the sum of the densities. -/
theorem hamiltonian_eq_sum (D : X → Matrix n n ℂ) : hamiltonian D = ∑ x, D x := by
  simp [hamiltonian, weighted]

/-- The uniform weighting, one over the number of places at each. -/
noncomputable def uniform (X : Type*) [Fintype X] : X → ℝ :=
  fun _ => (Fintype.card X : ℝ)⁻¹

/-- The uniform weighting is a weighting. -/
theorem uniform_mem_simplex [Nonempty X] : uniform X ∈ simplex X := by
  refine ⟨fun _ => inv_nonneg.mpr (Nat.cast_nonneg _), ?_⟩
  simp only [uniform, Finset.sum_const, Finset.card_univ, nsmul_eq_mul]
  exact mul_inv_cancel₀ (Nat.cast_ne_zero.mpr Fintype.card_ne_zero)

/-- The value a state gives the uniformly weighted sum is its value on the Hamiltonian over the number of places. -/
theorem expectation_uniform (ρ : Matrix n n ℂ) (D : X → Matrix n n ℂ) :
    expectation ρ (weighted D (uniform X)) = expectation ρ (hamiltonian D) / Fintype.card X := by
  rw [expectation_weighted, hamiltonian, expectation_weighted, div_eq_inv_mul, Finset.mul_sum]
  simp [uniform]

/-- Uniform weights certify the ground energy per place: `f(U) = λ_min(H) / n`. -/
theorem floor_uniform [Nonempty n] [Nonempty X] {D : X → Matrix n n ℂ} (hD : ∀ x, (D x).IsHermitian) :
    floor hD (uniform X) = lowest (weighted_isHermitian hD fun _ => 1) / Fintype.card X := by
  have hcard : (0 : ℝ) < Fintype.card X := Nat.cast_pos.mpr Fintype.card_pos
  obtain ⟨⟨σ, hσ, hσeq⟩, hσlow⟩ := lowest_isLeast (weighted_isHermitian hD fun _ => 1)
  obtain ⟨⟨τ, hτ, hτeq⟩, hτlow⟩ := lowest_isLeast (weighted_isHermitian hD (uniform X))
  have hσeq' : expectation σ (hamiltonian D) = lowest (weighted_isHermitian hD fun _ => 1) := hσeq
  have hτeq' : expectation τ (weighted D (uniform X)) = floor hD (uniform X) := hτeq
  apply le_antisymm
  · calc floor hD (uniform X) ≤ expectation σ (weighted D (uniform X)) := hτlow ⟨σ, hσ, rfl⟩
      _ = _ := by rw [expectation_uniform, hσeq']
  · calc lowest (weighted_isHermitian hD fun _ => 1) / Fintype.card X
        ≤ expectation τ (hamiltonian D) / Fintype.card X :=
          div_le_div_of_nonneg_right (hσlow ⟨τ, hτ, rfl⟩) hcard.le
      _ = floor hD (uniform X) := by rw [← expectation_uniform, hτeq']

/-- **Corollary `cor:gain` (the gain of reweighting).** For an optimal weighting `N`, whose floor is the common value
`e*`, the gain `Γ = e* - λ_min(H)/n` of Eq. `eq:gain` is never negative, and it vanishes exactly when uniform weights
are optimal too. -/
theorem saddle_gain [Nonempty n] [Nonempty X] {D : X → Matrix n n ℂ} (hD : ∀ x, (D x).IsHermitian)
    {N : X → ℝ} (hN : N ∈ simplex X) (hmax : IsMaxOn (floor hD) (simplex X) N) :
    0 ≤ floor hD N - lowest (weighted_isHermitian hD fun _ => 1) / Fintype.card X ∧
      (floor hD N - lowest (weighted_isHermitian hD fun _ => 1) / Fintype.card X = 0 ↔
        IsMaxOn (floor hD) (simplex X) (uniform X)) := by
  rw [← floor_uniform hD]
  have hle := isMaxOn_iff.mp hmax _ uniform_mem_simplex
  refine ⟨sub_nonneg.mpr hle, ⟨fun h => ?_, fun h => ?_⟩⟩
  · exact isMaxOn_iff.mpr fun M hM => (isMaxOn_iff.mp hmax M hM).trans (sub_eq_zero.mp h).le
  · exact sub_eq_zero.mpr (le_antisymm (isMaxOn_iff.mp h N hN) hle)

/-! ### Uniform rings -/

omit [DecidableEq n] in
/-- The value a state conjugated by `B` gives `A` is the value the state gives `Bᴴ A B`. -/
theorem expectation_conj (ρ A B : Matrix n n ℂ) :
    expectation (B * ρ * Bᴴ) A = expectation ρ (Bᴴ * A * B) := by
  unfold expectation
  rw [show B * ρ * Bᴴ * A = B * (ρ * Bᴴ * A) by simp only [Matrix.mul_assoc], trace_mul_comm]
  simp only [Matrix.mul_assoc]

/-- A state conjugated by a unitary is a state. -/
theorem conj_mem_states {ρ B : Matrix n n ℂ} (hρ : ρ ∈ states n) (hB : Bᴴ * B = 1) : B * ρ * Bᴴ ∈ states n := by
  refine ⟨hρ.1.mul_mul_conjTranspose_same B, ?_⟩
  rw [Matrix.mul_assoc, trace_mul_comm, Matrix.mul_assoc, hB, Matrix.mul_one, hρ.2]

/-- Undoing the rotation `k` times carries each density back `k` places: `(R†)ᵏ Dₓ Rᵏ = D_{x-k}`. -/
theorem rotate_back {m : ℕ} {D : ZMod m → Matrix n n ℂ} {R : Matrix n n ℂ} (hR : Rᴴ * R = 1)
    (hrot : ∀ x, R * D x * Rᴴ = D (x + 1)) (x : ZMod m) (k : ℕ) :
    Rᴴ ^ k * D x * R ^ k = D (x - k) := by
  induction k generalizing x with
  | zero => simp
  | succ k ih =>
    have hback : Rᴴ * D x * R = D (x - 1) := by
      conv_lhs => rw [← sub_add_cancel x 1, ← hrot (x - 1)]
      simp only [← Matrix.mul_assoc, hR, Matrix.one_mul]
      rw [Matrix.mul_assoc, hR, Matrix.mul_one]
    rw [pow_succ, pow_succ', show Rᴴ ^ k * Rᴴ * D x * (R * R ^ k) = Rᴴ ^ k * (Rᴴ * D x * R) * R ^ k by
      simp only [Matrix.mul_assoc], hback, ih]
    congr 1
    push_cast
    ring

/-- **Proposition `prop:uniform`.** On a ring of `m` cells whose rotation, a unitary `R`, carries each density to the
next place's, `R Dₓ R† = D_{x+1}`, uniform weights are optimal and the common value is the ground energy per cell,
`e* = λ_min(H) / m` (Eq. `eq:uniform`). The proof is the paper's: a ground state averaged over one turn of the ring
holds every place at `λ_min(H) / m`, which no weighting's floor can exceed, and uniform weights reach it. -/
theorem uniform_exact [Nonempty n] {m : ℕ} [NeZero m] {D : ZMod m → Matrix n n ℂ} (hD : ∀ x, (D x).IsHermitian)
    {R : Matrix n n ℂ} (hR : Rᴴ * R = 1) (hrot : ∀ x, R * D x * Rᴴ = D (x + 1)) :
    IsMaxOn (floor hD) (simplex (ZMod m)) (uniform (ZMod m)) ∧
      floor hD (uniform (ZMod m)) = lowest (weighted_isHermitian hD fun _ => 1) / m := by
  have hm : (Fintype.card (ZMod m) : ℝ) = m := by rw [ZMod.card]
  have hfloor : floor hD (uniform (ZMod m)) = lowest (weighted_isHermitian hD fun _ => 1) / m := by
    rw [floor_uniform, hm]
  refine ⟨isMaxOn_iff.mpr fun N hN => ?_, hfloor⟩
  -- a ground state of `H`, and its average over one turn of the ring
  obtain ⟨⟨σ, hσ, hσeq⟩, -⟩ := lowest_isLeast (weighted_isHermitian hD fun _ => 1)
  have hσeq' : expectation σ (hamiltonian D) = lowest (weighted_isHermitian hD fun _ => 1) := hσeq
  have hpow : ∀ k : ℕ, (R ^ k)ᴴ * R ^ k = 1 := fun k => by
    rw [conjTranspose_pow]
    induction k with
    | zero => simp
    | succ k ih => rw [pow_succ', pow_succ, Matrix.mul_assoc, ← Matrix.mul_assoc (Rᴴ ^ k), ih, Matrix.one_mul, hR]
  set turns : ZMod m → Matrix n n ℂ := fun k => R ^ k.val * σ * (R ^ k.val)ᴴ
  set average := ∑ k, uniform (ZMod m) k • turns k
  have haverage : average ∈ states n :=
    states_convex.sum_mem (fun k _ => uniform_mem_simplex.1 k) uniform_mem_simplex.2
      (fun k _ => conj_mem_states hσ (hpow _))
  -- every place's density in the average is the ground energy per cell
  have hevery : ∀ x, expectation average (D x) = lowest (weighted_isHermitian hD fun _ => 1) / m := fun x => by
    calc expectation average (D x) = ∑ k, uniform (ZMod m) k * expectation (turns k) (D x) := by
          simp only [average, expectation_sum_left, expectation_smul]
      _ = ∑ k : ZMod m, uniform (ZMod m) k * expectation σ (D (x - k)) := by
          refine Finset.sum_congr rfl fun k _ => ?_
          simp only [turns]
          rw [expectation_conj, conjTranspose_pow, rotate_back hR hrot, ZMod.natCast_zmod_val]
      _ = (m : ℝ)⁻¹ * ∑ y, expectation σ (D y) := by
          simp only [uniform, hm, ← Finset.mul_sum]
          congr 1
          exact Equiv.sum_comp (Equiv.subLeft x) (fun y => expectation σ (D y))
      _ = lowest (weighted_isHermitian hD fun _ => 1) / m := by
          rw [← expectation_sum_right, ← hamiltonian_eq_sum, hσeq', inv_mul_eq_div]
  calc floor hD N ≤ peak D average := floor_le_peak hD haverage hN
    _ = floor hD (uniform (ZMod m)) := by
        rw [hfloor]
        exact le_antisymm (ciSup_le fun x => (hevery x).le) ((hevery 0).symm.le.trans (expectation_le_peak D average 0))

end Saddle
