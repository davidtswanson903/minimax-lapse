import Mathlib.Analysis.Matrix.Order
import Mathlib.Analysis.Matrix.PosDef
import Mathlib.Analysis.Matrix.Normed
import Mathlib.Analysis.RCLike.Lemmas

/-!
# States, densities and the weighted sum

The objects of the appendix on the saddle (`app:saddle`): the density matrices on a finite-dimensional
space, the value a state gives an operator, the weighted sum `A_N = ∑ₓ Nₓ Dₓ` of the densities, its lowest
eigenvalue, the peak of a state and the floor a weighting certifies. The set of states is compact and
convex, which is what lets both sides of the saddle be attained.
-/

open Matrix
open scoped ComplexOrder

namespace Saddle

/-- The states: density matrices on the space indexed by `n`, positive semidefinite with trace one. -/
def states (n : Type*) [Fintype n] [DecidableEq n] : Set (Matrix n n ℂ) :=
  {ρ | ρ.PosSemidef ∧ ρ.trace = 1}

/-- The weightings: nonnegative weights on the places that sum to one, the simplex `Δ`. -/
def simplex (X : Type*) [Fintype X] : Set (X → ℝ) :=
  {N | (∀ x, 0 ≤ N x) ∧ ∑ x, N x = 1}

variable {n X : Type*} [Fintype n] [DecidableEq n] [Fintype X]

/-- The value the state `ρ` gives the operator `A`: the real part of `Tr(ρ A)`. -/
noncomputable def expectation (ρ A : Matrix n n ℂ) : ℝ :=
  (ρ * A).trace.re

/-- The weighted sum `A_N = ∑ₓ Nₓ Dₓ` of the densities. -/
noncomputable def weighted (D : X → Matrix n n ℂ) (N : X → ℝ) : Matrix n n ℂ :=
  ∑ x, (N x : ℂ) • D x

/-- The lowest eigenvalue `λ_min(A)` of a Hermitian matrix. -/
noncomputable def lowest {A : Matrix n n ℂ} (hA : A.IsHermitian) : ℝ :=
  ⨅ i, hA.eigenvalues i

omit [Fintype n] [DecidableEq n] in
/-- A weighted sum of Hermitian densities with real weights is Hermitian. -/
theorem weighted_isHermitian {D : X → Matrix n n ℂ} (hD : ∀ x, (D x).IsHermitian) (N : X → ℝ) :
    (weighted D N).IsHermitian := by
  unfold Matrix.IsHermitian weighted
  rw [conjTranspose_sum]
  refine Finset.sum_congr rfl fun x _ => ?_
  rw [conjTranspose_smul, (hD x).eq, Complex.star_def, Complex.conj_ofReal]

/-- The peak of a state: its greatest density, `maxₓ Tr(ρ Dₓ)`. -/
noncomputable def peak (D : X → Matrix n n ℂ) (ρ : Matrix n n ℂ) : ℝ :=
  ⨆ x, expectation ρ (D x)

/-- The floor a weighting certifies: `f(N) = λ_min(A_N)`. -/
noncomputable def floor {D : X → Matrix n n ℂ} (hD : ∀ x, (D x).IsHermitian) (N : X → ℝ) : ℝ :=
  lowest (weighted_isHermitian hD N)

/-! ### The value a state gives an operator -/

/-- The value a state gives a weighted sum is the weighted sum of the values it gives the densities. -/
theorem expectation_weighted (ρ : Matrix n n ℂ) (D : X → Matrix n n ℂ) (N : X → ℝ) :
    expectation ρ (weighted D N) = ∑ x, N x * expectation ρ (D x) := by
  simp only [expectation, weighted, Matrix.mul_sum, trace_sum, Complex.re_sum, mul_smul_comm, trace_smul,
    smul_eq_mul, Complex.re_ofReal_mul]

omit [DecidableEq n] in
/-- The value a state gives an operator is additive in the state. -/
theorem expectation_add (ρ σ A : Matrix n n ℂ) :
    expectation (ρ + σ) A = expectation ρ A + expectation σ A := by
  simp only [expectation, Matrix.add_mul, trace_add, Complex.add_re]

omit [DecidableEq n] in
/-- The value a state gives an operator is homogeneous in the state, for real factors. -/
theorem expectation_smul (c : ℝ) (ρ A : Matrix n n ℂ) :
    expectation (c • ρ) A = c * expectation ρ A := by
  simp only [expectation, Matrix.smul_mul, trace_smul, Complex.real_smul, Complex.re_ofReal_mul]

omit [DecidableEq n] in
/-- The value a sum of states gives an operator is the sum of the values. -/
theorem expectation_sum_left {ι : Type*} (s : Finset ι) (ρ : ι → Matrix n n ℂ) (A : Matrix n n ℂ) :
    expectation (∑ i ∈ s, ρ i) A = ∑ i ∈ s, expectation (ρ i) A := by
  simp only [expectation, Finset.sum_mul, trace_sum, Complex.re_sum]

omit [DecidableEq n] in
/-- The value a state gives a sum of operators is the sum of the values. -/
theorem expectation_sum_right {ι : Type*} (s : Finset ι) (ρ : Matrix n n ℂ) (A : ι → Matrix n n ℂ) :
    expectation ρ (∑ i ∈ s, A i) = ∑ i ∈ s, expectation ρ (A i) := by
  simp only [expectation, Matrix.mul_sum, trace_sum, Complex.re_sum]

omit [DecidableEq n] in
/-- The value a state gives an operator is continuous in the state. -/
theorem continuous_expectation (A : Matrix n n ℂ) : Continuous fun ρ : Matrix n n ℂ => expectation ρ A := by
  unfold expectation
  fun_prop

/-! ### The set of states is compact and convex -/

/-- The states form a closed set. -/
theorem states_isClosed : IsClosed (states n) :=
  posSemidef_is_closed.inter (isClosed_eq (by fun_prop) continuous_const)

/-- The states form a convex set. -/
theorem states_convex : Convex ℝ (states n) := by
  rintro ρ ⟨hρ, hρ1⟩ σ ⟨hσ, hσ1⟩ a b ha hb hab
  refine ⟨by simpa using (hρ.smul ha).add (hσ.smul hb), ?_⟩
  rw [trace_add, trace_smul, trace_smul, hρ1, hσ1, Complex.real_smul, Complex.real_smul, mul_one, mul_one]
  exact_mod_cast hab

/-- There is a state: the maximally mixed one, the identity divided by the dimension. -/
theorem states_nonempty [Nonempty n] : (states n).Nonempty := by
  have hcard : (Fintype.card n : ℂ) ≠ 0 := Nat.cast_ne_zero.mpr Fintype.card_ne_zero
  refine ⟨(Fintype.card n : ℝ)⁻¹ • (1 : Matrix n n ℂ), PosSemidef.one.smul (by positivity), ?_⟩
  rw [trace_smul, trace_one, Complex.real_smul]
  push_cast
  exact inv_mul_cancel₀ hcard

/-- Every entry of a state is at most one in size. The two-by-two principal minor on rows `i` and `j` of a
positive semidefinite matrix has a nonnegative determinant, so `|ρᵢⱼ|² ≤ ρᵢᵢ ρⱼⱼ`, and each diagonal entry
lies between zero and the trace, which is one. -/
theorem norm_entry_le_one {ρ : Matrix n n ℂ} (hρ : ρ ∈ states n) (i j : n) : ‖ρ i j‖ ≤ 1 := by
  obtain ⟨hpsd, htr⟩ := hρ
  have hdiag : ∀ k, 0 ≤ (ρ k k).re ∧ (ρ k k).im = 0 := fun k =>
    ⟨(Complex.le_def.mp hpsd.diag_nonneg).1, (Complex.le_def.mp hpsd.diag_nonneg).2.symm⟩
  have hle : ∀ k, (ρ k k).re ≤ 1 := fun k => by
    have hsum : ∑ l, (ρ l l).re = 1 := by simpa [Matrix.trace, Complex.re_sum] using congrArg Complex.re htr
    exact hsum ▸ Finset.single_le_sum (f := fun l => (ρ l l).re) (fun l _ => (hdiag l).1) (Finset.mem_univ k)
  -- the minor on rows `i` and `j`
  have hminor := (hpsd.submatrix ![i, j]).det_nonneg
  rw [det_fin_two] at hminor
  simp only [submatrix_apply, Matrix.cons_val_zero, Matrix.cons_val_one] at hminor
  have hji : ρ j i = star (ρ i j) := by rw [← hpsd.isHermitian.apply i j, star_star]
  rw [hji, ← starRingEnd_apply, Complex.mul_conj] at hminor
  have hre := (Complex.le_def.mp hminor).1
  simp only [Complex.zero_re, Complex.sub_re, Complex.mul_re, (hdiag i).2, (hdiag j).2, mul_zero, sub_zero,
    Complex.ofReal_re] at hre
  have hsq : Complex.normSq (ρ i j) ≤ 1 := by
    nlinarith [hdiag i, hdiag j, hle i, hle j, mul_le_mul (hle i) (hle j) (hdiag j).1 zero_le_one]
  rw [← sq_le_one_iff₀ (norm_nonneg (ρ i j)), Complex.sq_norm]
  exact hsq

open scoped Matrix.Norms.Elementwise in
/-- The states form a compact set: closed, and bounded because every entry is at most one in size. -/
theorem states_isCompact : IsCompact (states n) := by
  have := FiniteDimensional.proper_rclike ℂ (Matrix n n ℂ)
  refine Metric.isCompact_of_isClosed_isBounded states_isClosed ?_
  rw [isBounded_iff_forall_norm_le]
  exact ⟨1, fun ρ hρ => (Matrix.norm_le_iff zero_le_one).mpr (norm_entry_le_one hρ)⟩

/-! ### The simplex is compact -/

omit [DecidableEq n] in
/-- The simplex is compact: closed, and bounded because every weight lies between zero and one. -/
theorem simplex_isCompact : IsCompact (simplex X) := by
  refine Metric.isCompact_of_isClosed_isBounded ?_ ?_
  · have hpos : IsClosed {N : X → ℝ | ∀ x, 0 ≤ N x} := by
      simp only [Set.ofPred_forall]
      exact isClosed_iInter fun x => isClosed_le continuous_const (continuous_apply x)
    exact hpos.inter (isClosed_eq (by fun_prop) continuous_const)
  · rw [isBounded_iff_forall_norm_le]
    refine ⟨1, fun N hN => (pi_norm_le_iff_of_nonneg zero_le_one).mpr fun x => ?_⟩
    rw [Real.norm_of_nonneg (hN.1 x), ← hN.2]
    exact Finset.single_le_sum (fun y _ => hN.1 y) (Finset.mem_univ x)

omit [DecidableEq n] in
/-- There is a weighting: all the weight on one place. -/
theorem simplex_nonempty [Nonempty X] : (simplex X).Nonempty := by
  classical
  refine ⟨Pi.single (Classical.arbitrary X) 1, fun y => ?_, by simp⟩
  by_cases h : y = Classical.arbitrary X <;> simp [h]

end Saddle
