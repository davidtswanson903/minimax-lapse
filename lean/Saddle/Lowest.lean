import Saddle.States
import Mathlib.Analysis.Matrix.Spectrum
import Mathlib.Analysis.SpecialFunctions.ContinuousFunctionalCalculus.Rpow.Basic

/-!
# The lowest value a state gives a Hermitian matrix

Lemma `lem:lowest` of the appendix on the saddle: over the states, the least value of `Tr(ρ A)` is the
lowest eigenvalue of `A`, and a state attains it exactly when `A ρ = λ_min(A) ρ`. The proof is the paper's:
`B = A - λ_min(A) I` is positive semidefinite, so `Tr(ρ B) = Tr(B^{1/2} ρ B^{1/2}) ≥ 0`; equality forces
`B^{1/2} ρ^{1/2} = 0` and so `B ρ = 0`; and a state along an eigenvector of the lowest eigenvalue attains it.
-/

open Matrix
open Unitary
open scoped ComplexOrder MatrixOrder

namespace Saddle

variable {n : Type*} [Fintype n] [DecidableEq n]

/-- No eigenvalue of a Hermitian matrix lies below its lowest. -/
theorem lowest_le_eigenvalues {A : Matrix n n ℂ} (hA : A.IsHermitian) (i : n) :
    lowest hA ≤ hA.eigenvalues i :=
  ciInf_le (Set.finite_range _).bddBelow i

/-- The lowest eigenvalue of a Hermitian matrix is one of its eigenvalues. -/
theorem exists_eigenvalues_eq_lowest [Nonempty n] {A : Matrix n n ℂ} (hA : A.IsHermitian) :
    ∃ i, hA.eigenvalues i = lowest hA :=
  exists_eq_ciInf_of_finite

/-- `A - λ_min(A) I` is positive semidefinite: in the basis of eigenvectors it is diagonal, with the
nonnegative entries `λᵢ - λ_min(A)`. -/
theorem sub_lowest_posSemidef {A : Matrix n n ℂ} (hA : A.IsHermitian) :
    (A - (lowest hA : ℂ) • 1).PosSemidef := by
  have hle := lowest_le_eigenvalues hA
  generalize hc : ((lowest hA : ℝ) : ℂ) = c
  have hshift : A - c • 1 = conjStarAlgAut ℂ _ hA.eigenvectorUnitary
      (diagonal (RCLike.ofReal ∘ hA.eigenvalues) - c • 1) := by
    conv_lhs => rw [hA.spectral_theorem]
    rw [map_sub, map_smul, map_one]
  rw [hshift]
  simp only [conjStarAlgAut_apply, isUnit_coe.posSemidef_star_right_conjugate_iff, smul_one_eq_diagonal,
    diagonal_sub, posSemidef_diagonal_iff, ← hc]
  intro i
  simpa using sub_nonneg.mpr (Complex.real_le_real.mpr (hle i))

/-- For positive semidefinite `ρ` and `B`, `Tr(ρ B)` is real and nonnegative, being `Tr(B^{1/2} ρ B^{1/2})`;
and when it vanishes, `B ρ = 0`. -/
theorem trace_mul_re_nonneg {ρ B : Matrix n n ℂ} (hρ : ρ.PosSemidef) (hB : B.PosSemidef) :
    0 ≤ (ρ * B).trace.re ∧ ((ρ * B).trace.re = 0 → B * ρ = 0) := by
  set S := CFC.sqrt B
  set T := CFC.sqrt ρ
  have hS : S * S = B := CFC.sqrt_mul_sqrt_self B hB.nonneg
  have hT : T * T = ρ := CFC.sqrt_mul_sqrt_self ρ hρ.nonneg
  have hSh : Sᴴ = S := (CFC.sqrt_nonneg B).posSemidef.isHermitian.eq
  have hTh : Tᴴ = T := (CFC.sqrt_nonneg ρ).posSemidef.isHermitian.eq
  have htrace : (ρ * B).trace = ((S * T) * (S * T)ᴴ).trace := by
    rw [conjTranspose_mul, hSh, hTh, ← hS, ← hT, show T * T * (S * S) = T * T * S * S by simp only [mul_assoc],
      trace_mul_comm]
    simp only [mul_assoc]
  have hpsd := posSemidef_self_mul_conjTranspose (S * T)
  have hnonneg := Complex.le_def.mp hpsd.trace_nonneg
  refine ⟨by rw [htrace]; exact hnonneg.1, fun h0 => ?_⟩
  have hST : S * T = 0 := by
    refine trace_mul_conjTranspose_self_eq_zero_iff.mp ?_
    rw [← htrace]
    exact Complex.ext (by simpa using h0) (by rw [htrace]; simpa using hnonneg.2.symm)
  rw [← hS, ← hT, show S * S * (T * T) = S * (S * T) * T by simp only [mul_assoc], hST, mul_zero, zero_mul]

/-- The state along a unit eigenvector of a Hermitian matrix: the projection `v v†` onto it. -/
noncomputable def eigenstate {A : Matrix n n ℂ} (hA : A.IsHermitian) (i : n) : Matrix n n ℂ :=
  vecMulVec ⇑(hA.eigenvectorBasis i) (star ⇑(hA.eigenvectorBasis i))

omit [DecidableEq n] in
/-- The value the projection `v v†` gives `A` is `v† A v`. -/
theorem expectation_vecMulVec (v : n → ℂ) (A : Matrix n n ℂ) :
    expectation (vecMulVec v (star v)) A = (star v ⬝ᵥ (A *ᵥ v)).re := by
  rw [expectation, vecMulVec_mul, trace_vecMulVec, dotProduct_comm, ← dotProduct_mulVec]

/-- The state along a unit eigenvector is a state. -/
theorem eigenstate_mem_states {A : Matrix n n ℂ} (hA : A.IsHermitian) (i : n) : eigenstate hA i ∈ states n := by
  refine ⟨posSemidef_vecMulVec_self_star _, ?_⟩
  have hunit : inner ℂ (hA.eigenvectorBasis i) (hA.eigenvectorBasis i) = 1 := by
    rw [inner_self_eq_norm_sq_to_K, hA.eigenvectorBasis.orthonormal.1 i]
    simp
  rw [EuclideanSpace.inner_eq_star_dotProduct] at hunit
  rw [eigenstate, trace_vecMulVec]
  exact hunit

/-- The state along the eigenvector of `λᵢ` gives `A` the value `λᵢ`. -/
theorem expectation_eigenstate {A : Matrix n n ℂ} (hA : A.IsHermitian) (i : n) :
    expectation (eigenstate hA i) A = hA.eigenvalues i := by
  rw [hA.eigenvalues_eq i, eigenstate, expectation_vecMulVec]
  rfl

/-- `Tr(ρ (A - c I)) = Tr(ρ A) - c` for a state `ρ`. -/
theorem expectation_sub_smul_one {ρ : Matrix n n ℂ} (hρ : ρ ∈ states n) (A : Matrix n n ℂ) (c : ℝ) :
    (ρ * (A - (c : ℂ) • 1)).trace.re = expectation ρ A - c := by
  rw [Matrix.mul_sub, trace_sub, Matrix.mul_smul, Matrix.mul_one, trace_smul, hρ.2, smul_eq_mul, mul_one,
    Complex.sub_re, Complex.ofReal_re, expectation]

/-- **Lemma `lem:lowest`, first part.** For every Hermitian `A`, the least value a state gives `A` is its
lowest eigenvalue: some state attains it, and no state goes below it. -/
theorem lowest_isLeast [Nonempty n] {A : Matrix n n ℂ} (hA : A.IsHermitian) :
    IsLeast ((fun ρ => expectation ρ A) '' states n) (lowest hA) := by
  refine ⟨?_, ?_⟩
  · obtain ⟨i, hi⟩ := exists_eigenvalues_eq_lowest hA
    exact ⟨eigenstate hA i, eigenstate_mem_states hA i, by simp only [expectation_eigenstate, hi]⟩
  · rintro _ ⟨ρ, hρ, rfl⟩
    have h := (trace_mul_re_nonneg hρ.1 (sub_lowest_posSemidef hA)).1
    rw [expectation_sub_smul_one hρ] at h
    linarith

/-- **Lemma `lem:lowest`, second part.** A state attains the lowest eigenvalue of a Hermitian `A` exactly
when `A ρ = λ_min(A) ρ`. -/
theorem expectation_eq_lowest_iff {A : Matrix n n ℂ} (hA : A.IsHermitian) {ρ : Matrix n n ℂ}
    (hρ : ρ ∈ states n) :
    expectation ρ A = lowest hA ↔ A * ρ = (lowest hA : ℂ) • ρ := by
  constructor
  · intro h
    have hzero := (trace_mul_re_nonneg hρ.1 (sub_lowest_posSemidef hA)).2
      (by rw [expectation_sub_smul_one hρ, h, sub_self])
    rw [Matrix.sub_mul, Matrix.smul_mul, Matrix.one_mul, sub_eq_zero] at hzero
    exact hzero
  · intro h
    rw [expectation, trace_mul_comm, h, trace_smul, hρ.2, smul_eq_mul, mul_one, Complex.ofReal_re]

/-- The lowest eigenvalue depends only on the matrix, not on the proof that it is Hermitian. -/
theorem lowest_congr {A B : Matrix n n ℂ} (hA : A.IsHermitian) (hB : B.IsHermitian) (h : A = B) :
    lowest hA = lowest hB := by
  subst h
  rfl

end Saddle
