import Saddle.Minimax

/-!
# Real states suffice for real densities

Lemma `lem:real` of the appendix on the two programs: when every density is real symmetric, real symmetric states
attain the least greatest local density. The proof is the paper's: a state's transpose is its complex conjugate,
also a state, and it gives each real symmetric density the same value; so the average of a state and its transpose
is a real symmetric state with the same densities.
-/

open Matrix
open scoped ComplexOrder

namespace Saddle

variable {n X : Type*} [Fintype n] [DecidableEq n] [Fintype X]

omit [DecidableEq n] in
/-- A state's transpose gives a symmetric operator the value the state gives it. -/
theorem expectation_transpose {ρ A : Matrix n n ℂ} (hA : Aᵀ = A) : expectation ρᵀ A = expectation ρ A := by
  unfold expectation
  rw [← trace_transpose, transpose_mul, transpose_transpose, hA, trace_mul_comm]

/-- A state's transpose is a state. -/
theorem transpose_mem_states {ρ : Matrix n n ℂ} (hρ : ρ ∈ states n) : ρᵀ ∈ states n :=
  ⟨hρ.1.transpose, by rw [trace_transpose, hρ.2]⟩

/-- **Lemma `lem:real`.** If every density is real symmetric, a real symmetric state attains the left side of
Eq. `eq:saddle`: some state that is its own transpose and its own complex conjugate has the least peak over all
states. -/
theorem real_states_attain [Nonempty n] [Nonempty X] {D : X → Matrix n n ℂ} (hD : ∀ x, (D x).IsHermitian)
    (hsymmetric : ∀ x, (D x)ᵀ = D x) :
    ∃ ρ ∈ states n, ρᵀ = ρ ∧ ρ.map star = ρ ∧ IsMinOn (peak D) (states n) ρ := by
  obtain ⟨ρ₀, hρ₀, -, -, hmin, -, -⟩ := saddle_value hD
  set ρ := (1 / 2 : ℝ) • ρ₀ + (1 / 2 : ℝ) • ρ₀ᵀ
  have hρ : ρ ∈ states n :=
    states_convex hρ₀ (transpose_mem_states hρ₀) (by norm_num) (by norm_num) (by norm_num)
  have hsame : ∀ x, expectation ρ (D x) = expectation ρ₀ (D x) := fun x => by
    rw [expectation_add, expectation_smul, expectation_smul, expectation_transpose (hsymmetric x)]
    ring
  have hpeak : peak D ρ = peak D ρ₀ := by simp only [peak, hsame]
  have hsym : ρᵀ = ρ := by
    simp only [ρ, transpose_add, transpose_smul, transpose_transpose]
    exact add_comm _ _
  refine ⟨ρ, hρ, hsym, ?_, isMinOn_iff.mpr fun σ hσ => hpeak ▸ isMinOn_iff.mp hmin σ hσ⟩
  -- Hermitian and symmetric, so real: `star ρᵢⱼ = ρⱼᵢ = ρᵢⱼ`
  ext i j
  have hji : ρ j i = ρ i j := by simpa [transpose_apply] using congrFun (congrFun hsym i) j
  calc (ρ.map star) i j = star (ρ j i) := by rw [map_apply, hji]
    _ = ρ i j := hρ.1.isHermitian.apply i j

end Saddle
