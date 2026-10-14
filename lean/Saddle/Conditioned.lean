import Mathlib.Analysis.Normed.Algebra.MatrixExponential
import Mathlib.LinearAlgebra.Matrix.Kronecker

/-!
# Conditioned evolution

Proposition `prop:conditioned` of the appendix on probes. The whole is the rest, indexed by `r`, beside a clock,
indexed by `c`; the rest's generator `A` and the clock's energy `H` act on different factors, `A ⊗ I` and `I ⊗ H`, so
they commute. If the whole's state `Ψ` is an eigenvector of their sum with eigenvalue `λ`, then running the rest
forward is running the clock back, and the rest conditioned on a clock bra `β` read after the clock has run,
`β e^{-sH}`, is the rest conditioned on `β` evolved by `e^{sA}`, up to the factor `e^{sλ}` (Eq. `eq:conditioned`, with
`s = -it`). The proof is the paper's: the exponential of the sum is the product of the exponentials, and an
eigenvector of the sum stays one of its exponential.
-/

open Matrix
open scoped Kronecker Matrix.Norms.Operator

namespace Saddle

variable {r c : Type*} [Fintype r] [DecidableEq r] [Fintype c] [DecidableEq c]

/-- The rest conditioned on a clock bra `β`: `ψ(i) = ∑ₖ βₖ Ψ(i, k)`. -/
def conditioned (β : c → ℂ) (Ψ : r × c → ℂ) : r → ℂ :=
  fun i => ∑ k, β k * Ψ (i, k)

omit [Fintype r] [DecidableEq r] [DecidableEq c] in
/-- Conditioning is linear in the whole's state. -/
theorem conditioned_smul (β : c → ℂ) (a : ℂ) (Ψ : r × c → ℂ) :
    conditioned β (a • Ψ) = a • conditioned β Ψ := by
  funext i
  simp only [conditioned, Pi.smul_apply, smul_eq_mul, Finset.mul_sum]
  exact Finset.sum_congr rfl fun k _ => by ring

omit [DecidableEq r] in
/-- Conditioning commutes with what acts on the rest alone. -/
theorem conditioned_kronecker_one (β : c → ℂ) (B : Matrix r r ℂ) (Ψ : r × c → ℂ) :
    conditioned β ((B ⊗ₖ (1 : Matrix c c ℂ)) *ᵥ Ψ) = B *ᵥ conditioned β Ψ := by
  funext i
  simp only [conditioned, mulVec, dotProduct, Fintype.sum_prod_type, kroneckerMap_apply, one_apply, mul_ite,
    mul_one, mul_zero, ite_mul, zero_mul, Finset.sum_ite_eq, Finset.mem_univ, ite_true, Finset.mul_sum]
  rw [Finset.sum_comm]
  exact Finset.sum_congr rfl fun j _ => Finset.sum_congr rfl fun k _ => by ring

omit [DecidableEq c] in
/-- What acts on the clock alone acts at each of the rest's levels: `((I ⊗ C) Ψ)(i, k) = ∑ₗ C_{kl} Ψ(i, l)`. -/
theorem one_kronecker_mulVec (C : Matrix c c ℂ) (Ψ : r × c → ℂ) (i : r) (k : c) :
    (((1 : Matrix r r ℂ) ⊗ₖ C) *ᵥ Ψ) (i, k) = ∑ l, C k l * Ψ (i, l) := by
  simp only [mulVec, dotProduct, Fintype.sum_prod_type, kroneckerMap_apply, one_apply, ite_mul, one_mul, zero_mul]
  rw [Finset.sum_comm]
  simp only [Finset.sum_ite_eq, Finset.mem_univ, ite_true]

omit [DecidableEq c] in
/-- What acts on the clock alone moves to the bra. -/
theorem conditioned_one_kronecker (β : c → ℂ) (C : Matrix c c ℂ) (Ψ : r × c → ℂ) :
    conditioned β (((1 : Matrix r r ℂ) ⊗ₖ C) *ᵥ Ψ) = conditioned (β ᵥ* C) Ψ := by
  funext i
  simp only [conditioned, one_kronecker_mulVec, vecMul, dotProduct, Finset.mul_sum, Finset.sum_mul]
  rw [Finset.sum_comm]
  exact Finset.sum_congr rfl fun l _ => Finset.sum_congr rfl fun k _ => by ring

/-- `A ↦ A ⊗ I`, a ring homomorphism. -/
def kroneckerOne : Matrix r r ℂ →+* Matrix (r × c) (r × c) ℂ where
  toFun A := A ⊗ₖ (1 : Matrix c c ℂ)
  map_one' := one_kronecker_one
  map_mul' A B := by rw [← mul_kronecker_mul, Matrix.mul_one]
  map_zero' := by ext; simp
  map_add' A B := add_kronecker A B 1

/-- `C ↦ I ⊗ C`, a ring homomorphism. -/
def oneKronecker : Matrix c c ℂ →+* Matrix (r × c) (r × c) ℂ where
  toFun C := (1 : Matrix r r ℂ) ⊗ₖ C
  map_one' := one_kronecker_one
  map_mul' A B := by rw [← mul_kronecker_mul, Matrix.mul_one]
  map_zero' := by ext; simp
  map_add' A B := kronecker_add 1 A B

/-- The exponential of what acts on the rest alone acts on the rest alone: `e^{A ⊗ I} = e^A ⊗ I`. -/
theorem exp_kronecker_one (A : Matrix r r ℂ) :
    NormedSpace.exp (A ⊗ₖ (1 : Matrix c c ℂ)) = NormedSpace.exp A ⊗ₖ (1 : Matrix c c ℂ) := by
  have hc : Continuous (kroneckerOne (r := r) (c := c)) :=
    continuous_pi fun i => continuous_pi fun j =>
      ((continuous_apply j.1).comp (continuous_apply i.1)).mul continuous_const
  exact (NormedSpace.map_exp (kroneckerOne (r := r) (c := c)) hc A).symm

/-- The exponential of what acts on the clock alone acts on the clock alone: `e^{I ⊗ C} = I ⊗ e^C`. -/
theorem exp_one_kronecker (C : Matrix c c ℂ) :
    NormedSpace.exp ((1 : Matrix r r ℂ) ⊗ₖ C) = (1 : Matrix r r ℂ) ⊗ₖ NormedSpace.exp C := by
  have hc : Continuous (oneKronecker (r := r) (c := c)) :=
    continuous_pi fun i => continuous_pi fun j =>
      continuous_const.mul ((continuous_apply j.2).comp (continuous_apply i.2))
  exact (NormedSpace.map_exp (oneKronecker (r := r) (c := c)) hc C).symm

/-- An eigenvector of a matrix is one of its exponential: `M v = μ v` gives `e^M v = e^μ v`. -/
theorem exp_mulVec_of_eigen {m : Type*} [Fintype m] [DecidableEq m] {M : Matrix m m ℂ} {v : m → ℂ} {μ : ℂ}
    (h : M *ᵥ v = μ • v) : NormedSpace.exp M *ᵥ v = NormedSpace.exp μ • v := by
  rcases isEmpty_or_nonempty m with hm | ⟨⟨j⟩⟩
  · exact funext fun i => (IsEmpty.false i).elim
  -- every column of `X` is `v`, so `X (μ I) = M X`, and the exponential carries the relation over
  set X : Matrix m m ℂ := Matrix.of fun i _ => v i
  have hX : X * (μ • (1 : Matrix m m ℂ)) = M * X := by
    rw [Matrix.mul_smul, Matrix.mul_one]
    ext i k
    have hi := congrFun h i
    simp only [mulVec, dotProduct, Pi.smul_apply, smul_eq_mul] at hi
    simp only [X, Matrix.smul_apply, of_apply, smul_eq_mul, mul_apply, hi]
  have hexp : X * NormedSpace.exp (μ • (1 : Matrix m m ℂ)) = NormedSpace.exp M * X :=
    (SemiconjBy.exp_right (x := X) hX)
  have hscalar : NormedSpace.exp (μ • (1 : Matrix m m ℂ)) = NormedSpace.exp μ • (1 : Matrix m m ℂ) := by
    rw [← Algebra.algebraMap_eq_smul_one, ← Algebra.algebraMap_eq_smul_one]
    exact (NormedSpace.algebraMap_exp_comm μ).symm
  rw [hscalar, Matrix.mul_smul, Matrix.mul_one] at hexp
  funext i
  have hij := congrFun (congrFun hexp.symm i) j
  simp only [mul_apply, X, of_apply, Matrix.smul_apply, smul_eq_mul] at hij
  simpa [mulVec, dotProduct] using hij

/-- **Proposition `prop:conditioned`.** If the whole is an eigenvector of the rest's generator plus the clock's
energy, `(A ⊗ I + I ⊗ H) Ψ = λ Ψ`, then for every `s`, `e^{sλ} ψ(β e^{-sH}) = e^{sA} ψ(β)`: the rest conditioned on
the clock read after it has run is the rest conditioned before, evolved by the rest's generator (Eq.
`eq:conditioned` with `s = -it`). -/
theorem conditioned_evolution (A : Matrix r r ℂ) (H : Matrix c c ℂ) (Ψ : r × c → ℂ) (ev : ℂ)
    (hΨ : (A ⊗ₖ (1 : Matrix c c ℂ) + (1 : Matrix r r ℂ) ⊗ₖ H) *ᵥ Ψ = ev • Ψ) (β : c → ℂ) (s : ℂ) :
    NormedSpace.exp (s * ev) • conditioned (β ᵥ* NormedSpace.exp (-(s • H))) Ψ =
      NormedSpace.exp (s • A) *ᵥ conditioned β Ψ := by
  set P := A ⊗ₖ (1 : Matrix c c ℂ)
  set Q := (1 : Matrix r r ℂ) ⊗ₖ H
  have hPQ : Commute P Q := by
    simp only [Commute, SemiconjBy, P, Q, ← mul_kronecker_mul, Matrix.mul_one, Matrix.one_mul]
  -- the clock run back keeps the whole an eigenvector of the sum, and so of its exponential
  have hcomm : Commute (P + Q) (NormedSpace.exp (-(s • Q))) :=
    (((hPQ.add_left (Commute.refl Q)).smul_right s).neg_right).exp_right
  have hback : (s • (P + Q)) *ᵥ (NormedSpace.exp (-(s • Q)) *ᵥ Ψ) =
      (s * ev) • (NormedSpace.exp (-(s • Q)) *ᵥ Ψ) := by
    rw [smul_mulVec, mulVec_mulVec, hcomm.eq, ← mulVec_mulVec, hΨ, mulVec_smul, smul_smul]
  -- running the rest forward is running the clock back
  have hsplit : NormedSpace.exp (s • P) = NormedSpace.exp (s • (P + Q)) * NormedSpace.exp (-(s • Q)) := by
    rw [smul_add, Matrix.exp_add_of_commute _ _ ((hPQ.smul_left s).smul_right s), Matrix.mul_assoc,
      ← Matrix.exp_add_of_commute _ _ (Commute.refl (s • Q)).neg_right, add_neg_cancel, NormedSpace.exp_zero,
      Matrix.mul_one]
  have hrun : NormedSpace.exp (s • P) *ᵥ Ψ = NormedSpace.exp (s * ev) • (NormedSpace.exp (-(s • Q)) *ᵥ Ψ) := by
    rw [hsplit, ← mulVec_mulVec, exp_mulVec_of_eigen hback]
  -- condition both sides on the clock bra
  have hP : NormedSpace.exp (s • P) = NormedSpace.exp (s • A) ⊗ₖ (1 : Matrix c c ℂ) := by
    rw [← exp_kronecker_one, smul_kronecker]
  have hQ : NormedSpace.exp (-(s • Q)) = (1 : Matrix r r ℂ) ⊗ₖ NormedSpace.exp (-(s • H)) := by
    have hneg : -(s • Q) = (1 : Matrix r r ℂ) ⊗ₖ (-(s • H)) := by
      ext ⟨i, k⟩ ⟨j, l⟩
      simp only [Q, Matrix.neg_apply, Matrix.smul_apply, kroneckerMap_apply, smul_eq_mul]
      ring
    rw [hneg, exp_one_kronecker]
  have hcond := congrArg (conditioned β) hrun
  rw [hP, conditioned_kronecker_one, conditioned_smul, hQ, conditioned_one_kronecker] at hcond
  exact hcond.symm

end Saddle
