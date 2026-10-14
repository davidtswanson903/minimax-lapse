import Saddle.Conditions
import Mathlib.Analysis.Normed.Algebra.MatrixExponential

/-!
# The weighted evolution

Corollary `cor:evolution` of the appendix on the saddle. At a saddle pair the state commutes with the weighted
Hamiltonian `A_{N*}`, so it is stationary under `e^{-iA_{N*}t}`: the first condition gives `A ρ = e* ρ`, and its
adjoint `ρ A = e* ρ`. And when each density is its place's own piece plus half of every link piece at it, the
weighted Hamiltonian is each own piece at its place's weight plus each link at the mean weight of its two ends
(Eq. `eq:evolution`), so a place with no weight has its own piece frozen.
-/

open Matrix NormedSpace
open scoped ComplexOrder

namespace Saddle

variable {n X : Type*} [Fintype n] [DecidableEq n] [Fintype X]

/-- **Corollary `cor:evolution`, first part.** At a saddle pair the state commutes with the weighted Hamiltonian,
and so is stationary under the evolution it generates, `e^{-iA_{N*}t} ρ e^{iA_{N*}t} = ρ` for every `t`. -/
theorem saddle_evolution [Nonempty n] [Nonempty X] {D : X → Matrix n n ℂ} (hD : ∀ x, (D x).IsHermitian)
    {ρ : Matrix n n ℂ} {N : X → ℝ} (hρ : ρ ∈ states n) (hN : N ∈ simplex X)
    (hmin : IsMinOn (peak D) (states n) ρ) (hmax : IsMaxOn (floor hD) (simplex X) N) :
    Commute (weighted D N) ρ ∧
      ∀ t : ℝ, exp (-((t : ℂ) * Complex.I) • weighted D N) * ρ * exp (((t : ℂ) * Complex.I) • weighted D N) = ρ := by
  have held := (saddle_conditions hD hρ hN hmin hmax).1
  -- the adjoint of `A ρ = e* ρ` is `ρ A = e* ρ`, since `A`, `ρ` and `e*` are self-adjoint
  have hadjoint : ρ * weighted D N = (peak D ρ : ℂ) • ρ := by
    have h := congrArg conjTranspose held
    rwa [conjTranspose_mul, hρ.1.isHermitian.eq, (weighted_isHermitian hD N).eq, conjTranspose_smul,
      Complex.star_def, Complex.conj_ofReal, hρ.1.isHermitian.eq] at h
  have hcommute : Commute (weighted D N) ρ := held.trans hadjoint.symm
  refine ⟨hcommute, fun t => ?_⟩
  set B := ((t : ℂ) * Complex.I) • weighted D N
  have hB : Commute ρ (exp B) := (hcommute.smul_left ((t : ℂ) * Complex.I)).symm.exp_right
  have hone : exp (-B) * exp B = 1 := by
    rw [← Matrix.exp_add_of_commute _ _ (Commute.neg_left (Commute.refl B)), neg_add_cancel, exp_zero]
  rw [neg_smul, mul_assoc, hB.eq, ← mul_assoc, hone, one_mul]

omit [Fintype n] [DecidableEq n] in
/-- **Corollary `cor:evolution`, second part (Eq. `eq:evolution`).** When each density is its place's own piece
plus half of every link piece at it, the weighted sum is each own piece at its place's weight plus each link piece
at the mean weight of its two ends; so a place with no weight contributes no piece of its own. -/
theorem weighted_of_pieces [DecidableEq X] (own : X → Matrix n n ℂ) (links : Finset (X × X))
    (link : X × X → Matrix n n ℂ) (hlinks : ∀ l ∈ links, l.1 ≠ l.2) (N : X → ℝ) :
    weighted (fun x => own x + (1 / 2 : ℂ) • ∑ l ∈ links.filter (fun l => l.1 = x ∨ l.2 = x), link l) N =
      ∑ x, (N x : ℂ) • own x + ∑ l ∈ links, (((N l.1 + N l.2) / 2 : ℝ) : ℂ) • link l := by
  unfold weighted
  simp only [smul_add, Finset.sum_add_distrib, Finset.sum_filter, Finset.smul_sum]
  congr 1
  rw [Finset.sum_comm]
  refine Finset.sum_congr rfl fun l hl => ?_
  have hpair : (Finset.univ.filter fun x : X => l.1 = x ∨ l.2 = x) = {l.1, l.2} := by
    ext x
    simp only [Finset.mem_filter, Finset.mem_univ, true_and, Finset.mem_insert, Finset.mem_singleton]
    tauto
  calc ∑ x, (N x : ℂ) • (1 / 2 : ℂ) • (if l.1 = x ∨ l.2 = x then link l else 0)
      = ∑ x ∈ Finset.univ.filter (fun x : X => l.1 = x ∨ l.2 = x), (N x : ℂ) • (1 / 2 : ℂ) • link l := by
        rw [Finset.sum_filter]
        exact Finset.sum_congr rfl fun x _ => by split_ifs <;> simp
    _ = (((N l.1 + N l.2) / 2 : ℝ) : ℂ) • link l := by
        rw [hpair, Finset.sum_pair (hlinks l hl), ← add_smul, smul_smul]
        congr 1
        push_cast
        ring

end Saddle
