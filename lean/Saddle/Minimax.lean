import Saddle.Lowest
import Mathlib.Analysis.LocallyConvex.Separation
import Mathlib.Topology.Order.Compact

/-!
# The saddle

Proposition `prop:exchange` of the appendix on the saddle. The least over states of
the greatest local density equals the greatest over weightings of the lowest eigenvalue of the weighted sum
(Eq. `eq:saddle`), and both sides are attained. The proofs are the paper's: every floor lies below every peak; the floor, being the
least value over a compact set of states, is continuous and attains its greatest value `e*` on the simplex;
and if no state held every density at most `e*`, a hyperplane strictly separating the values the states give
the densities from the points at most `e*` everywhere would yield a weighting whose floor exceeds `e*`.
-/

open Matrix
open scoped ComplexOrder

namespace Saddle

variable {n X : Type*} [Fintype n] [DecidableEq n] [Fintype X]

/-! ### Peaks and floors -/

omit [DecidableEq n] in
/-- No density of a state exceeds its peak. -/
theorem expectation_le_peak (D : X → Matrix n n ℂ) (ρ : Matrix n n ℂ) (x : X) :
    expectation ρ (D x) ≤ peak D ρ :=
  le_ciSup (f := fun x => expectation ρ (D x)) (Set.finite_range _).bddAbove x

/-- Every floor lies below every peak: `f(N) ≤ Tr(ρ A_N) = ∑ₓ Nₓ Tr(ρ Dₓ) ≤ maxₓ Tr(ρ Dₓ)`. -/
theorem floor_le_peak [Nonempty n] {D : X → Matrix n n ℂ} (hD : ∀ x, (D x).IsHermitian)
    {ρ : Matrix n n ℂ} {N : X → ℝ} (hρ : ρ ∈ states n) (hN : N ∈ simplex X) :
    floor hD N ≤ peak D ρ :=
  calc floor hD N ≤ expectation ρ (weighted D N) := (lowest_isLeast _).2 ⟨ρ, hρ, rfl⟩
    _ = ∑ x, N x * expectation ρ (D x) := expectation_weighted ρ D N
    _ ≤ ∑ x, N x * peak D ρ :=
        Finset.sum_le_sum fun x _ => mul_le_mul_of_nonneg_left (expectation_le_peak D ρ x) (hN.1 x)
    _ = peak D ρ := by rw [← Finset.sum_mul, hN.2, one_mul]

/-- The floor is the least value the states give the weighted sum, so it is continuous in the weights:
the least value over a compact set of a jointly continuous function. -/
theorem continuous_floor [Nonempty n] {D : X → Matrix n n ℂ} (hD : ∀ x, (D x).IsHermitian) :
    Continuous (floor hD) := by
  have h : floor hD = fun N => sInf ((fun ρ => expectation ρ (weighted D N)) '' states n) :=
    funext fun N => ((lowest_isLeast (weighted_isHermitian hD N)).csInf_eq).symm
  rw [h]
  refine states_isCompact.continuous_sInf (f := fun N ρ => expectation ρ (weighted D N)) ?_
  unfold expectation weighted
  fun_prop

/-! ### The values the states give the densities -/

/-- The values the states give the densities, `V = {(Tr ρ Dₓ)ₓ : ρ a state}`. -/
def values (D : X → Matrix n n ℂ) : Set (X → ℝ) :=
  (fun ρ x => expectation ρ (D x)) '' states n

omit [Fintype X] in
/-- `V` is compact, the continuous image of the compact set of states. -/
theorem values_isCompact (D : X → Matrix n n ℂ) : IsCompact (values D) :=
  states_isCompact.image (continuous_pi fun x => continuous_expectation (D x))

omit [Fintype X] in
/-- `V` is convex, the linear image of the convex set of states. -/
theorem values_convex (D : X → Matrix n n ℂ) : Convex ℝ (values D) :=
  states_convex.is_linear_image
    ⟨fun ρ σ => funext fun x => expectation_add ρ σ (D x), fun c ρ => funext fun x => expectation_smul c ρ (D x)⟩

omit [Fintype n] [DecidableEq n] [Fintype X] in
/-- The points at most `e` in every place form a closed set. -/
theorem below_isClosed (e : ℝ) : IsClosed {t : X → ℝ | ∀ x, t x ≤ e} := by
  simp only [Set.ofPred_forall]
  exact isClosed_iInter fun x => isClosed_le (continuous_apply x) continuous_const

omit [Fintype n] [DecidableEq n] [Fintype X] in
/-- The points at most `e` in every place form a convex set. -/
theorem below_convex (e : ℝ) : Convex ℝ {t : X → ℝ | ∀ x, t x ≤ e} := by
  intro s hs t ht a b ha hb hab x
  calc (a • s + b • t) x = a * s x + b * t x := rfl
    _ ≤ a * e + b * e := add_le_add (mul_le_mul_of_nonneg_left (hs x) ha) (mul_le_mul_of_nonneg_left (ht x) hb)
    _ = e := by rw [← add_mul, hab, one_mul]

/-- **The separation step.** If the weighting `N` certifies the greatest floor `e*`, some state holds every
density at most `e*`. Otherwise `V` and the points at most `e*` everywhere would be disjoint, a hyperplane
would separate them strictly, its normal would be nonnegative, and normalized it would be a weighting whose
floor exceeds `e*`. -/
theorem exists_state_le_floor [Nonempty n] {D : X → Matrix n n ℂ} (hD : ∀ x, (D x).IsHermitian)
    {N : X → ℝ} (hmax : IsMaxOn (floor hD) (simplex X) N) :
    ∃ ρ ∈ states n, ∀ x, expectation ρ (D x) ≤ floor hD N := by
  classical
  set e := floor hD N
  by_contra hnone
  push Not at hnone
  have hdisj : Disjoint (values D) {t : X → ℝ | ∀ x, t x ≤ e} := by
    rw [Set.disjoint_left]
    rintro _ ⟨ρ, hρ, rfl⟩ ht
    obtain ⟨x, hx⟩ := hnone ρ hρ
    exact absurd (ht x) (not_le.mpr hx)
  obtain ⟨f, u, v, hV, huv, hC⟩ :=
    geometric_hahn_banach_compact_closed (values_convex D) (values_isCompact D) (below_convex e) (below_isClosed e) hdisj
  -- the normal of the hyperplane, as weights: `f t = -∑ₓ tₓ wₓ`
  set w : X → ℝ := fun x => -f (fun j => if x = j then 1 else 0)
  have hf : ∀ t, f t = -∑ x, t x * w x := fun t => by
    have hsum := LinearMap.pi_apply_eq_sum_univ (f : (X → ℝ) →ₗ[ℝ] ℝ) t
    simp only [ContinuousLinearMap.coe_coe] at hsum
    rw [hsum, ← Finset.sum_neg_distrib]
    refine Finset.sum_congr rfl fun x _ => ?_
    simp [w, smul_eq_mul]
  set W := ∑ x, w x
  have hone : f (fun _ => e) = -(e * W) := by simp only [hf, W, Finset.mul_sum]
  have hv : v < -(e * W) := hone ▸ hC _ fun _ => le_rfl
  -- each weight is nonnegative: pushing one coordinate down inside the set would cross the hyperplane
  have hw : ∀ x, 0 ≤ w x := by
    intro x
    by_contra hneg
    push Not at hneg
    set s := (v + e * W) / w x
    have hs : 0 ≤ s := div_nonneg_of_nonpos (by linarith) hneg.le
    have hmem : (fun y => e - if y = x then s else 0) ∈ {t : X → ℝ | ∀ x, t x ≤ e} := fun y => by
      dsimp only
      split_ifs <;> linarith
    have hsum : ∑ y, (e - if y = x then s else 0) * w y = e * W - s * w x := by
      simp only [sub_mul, Finset.sum_sub_distrib, ← Finset.mul_sum, ite_mul, zero_mul, Finset.sum_ite_eq',
        Finset.mem_univ, ite_true, W]
    have hsw : s * w x = v + e * W := div_mul_cancel₀ _ hneg.ne
    have hval : f (fun y => e - if y = x then s else 0) = v := by
      rw [hf, hsum, hsw]
      ring
    exact (hC _ hmem).ne' hval
  -- they do not all vanish, since the hyperplane separates something
  have hW : 0 < W := by
    by_contra hW
    have hzero : ∀ x, w x = 0 := fun x =>
      le_antisymm ((Finset.single_le_sum (fun y _ => hw y) (Finset.mem_univ x)).trans (not_lt.mp hW)) (hw x)
    have hW0 : W = 0 := by simp only [W, hzero, Finset.sum_const_zero]
    obtain ⟨ρ, hρ⟩ := states_nonempty (n := n)
    have := hV _ ⟨ρ, hρ, rfl⟩
    simp only [hf, hzero, mul_zero, Finset.sum_const_zero, neg_zero] at this
    rw [hW0, mul_zero, neg_zero] at hv
    linarith
  -- normalized, the normal is a weighting whose floor exceeds `e*`
  set M : X → ℝ := fun x => w x / W
  have hM : M ∈ simplex X := ⟨fun x => div_nonneg (hw x) hW.le, by simp [M, ← Finset.sum_div, W, hW.ne']⟩
  obtain ⟨⟨σ, hσ, hσeq⟩, -⟩ := lowest_isLeast (weighted_isHermitian hD M)
  have hσV := hV _ ⟨σ, hσ, rfl⟩
  rw [hf] at hσV
  have hfloor : floor hD M = (∑ x, expectation σ (D x) * w x) / W := by
    rw [floor, ← hσeq]
    dsimp only
    rw [expectation_weighted, Finset.sum_div]
    exact Finset.sum_congr rfl fun x _ => by simp only [M]; ring
  have hgt : e < floor hD M := by
    rw [hfloor, lt_div_iff₀ hW]
    linarith
  exact absurd (isMaxOn_iff.mp hmax M hM) (not_le.mpr hgt)

/-! ### The propositions -/

/-- **Proposition `prop:exchange` (the exchange).** Both sides of Eq. `eq:saddle` are attained, and they are
equal: some state has the least peak over all states, some weighting on the simplex certifies the greatest
floor over all weightings, and the least peak equals the greatest floor. -/
theorem saddle_value [Nonempty n] [Nonempty X] {D : X → Matrix n n ℂ} (hD : ∀ x, (D x).IsHermitian) :
    ∃ ρ ∈ states n, ∃ N ∈ simplex X,
      IsMinOn (peak D) (states n) ρ ∧ IsMaxOn (floor hD) (simplex X) N ∧ peak D ρ = floor hD N := by
  obtain ⟨N, hN, hmax⟩ := simplex_isCompact.exists_isMaxOn simplex_nonempty (continuous_floor hD).continuousOn
  obtain ⟨ρ, hρ, hle⟩ := exists_state_le_floor hD hmax
  have hpeak : peak D ρ ≤ floor hD N := ciSup_le hle
  refine ⟨ρ, hρ, N, hN, isMinOn_iff.mpr fun σ hσ => hpeak.trans (floor_le_peak hD hσ hN), hmax,
    le_antisymm hpeak (floor_le_peak hD hρ hN)⟩

end Saddle
