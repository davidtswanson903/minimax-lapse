import Saddle.Forbids

/-!
# The drop-out depth

Proposition `prop:dropout` of the appendix on the drop-out depth. Deepen place `y`: replace `D_y` by `D_y - μ P`
for a Hermitian `P` and a depth `μ`, leaving every other density unchanged. The other places' problem does not
involve `μ`; let `e'` be its value, the least over states of the greatest density among the other places, and `S'`
the states holding every other place at or below `e'`, which are its optimal states.

Some optimal weighting gives `y` no weight exactly when some state of `S'` holds the deepened `y` at or below `e'`,
and then the common value is `e'`. If every optimal weighting gives `y` weight at depth zero, `y` carries none
exactly from the depth `μ_c` of Eq. `eq:dropout` on, the least ratio `(Tr(ρ D_y) - e') / Tr(ρ P)` over the states of
`S'` with `Tr(ρ P) > 0`; beyond it every optimal weighting gives `y` none, and the common value stays `e'`. If no
state of `S'` has `Tr(ρ P) > 0`, `y` never drops out. The proofs are the paper's.
-/

open Matrix Set
open scoped ComplexOrder

namespace Saddle

variable {n X : Type*} [Fintype n] [DecidableEq n] [Fintype X] [DecidableEq X]

/-- The densities with place `y` deepened by `μ` along `P`: `D_y - μ P`, every other place's unchanged. -/
noncomputable def deepened (D : X → Matrix n n ℂ) (P : Matrix n n ℂ) (y : X) (μ : ℝ) : X → Matrix n n ℂ :=
  fun x => if x = y then D x - (μ : ℂ) • P else D x

/-- The other places' densities, place `y` left out. -/
def others (D : X → Matrix n n ℂ) (y : X) : {x // x ≠ y} → Matrix n n ℂ :=
  fun x => D x

/-- The states holding every place other than `y` at or below `e`. -/
def held (D : X → Matrix n n ℂ) (y : X) (e : ℝ) : Set (Matrix n n ℂ) :=
  {ρ | ρ ∈ states n ∧ ∀ x, x ≠ y → expectation ρ (D x) ≤ e}

omit [Fintype n] [DecidableEq n] [Fintype X] in
/-- The deepened densities are Hermitian. -/
theorem deepened_isHermitian {D : X → Matrix n n ℂ} (hD : ∀ x, (D x).IsHermitian) {P : Matrix n n ℂ}
    (hP : P.IsHermitian) (y : X) (μ : ℝ) : ∀ x, (deepened D P y μ x).IsHermitian := fun x => by
  unfold deepened
  split_ifs
  · unfold Matrix.IsHermitian
    rw [conjTranspose_sub, (hD x).eq, conjTranspose_smul, hP.eq, Complex.star_def, Complex.conj_ofReal]
  · exact hD x

omit [Fintype n] [DecidableEq n] [Fintype X] [DecidableEq X] in
/-- The other places' densities are Hermitian. -/
theorem others_isHermitian {D : X → Matrix n n ℂ} (hD : ∀ x, (D x).IsHermitian) (y : X) :
    ∀ x, (others D y x).IsHermitian := fun x => hD x

omit [DecidableEq n] [Fintype X] in
/-- The value a state gives the deepened place: `Tr(ρ D_y) - μ Tr(ρ P)`. -/
theorem expectation_deepened_self (ρ : Matrix n n ℂ) (D : X → Matrix n n ℂ) (P : Matrix n n ℂ) (y : X) (μ : ℝ) :
    expectation ρ (deepened D P y μ y) = expectation ρ (D y) - μ * expectation ρ P := by
  simp only [deepened, ↓reduceIte, expectation, Matrix.mul_sub, trace_sub, Matrix.mul_smul, trace_smul,
    Complex.sub_re, smul_eq_mul, Complex.re_ofReal_mul]

omit [DecidableEq n] [Fintype X] in
/-- Every other place keeps its value. -/
theorem expectation_deepened_ne (ρ : Matrix n n ℂ) (D : X → Matrix n n ℂ) (P : Matrix n n ℂ) {x y : X} (μ : ℝ)
    (hx : x ≠ y) : expectation ρ (deepened D P y μ x) = expectation ρ (D x) := by
  simp only [deepened, hx, ↓reduceIte]

/-! ### Weightings that give `y` no weight -/

omit [Fintype n] [DecidableEq n] in
/-- A weighting that gives `y` no weight, read on the other places, is a weighting of them. -/
theorem restrict_mem_simplex {N : X → ℝ} (hN : N ∈ simplex X) {y : X} (hy : N y = 0) :
    (fun x : {x // x ≠ y} => N x) ∈ simplex {x // x ≠ y} := by
  refine ⟨fun x => hN.1 x, ?_⟩
  have h := hN.2
  rwa [Fintype.sum_eq_add_sum_subtype_ne _ y, hy, zero_add] at h

omit [Fintype n] [DecidableEq n] in
/-- Where `y` has no weight, the weighted sum is the other places' weighted sum. -/
theorem weighted_deepened_of_zero (D : X → Matrix n n ℂ) (P : Matrix n n ℂ) (y : X) (μ : ℝ) {N : X → ℝ}
    (hy : N y = 0) : weighted (deepened D P y μ) N = weighted (others D y) fun x => N x := by
  rw [weighted, Fintype.sum_eq_add_sum_subtype_ne _ y, hy, Complex.ofReal_zero, zero_smul, zero_add, weighted]
  exact Finset.sum_congr rfl fun x _ => by simp only [deepened, x.2, ↓reduceIte, others]

/-- A weighting of the other places, extended by no weight on `y`. -/
noncomputable def extend (y : X) (N : {x // x ≠ y} → ℝ) : X → ℝ :=
  fun x => if h : x = y then 0 else N ⟨x, h⟩

omit [Fintype X] in
/-- Extended, the weighting gives `y` no weight. -/
theorem extend_self (y : X) (N : {x // x ≠ y} → ℝ) : extend y N y = 0 := by
  simp [extend]

omit [Fintype X] in
/-- Extended and read on the other places, the weighting is unchanged. -/
theorem restrict_extend (y : X) (N : {x // x ≠ y} → ℝ) : (fun x : {x // x ≠ y} => extend y N x) = N :=
  funext fun x => by simp [extend, x.2]

omit [Fintype n] [DecidableEq n] in
/-- Extended, a weighting of the other places is a weighting. -/
theorem extend_mem_simplex {y : X} {N : {x // x ≠ y} → ℝ} (hN : N ∈ simplex {x // x ≠ y}) :
    extend y N ∈ simplex X := by
  refine ⟨fun x => ?_, ?_⟩
  · unfold extend
    split_ifs
    · exact le_rfl
    · exact hN.1 _
  · rw [Fintype.sum_eq_add_sum_subtype_ne _ y, extend_self, zero_add, ← hN.2]
    exact Finset.sum_congr rfl fun x _ => by simp [extend, x.2]

/-! ### Some optimal weighting gives `y` no weight -/

variable [Nonempty n]

omit [DecidableEq n] [Nonempty n] in
/-- No state's greatest density over the other places exceeds its greatest density over all of them. -/
theorem peak_others_le {y : X} [Nonempty {x // x ≠ y}] (D : X → Matrix n n ℂ) (P : Matrix n n ℂ) (μ : ℝ)
    (ρ : Matrix n n ℂ) : peak (others D y) ρ ≤ peak (deepened D P y μ) ρ :=
  ciSup_le fun x => by
    have h := expectation_le_peak (deepened D P y μ) ρ x
    rwa [expectation_deepened_ne ρ D P μ x.2] at h

/-- If a state of `S'` holds the deepened `y` at or below `e'`, it is an optimal state, the common value is `e'`,
and a weighting optimal for the other places, extended by no weight on `y`, is an optimal weighting. -/
theorem dropout_of_held {y : X} [Nonempty {x // x ≠ y}] {D : X → Matrix n n ℂ} (hD : ∀ x, (D x).IsHermitian)
    {P : Matrix n n ℂ} (hP : P.IsHermitian) (μ : ℝ) {e' : ℝ} (he' : IsLeast (peak (others D y) '' states n) e')
    {ρ : Matrix n n ℂ} (hρ : ρ ∈ held D y e') (hρy : expectation ρ (D y) - μ * expectation ρ P ≤ e') :
    IsMinOn (peak (deepened D P y μ)) (states n) ρ ∧ peak (deepened D P y μ) ρ = e' ∧
      ∃ N ∈ simplex X, IsMaxOn (floor (deepened_isHermitian hD hP y μ)) (simplex X) N ∧ N y = 0 ∧
        floor (deepened_isHermitian hD hP y μ) N = e' := by
  have : Nonempty X := ⟨y⟩
  have hpeak : peak (deepened D P y μ) ρ ≤ e' := by
    refine ciSup_le fun x => ?_
    by_cases hx : x = y
    · subst hx
      rwa [expectation_deepened_self]
    · rw [expectation_deepened_ne ρ D P μ hx]
      exact hρ.2 x hx
  -- the full value is at least `e'`
  have hlow : ∀ σ ∈ states n, e' ≤ peak (deepened D P y μ) σ := fun σ hσ =>
    (he'.2 ⟨σ, hσ, rfl⟩).trans (peak_others_le D P μ σ)
  have heq : peak (deepened D P y μ) ρ = e' := le_antisymm hpeak (hlow ρ hρ.1)
  refine ⟨isMinOn_iff.mpr fun σ hσ => heq.le.trans (hlow σ hσ), heq, ?_⟩
  -- a weighting optimal for the other places, extended by no weight on `y`
  obtain ⟨ρ', hρ', N', hN', hmin', -, hval'⟩ := saddle_value (others_isHermitian hD y)
  obtain ⟨σ, hσ, hσe⟩ := he'.1
  have hfloor' : floor (others_isHermitian hD y) N' = e' := by
    rw [← hval']
    exact le_antisymm ((isMinOn_iff.mp hmin' σ hσ).trans hσe.le) (he'.2 ⟨ρ', hρ', rfl⟩)
  have hweighted : weighted (deepened D P y μ) (extend y N') = weighted (others D y) N' := by
    rw [weighted_deepened_of_zero D P y μ (extend_self y N'), restrict_extend]
  have hfloor : floor (deepened_isHermitian hD hP y μ) (extend y N') = e' := by
    rw [floor, lowest_congr _ (weighted_isHermitian (others_isHermitian hD y) N') hweighted]
    exact hfloor'
  refine ⟨extend y N', extend_mem_simplex hN', isMaxOn_iff.mpr fun M hM => ?_, extend_self y N', hfloor⟩
  rw [hfloor, ← heq]
  exact floor_le_peak _ hρ.1 hM

/-- **Proposition `prop:dropout`, first part.** At depth `μ`, some optimal weighting gives `y` no weight exactly
when some state of `S'` holds the deepened `y` at or below `e'`: `Tr(ρ D_y) - μ Tr(ρ P) ≤ e'`. -/
theorem dropout_iff {y : X} [Nonempty {x // x ≠ y}] {D : X → Matrix n n ℂ} (hD : ∀ x, (D x).IsHermitian)
    {P : Matrix n n ℂ} (hP : P.IsHermitian) (μ : ℝ) {e' : ℝ} (he' : IsLeast (peak (others D y) '' states n) e') :
    (∃ N ∈ simplex X, IsMaxOn (floor (deepened_isHermitian hD hP y μ)) (simplex X) N ∧ N y = 0) ↔
      ∃ ρ ∈ held D y e', expectation ρ (D y) - μ * expectation ρ P ≤ e' := by
  have : Nonempty X := ⟨y⟩
  constructor
  · rintro ⟨N, hN, hmax, hy⟩
    -- its floor is the other places' floor, at most `e'`
    obtain ⟨σ, hσ, hσe⟩ := he'.1
    have hfloor : floor (deepened_isHermitian hD hP y μ) N ≤ e' := by
      rw [floor, lowest_congr _ (weighted_isHermitian (others_isHermitian hD y) _)
        (weighted_deepened_of_zero D P y μ hy)]
      exact (floor_le_peak (others_isHermitian hD y) hσ (restrict_mem_simplex hN hy)).trans hσe.le
    -- so an optimal state holds every place at or below `e'`
    obtain ⟨ρ, hρ, N₀, hN₀, -, -, hval⟩ := saddle_value (deepened_isHermitian hD hP y μ)
    have hpeak : peak (deepened D P y μ) ρ ≤ e' := by
      rw [hval]
      exact (isMaxOn_iff.mp hmax N₀ hN₀).trans hfloor
    refine ⟨ρ, ⟨hρ, fun x hx => ?_⟩, ?_⟩
    · rw [← expectation_deepened_ne ρ D P μ hx]
      exact (expectation_le_peak _ ρ x).trans hpeak
    · rw [← expectation_deepened_self]
      exact (expectation_le_peak _ ρ y).trans hpeak
  · rintro ⟨ρ, hρ, hρy⟩
    obtain ⟨-, -, N, hN, hmax, hy, -⟩ := dropout_of_held hD hP μ he' hρ hρy
    exact ⟨N, hN, hmax, hy⟩

/-- **Proposition `prop:dropout`, the common value.** When some optimal weighting gives `y` no weight, the common
value is `e'`: every optimal weighting's floor is `e'`. -/
theorem dropout_value {y : X} [Nonempty {x // x ≠ y}] {D : X → Matrix n n ℂ} (hD : ∀ x, (D x).IsHermitian)
    {P : Matrix n n ℂ} (hP : P.IsHermitian) (μ : ℝ) {e' : ℝ} (he' : IsLeast (peak (others D y) '' states n) e')
    (h : ∃ N ∈ simplex X, IsMaxOn (floor (deepened_isHermitian hD hP y μ)) (simplex X) N ∧ N y = 0)
    {N : X → ℝ} (hN : N ∈ simplex X) (hmax : IsMaxOn (floor (deepened_isHermitian hD hP y μ)) (simplex X) N) :
    floor (deepened_isHermitian hD hP y μ) N = e' := by
  obtain ⟨ρ, hρ, hρy⟩ := (dropout_iff hD hP μ he').mp h
  obtain ⟨-, -, N₀, hN₀, hmax₀, -, hval⟩ := dropout_of_held hD hP μ he' hρ hρy
  rw [← hval]
  exact le_antisymm (isMaxOn_iff.mp hmax₀ N hN) (isMaxOn_iff.mp hmax N₀ hN₀)

/-! ### The depth -/

omit [Nonempty n] [Fintype X] [DecidableEq X] in
/-- `S'` is compact: the states, cut by finitely many closed half-spaces. -/
theorem held_isCompact (D : X → Matrix n n ℂ) (y : X) (e : ℝ) : IsCompact (held D y e) := by
  have hclosed : IsClosed {ρ : Matrix n n ℂ | ∀ x, x ≠ y → expectation ρ (D x) ≤ e} := by
    simp only [Set.ofPred_forall]
    exact isClosed_iInter fun x => isClosed_iInter fun _ =>
      isClosed_le (continuous_expectation (D x)) continuous_const
  exact states_isCompact.inter_right hclosed

/-- At depth zero, if every optimal weighting gives `y` weight, every state of `S'` holds `y` above `e'`. -/
theorem held_above {y : X} [Nonempty {x // x ≠ y}] {D : X → Matrix n n ℂ} (hD : ∀ x, (D x).IsHermitian)
    {P : Matrix n n ℂ} (hP : P.IsHermitian) {e' : ℝ} (he' : IsLeast (peak (others D y) '' states n) e')
    (h0 : ∀ N ∈ simplex X, IsMaxOn (floor (deepened_isHermitian hD hP y 0)) (simplex X) N → 0 < N y) :
    ∀ ρ ∈ held D y e', e' < expectation ρ (D y) := by
  intro ρ hρ
  by_contra hle
  push Not at hle
  obtain ⟨N, hN, hmax, hy⟩ := (dropout_iff hD hP 0 he').mpr ⟨ρ, hρ, by simpa using hle⟩
  exact (h0 N hN hmax).ne' hy

/-- **Proposition `prop:dropout`, the depth.** If every optimal weighting gives `y` weight at depth zero and some
state of `S'` has `Tr(ρ P) > 0`, the least ratio `μ_c` of Eq. `eq:dropout` is attained and positive; at a depth
`μ ≥ 0`, some optimal weighting gives `y` no weight exactly when `μ ≥ μ_c`; and beyond `μ_c` every optimal weighting
gives `y` none, with the common value `e'`. -/
theorem dropout_depth {y : X} [Nonempty {x // x ≠ y}] {D : X → Matrix n n ℂ} (hD : ∀ x, (D x).IsHermitian)
    {P : Matrix n n ℂ} (hP : P.IsHermitian) {e' : ℝ} (he' : IsLeast (peak (others D y) '' states n) e')
    (h0 : ∀ N ∈ simplex X, IsMaxOn (floor (deepened_isHermitian hD hP y 0)) (simplex X) N → 0 < N y)
    (hsome : ∃ ρ ∈ held D y e', 0 < expectation ρ P) :
    ∃ μc : ℝ, 0 < μc ∧
      IsLeast {r | ∃ ρ ∈ held D y e', 0 < expectation ρ P ∧ r = (expectation ρ (D y) - e') / expectation ρ P} μc ∧
      (∀ μ, 0 ≤ μ →
        ((∃ N ∈ simplex X, IsMaxOn (floor (deepened_isHermitian hD hP y μ)) (simplex X) N ∧ N y = 0) ↔ μc ≤ μ)) ∧
      (∀ μ, μc < μ → ∀ N ∈ simplex X, IsMaxOn (floor (deepened_isHermitian hD hP y μ)) (simplex X) N →
        N y = 0 ∧ floor (deepened_isHermitian hD hP y μ) N = e') := by
  have : Nonempty X := ⟨y⟩
  have hpos := held_above hD hP he' h0
  set ratio : Matrix n n ℂ → ℝ := fun ρ => (expectation ρ (D y) - e') / expectation ρ P
  -- the least excess of `y` over `e'` on `S'`, positive
  obtain ⟨ρ₀, hρ₀, hP₀⟩ := hsome
  obtain ⟨ρa, hρa, hmina⟩ := (held_isCompact D y e').exists_isMinOn ⟨ρ₀, hρ₀⟩
    ((continuous_expectation (D y)).sub continuous_const).continuousOn
  set a := expectation ρa (D y) - e'
  have ha : 0 < a := sub_pos.mpr (hpos ρa hρa)
  have hexcess : ∀ ρ ∈ held D y e', a ≤ expectation ρ (D y) - e' := fun ρ hρ => isMinOn_iff.mp hmina ρ hρ
  set r₀ := ratio ρ₀
  have hr₀ : 0 < r₀ := div_pos (sub_pos.mpr (hpos ρ₀ hρ₀)) hP₀
  -- the ratio is least on the compact part of `S'` where `Tr(ρ P) ≥ a / r₀`, and larger off it
  set K₁ := held D y e' ∩ {ρ | a / r₀ ≤ expectation ρ P}
  have hK₁ : IsCompact K₁ :=
    (held_isCompact D y e').inter_right (isClosed_le continuous_const (continuous_expectation P))
  have hρ₀K : ρ₀ ∈ K₁ := by
    refine ⟨hρ₀, ?_⟩
    show a / r₀ ≤ expectation ρ₀ P
    rw [div_le_iff₀ hr₀]
    have hcancel : expectation ρ₀ P * r₀ = expectation ρ₀ (D y) - e' := by
      simp only [r₀, ratio]
      field_simp
    rw [hcancel]
    exact hexcess ρ₀ hρ₀
  have hden : ∀ ρ ∈ K₁, expectation ρ P ≠ 0 := fun ρ hρ => (lt_of_lt_of_le (div_pos ha hr₀) hρ.2).ne'
  obtain ⟨ρc, hρc, hminc⟩ := hK₁.exists_isMinOn ⟨ρ₀, hρ₀K⟩
    ((((continuous_expectation (D y)).sub continuous_const).continuousOn).div
      (continuous_expectation P).continuousOn hden)
  set μc := ratio ρc
  have hPc : 0 < expectation ρc P := lt_of_lt_of_le (div_pos ha hr₀) hρc.2
  have hμc : μc * expectation ρc P = expectation ρc (D y) - e' := div_mul_cancel₀ _ hPc.ne'
  have hleast : ∀ ρ ∈ held D y e', 0 < expectation ρ P → μc ≤ ratio ρ := by
    intro ρ hρ hPρ
    by_cases hin : a / r₀ ≤ expectation ρ P
    · exact isMinOn_iff.mp hminc ρ ⟨hρ, hin⟩
    · push Not at hin
      calc μc ≤ r₀ := isMinOn_iff.mp hminc ρ₀ hρ₀K
        _ = a / (a / r₀) := by field_simp
        _ ≤ a / expectation ρ P := div_le_div_of_nonneg_left ha.le hPρ hin.le
        _ ≤ ratio ρ := div_le_div_of_nonneg_right (hexcess ρ hρ) hPρ.le
  refine ⟨μc, div_pos (sub_pos.mpr (hpos ρc hρc.1)) hPc, ⟨⟨ρc, hρc.1, hPc, rfl⟩, ?_⟩, fun μ hμ => ?_,
    fun μ hμ N hN hmax => ?_⟩
  · rintro _ ⟨ρ, hρ, hPρ, rfl⟩
    exact hleast ρ hρ hPρ
  · rw [dropout_iff hD hP μ he']
    constructor
    · rintro ⟨ρ, hρ, hρy⟩
      have hex := hpos ρ hρ
      have hPρ : 0 < expectation ρ P := by
        by_contra hneg
        push Not at hneg
        nlinarith
      refine (hleast ρ hρ hPρ).trans ?_
      show (expectation ρ (D y) - e') / expectation ρ P ≤ μ
      rw [div_le_iff₀ hPρ]
      linarith
    · intro hle
      exact ⟨ρc, hρc.1, by linarith [mul_le_mul_of_nonneg_right hle hPc.le]⟩
  · -- beyond `μ_c`, the state `ρc` is optimal and holds the deepened `y` strictly below `e'`
    have hbelow : expectation ρc (D y) - μ * expectation ρc P < e' := by
      linarith [mul_lt_mul_of_pos_right hμ hPc]
    obtain ⟨hminρ, hpeakρ, -⟩ := dropout_of_held hD hP μ he' hρc.1 hbelow.le
    have hy : N y = 0 := (saddle_forbids _ hρc.1.1 hN hminρ hmax y).1 (by
      rw [expectation_deepened_self, hpeakρ]
      exact hbelow)
    exact ⟨hy, dropout_value hD hP μ he' ⟨N, hN, hmax, hy⟩ hN hmax⟩

/-- **Proposition `prop:dropout`, never.** If every optimal weighting gives `y` weight at depth zero and no state of
`S'` has `Tr(ρ P) > 0`, every optimal weighting gives `y` weight at every depth `μ ≥ 0`: `y` never drops out. -/
theorem dropout_never {y : X} [Nonempty {x // x ≠ y}] {D : X → Matrix n n ℂ} (hD : ∀ x, (D x).IsHermitian)
    {P : Matrix n n ℂ} (hP : P.IsHermitian) {e' : ℝ} (he' : IsLeast (peak (others D y) '' states n) e')
    (h0 : ∀ N ∈ simplex X, IsMaxOn (floor (deepened_isHermitian hD hP y 0)) (simplex X) N → 0 < N y)
    (hnone : ∀ ρ ∈ held D y e', expectation ρ P ≤ 0) {μ : ℝ} (hμ : 0 ≤ μ) :
    ∀ N ∈ simplex X, IsMaxOn (floor (deepened_isHermitian hD hP y μ)) (simplex X) N → 0 < N y := by
  intro N hN hmax
  refine lt_of_le_of_ne (hN.1 y) fun hy => ?_
  obtain ⟨ρ, hρ, hρy⟩ := (dropout_iff hD hP μ he').mp ⟨N, hN, hmax, hy.symm⟩
  have hex := held_above hD hP he' h0 ρ hρ
  nlinarith [hnone ρ hρ]

end Saddle
