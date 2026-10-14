import Saddle.Forbids
import Mathlib.Analysis.Calculus.Deriv.Slope

/-!
# The weights as rates

Proposition `prop:rates` of the appendix on rates. Shift one place's density by a constant, `D_y → D_y + ε I`. Every
floor moves by `ε N_y`, so the common value becomes `e*(ε) = max_N (f(N) + ε N_y)`, and its one-sided derivatives at
zero are the greatest and the least weight `y` takes over the optimal weightings (Eq. `eq:rates`); where that weight
is unique, `e*` is differentiable with that derivative, and by Corollary `cor:forbids` that derivative is zero when
some optimal state leaves `y` below `e*`. The proof is the paper's, in the general form of a maximum of
`f + ε φ` over the simplex: the optimal weighting with the greatest `φ` bounds the difference quotient below, the
shifted problem's maximizer bounds it above, and that maximizer's `φ` cannot stay above the optimal set's greatest,
since the floor is strictly below `e*` on the compact set where it would.
-/

open Matrix Filter Topology Set
open scoped ComplexOrder

namespace Saddle

variable {n X : Type*} [Fintype n] [DecidableEq n] [Fintype X]

/-- The common value `e*`: the greatest floor over the simplex. -/
noncomputable def commonValue {D : X → Matrix n n ℂ} (hD : ∀ x, (D x).IsHermitian) : ℝ :=
  sSup (floor hD '' simplex X)

/-- The optimal weightings, those whose floor is greatest over the simplex. -/
def optimal {D : X → Matrix n n ℂ} (hD : ∀ x, (D x).IsHermitian) : Set (X → ℝ) :=
  {N | N ∈ simplex X ∧ IsMaxOn (floor hD) (simplex X) N}

omit [DecidableEq n] in
/-- A weight is at most one. -/
theorem weight_le_one {N : X → ℝ} (hN : N ∈ simplex X) (y : X) : N y ≤ 1 :=
  hN.2 ▸ Finset.single_le_sum (fun x _ => hN.1 x) (Finset.mem_univ y)

omit [Fintype X] in
/-- On a compact set, a continuous function's greatest value is the supremum of its image. -/
theorem sSup_image_eq_of_isMaxOn {K : Set (X → ℝ)} {g : (X → ℝ) → ℝ} {N : X → ℝ} (hN : N ∈ K)
    (hmax : IsMaxOn g K N) : sSup (g '' K) = g N :=
  IsGreatest.csSup_eq ⟨mem_image_of_mem g hN, by rintro _ ⟨M, hM, rfl⟩; exact isMaxOn_iff.mp hmax M hM⟩

/-- **The rate of a maximum (Danskin).** For a continuous `φ` bounded by one on the simplex, the right derivative at
zero of `ε ↦ max_N (f(N) + ε φ(N))` is the greatest `φ` over the optimal weightings. -/
theorem danskin_right [Nonempty n] [Nonempty X] {D : X → Matrix n n ℂ} (hD : ∀ x, (D x).IsHermitian)
    (φ : (X → ℝ) → ℝ) (hφ : Continuous φ) (hbound : ∀ N ∈ simplex X, |φ N| ≤ 1) :
    Tendsto (fun ε => (sSup ((fun N => floor hD N + ε * φ N) '' simplex X) - commonValue hD) / ε)
      (𝓝[>] 0) (𝓝 (sSup (φ '' optimal hD))) := by
  set K := simplex X
  set f := floor hD
  have hK : IsCompact K := simplex_isCompact
  have hf : Continuous f := continuous_floor hD
  -- the common value, attained
  obtain ⟨N₀, hN₀, hmax₀⟩ := hK.exists_isMaxOn simplex_nonempty hf.continuousOn
  have he : commonValue hD = f N₀ := sSup_image_eq_of_isMaxOn hN₀ hmax₀
  -- the optimal set is compact and nonempty, and `φ` takes its greatest value there
  have hO : optimal hD = K ∩ ⋂ M ∈ K, {N | f M ≤ f N} := by
    ext N
    simp only [optimal, mem_ofPred_eq, mem_inter_iff, mem_iInter, isMaxOn_iff]
    rfl
  have hOc : IsCompact (optimal hD) := by
    rw [hO]
    exact hK.inter_right (isClosed_biInter fun M _ => isClosed_le continuous_const hf)
  obtain ⟨N₁, hN₁, hmax₁⟩ := hOc.exists_isMaxOn ⟨N₀, hN₀, hmax₀⟩ hφ.continuousOn
  set top := sSup (φ '' optimal hD)
  have htop : top = φ N₁ := sSup_image_eq_of_isMaxOn hN₁ hmax₁
  have hN₁e : f N₁ = f N₀ := le_antisymm (isMaxOn_iff.mp hmax₀ N₁ hN₁.1) (isMaxOn_iff.mp hN₁.2 N₀ hN₀)
  -- for each ε, the shifted problem's maximizer, and the two bounds on the difference quotient
  have hshift : ∀ ε : ℝ, ∃ N ∈ K, sSup ((fun N => f N + ε * φ N) '' K) = f N + ε * φ N ∧
      ∀ M ∈ K, f M + ε * φ M ≤ f N + ε * φ N := fun ε => by
    obtain ⟨N, hN, hmax⟩ := hK.exists_isMaxOn simplex_nonempty (hf.add (continuous_const.mul hφ)).continuousOn
    exact ⟨N, hN, sSup_image_eq_of_isMaxOn hN hmax, isMaxOn_iff.mp hmax⟩
  have hlow : ∀ ε > 0, top ≤ (sSup ((fun N => f N + ε * φ N) '' K) - commonValue hD) / ε := fun ε hε => by
    obtain ⟨N, -, hval, hbest⟩ := hshift ε
    rw [le_div_iff₀ hε, hval, he, htop]
    linarith [hbest N₁ hN₁.1]
  refine tendsto_order.2 ⟨fun a ha => ?_, fun b hb => ?_⟩
  · exact eventually_nhdsWithin_of_forall fun ε hε => ha.trans_le (hlow ε hε)
  -- above: where `φ` exceeds the optimal set's greatest by `δ`, the floor stays below `e*` by some `η`
  set δ := (b - top) / 2
  have hδ : 0 < δ := by simp only [δ]; linarith
  set far := K ∩ {N | top + δ ≤ φ N}
  have hfarc : IsCompact far := hK.inter_right (isClosed_le continuous_const hφ)
  have hgap : ∃ η > 0, ∀ N ∈ far, f N ≤ f N₀ - η := by
    rcases far.eq_empty_or_nonempty with hempty | hne
    · exact ⟨1, one_pos, by simp [hempty]⟩
    obtain ⟨N₂, hN₂, hmax₂⟩ := hfarc.exists_isMaxOn hne hf.continuousOn
    have hbelow : f N₂ < f N₀ := by
      refine lt_of_le_of_ne (isMaxOn_iff.mp hmax₀ N₂ hN₂.1) fun heq => ?_
      have hopt : N₂ ∈ optimal hD :=
        ⟨hN₂.1, isMaxOn_iff.mpr fun M hM => (isMaxOn_iff.mp hmax₀ M hM).trans heq.symm.le⟩
      have hle := isMaxOn_iff.mp hmax₁ N₂ hopt
      have hfar₂ : top + δ ≤ φ N₂ := hN₂.2
      linarith [htop, hδ]
    exact ⟨f N₀ - f N₂, sub_pos.mpr hbelow, fun N hN => by linarith [isMaxOn_iff.mp hmax₂ N hN]⟩
  obtain ⟨η, hη, hfar⟩ := hgap
  have hsmall : Ioo (0 : ℝ) (η / 3) ∈ 𝓝[>] (0 : ℝ) := Ioo_mem_nhdsGT (by positivity)
  filter_upwards [hsmall] with ε hε
  obtain ⟨N, hN, hval, hbest⟩ := hshift ε
  have hφN : φ N < top + δ := by
    by_contra hnot
    have hfN := hfar N ⟨hN, not_lt.mp hnot⟩
    have h₁ := hbest N₁ hN₁.1
    have h₂ := abs_le.mp (hbound N hN)
    have h₃ := abs_le.mp (hbound N₁ hN₁.1)
    nlinarith [hε.1, hε.2]
  have hup : (sSup ((fun N => f N + ε * φ N) '' K) - commonValue hD) / ε ≤ φ N := by
    rw [div_le_iff₀ hε.1, hval, he]
    linarith [isMaxOn_iff.mp hmax₀ N hN]
  simp only [δ] at hφN
  linarith

variable [DecidableEq X]

/-- The densities with place `y`'s shifted by `ε`, `D_y + ε I`, every other place's unchanged. -/
noncomputable def shifted (D : X → Matrix n n ℂ) (y : X) (ε : ℝ) : X → Matrix n n ℂ :=
  fun x => if x = y then D x + (ε : ℂ) • 1 else D x

omit [Fintype n] [Fintype X] in
/-- The shifted densities are Hermitian. -/
theorem shifted_isHermitian {D : X → Matrix n n ℂ} (hD : ∀ x, (D x).IsHermitian) (y : X) (ε : ℝ) :
    ∀ x, (shifted D y ε x).IsHermitian := fun x => by
  unfold shifted
  split_ifs
  · unfold Matrix.IsHermitian
    rw [conjTranspose_add, (hD x).eq, conjTranspose_smul, conjTranspose_one, Complex.star_def, Complex.conj_ofReal]
  · exact hD x

omit [Fintype n] in
/-- Shifting place `y`'s density by `ε` shifts every weighted sum by `ε N_y`. -/
theorem weighted_shifted (D : X → Matrix n n ℂ) (y : X) (ε : ℝ) (N : X → ℝ) :
    weighted (shifted D y ε) N = weighted D N + ((ε * N y : ℝ) : ℂ) • 1 := by
  have hterm : ∀ x, (N x : ℂ) • shifted D y ε x =
      (N x : ℂ) • D x + if x = y then ((ε * N y : ℝ) : ℂ) • (1 : Matrix n n ℂ) else 0 := fun x => by
    unfold shifted
    split_ifs with h
    · subst h
      rw [smul_add, smul_smul]
      congr 2
      push_cast
      ring
    · simp
  simp only [weighted, hterm, Finset.sum_add_distrib, Finset.sum_ite_eq', Finset.mem_univ, ite_true]

omit [Fintype X] [DecidableEq X] in
/-- Adding `c I` to a Hermitian matrix adds `c` to its lowest eigenvalue. -/
theorem lowest_add_smul_one [Nonempty n] {A : Matrix n n ℂ} (hA : A.IsHermitian) (c : ℝ)
    (hAc : (A + (c : ℂ) • 1).IsHermitian) : lowest hAc = lowest hA + c := by
  have hstate : ∀ ρ ∈ states n, expectation ρ (A + (c : ℂ) • 1) = expectation ρ A + c := fun ρ hρ => by
    simp only [expectation, Matrix.mul_add, trace_add, Complex.add_re, Matrix.mul_smul, trace_smul,
      hρ.2, smul_eq_mul, mul_one, Complex.ofReal_re]
  refine (lowest_isLeast hAc).unique ⟨?_, ?_⟩
  · obtain ⟨⟨σ, hσ, hσeq⟩, -⟩ := lowest_isLeast hA
    exact ⟨σ, hσ, by simp only at hσeq ⊢; rw [hstate σ hσ, hσeq]⟩
  · rintro _ ⟨ρ, hρ, rfl⟩
    simp only
    rw [hstate ρ hρ]
    linarith [(lowest_isLeast hA).2 ⟨ρ, hρ, rfl⟩]

/-- Shifting place `y`'s density by `ε` moves the floor of every weighting by `ε N_y`. -/
theorem floor_shifted [Nonempty n] {D : X → Matrix n n ℂ} (hD : ∀ x, (D x).IsHermitian) (y : X) (ε : ℝ) :
    floor (shifted_isHermitian hD y ε) = fun N => floor hD N + ε * N y := by
  funext N
  have h := weighted_shifted D y ε N
  have hherm : (weighted D N + ((ε * N y : ℝ) : ℂ) • 1).IsHermitian :=
    h ▸ weighted_isHermitian (shifted_isHermitian hD y ε) N
  rw [floor, lowest_congr _ hherm h, lowest_add_smul_one (weighted_isHermitian hD N), floor]

/-- The common value of the shifted densities is the greatest of the shifted floors. -/
theorem commonValue_shifted [Nonempty n] {D : X → Matrix n n ℂ} (hD : ∀ x, (D x).IsHermitian) (y : X) (ε : ℝ) :
    commonValue (shifted_isHermitian hD y ε) = sSup ((fun N => floor hD N + ε * N y) '' simplex X) := by
  unfold commonValue
  rw [floor_shifted]

omit [DecidableEq n] [DecidableEq X] in
/-- The weights `y` takes over the optimal weightings are bounded by one in size. -/
theorem abs_weight_le_one {N : X → ℝ} (hN : N ∈ simplex X) (y : X) : |N y| ≤ 1 :=
  abs_le.mpr ⟨by linarith [hN.1 y], weight_le_one hN y⟩

/-- **Proposition `prop:rates`, the right derivative.** Shifting place `y`'s density by `ε > 0` moves the common value
at the rate of the greatest weight `y` takes over the optimal weightings (Eq. `eq:rates`). -/
theorem saddle_rates_right [Nonempty n] [Nonempty X] {D : X → Matrix n n ℂ} (hD : ∀ x, (D x).IsHermitian) (y : X) :
    Tendsto (fun ε => (commonValue (shifted_isHermitian hD y ε) - commonValue hD) / ε) (𝓝[>] 0)
      (𝓝 (sSup ((fun N => N y) '' optimal hD))) := by
  refine (danskin_right hD (fun N => N y) (continuous_apply y) fun N hN => abs_weight_le_one hN y).congr
    fun ε => ?_
  rw [commonValue_shifted hD y ε]

/-- **Proposition `prop:rates`, the left derivative.** Shifting place `y`'s density by `ε < 0` moves the common value
at the rate of the least weight `y` takes over the optimal weightings (Eq. `eq:rates`). -/
theorem saddle_rates_left [Nonempty n] [Nonempty X] {D : X → Matrix n n ℂ} (hD : ∀ x, (D x).IsHermitian) (y : X) :
    Tendsto (fun ε => (commonValue (shifted_isHermitian hD y ε) - commonValue hD) / ε) (𝓝[<] 0)
      (𝓝 (sInf ((fun N => N y) '' optimal hD))) := by
  have hright := danskin_right hD (fun N => -N y) (continuous_apply y).neg fun N hN => by
    rw [abs_neg]; exact abs_weight_le_one hN y
  have hneg : Tendsto (fun ε : ℝ => -ε) (𝓝[<] 0) (𝓝[>] 0) :=
    tendsto_nhdsWithin_of_tendsto_nhds_of_eventually_within _
      ((continuous_neg.tendsto' (0 : ℝ) 0 neg_zero).mono_left nhdsWithin_le_nhds)
      (eventually_nhdsWithin_of_forall fun ε hε => mem_Ioi.mpr (neg_pos.mpr hε))
  have hcomp := hright.neg.comp hneg
  have hinf : -sSup ((fun N => -N y) '' optimal hD) = sInf ((fun N => N y) '' optimal hD) := by
    rw [Real.sInf_def, ← Set.image_neg_eq_neg, Set.image_image]
  rw [hinf] at hcomp
  refine hcomp.congr' (eventually_nhdsWithin_of_forall fun ε hε => ?_)
  simp only [Function.comp]
  rw [commonValue_shifted hD y ε]
  have hne : ε ≠ 0 := (mem_Iio.mp hε).ne
  rw [show (fun N => floor hD N + -ε * -N y) = fun N => floor hD N + ε * N y by funext N; ring]
  field_simp

/-- **Proposition `prop:rates`, where the weight is unique.** If every optimal weighting gives `y` the same weight
`c`, the common value is differentiable in the shift, with derivative `c`. -/
theorem saddle_rates_unique [Nonempty n] [Nonempty X] {D : X → Matrix n n ℂ} (hD : ∀ x, (D x).IsHermitian) (y : X)
    {c : ℝ} (hc : ∀ N ∈ optimal hD, N y = c) :
    HasDerivAt (fun ε => commonValue (shifted_isHermitian hD y ε)) c 0 := by
  obtain ⟨N₀, hN₀, hmax₀⟩ := simplex_isCompact.exists_isMaxOn simplex_nonempty (continuous_floor hD).continuousOn
  have himage : (fun N => N y) '' optimal hD = {c} :=
    Set.eq_singleton_iff_unique_mem.mpr ⟨⟨N₀, ⟨hN₀, hmax₀⟩, hc N₀ ⟨hN₀, hmax₀⟩⟩, by rintro _ ⟨N, hN, rfl⟩; exact hc N hN⟩
  have hzero : commonValue (shifted_isHermitian hD y 0) = commonValue hD := by
    rw [commonValue_shifted hD y 0]
    simp only [zero_mul, add_zero]
    rfl
  have hright := saddle_rates_right hD y
  have hleft := saddle_rates_left hD y
  rw [himage, csSup_singleton] at hright
  rw [himage, csInf_singleton] at hleft
  rw [hasDerivAt_iff_tendsto_slope_zero, ← Set.Iio_union_Ioi, nhdsWithin_union]
  refine (hleft.sup hright).congr' (Eventually.of_forall fun ε => ?_)
  simp only [zero_add, hzero, smul_eq_mul]
  ring

/-- **Proposition `prop:rates`, a place below the common value.** If some optimal state leaves `y` below the common
value, the common value does not move to first order when `y`'s density is shifted: by Corollary `cor:forbids`
every optimal weighting gives `y` no weight. -/
theorem saddle_rates_forbidden [Nonempty n] [Nonempty X] {D : X → Matrix n n ℂ} (hD : ∀ x, (D x).IsHermitian)
    (y : X) {ρ : Matrix n n ℂ} (hρ : ρ ∈ states n) (hmin : IsMinOn (peak D) (states n) ρ)
    (hlow : expectation ρ (D y) < peak D ρ) :
    HasDerivAt (fun ε => commonValue (shifted_isHermitian hD y ε)) 0 0 :=
  saddle_rates_unique hD y fun _ hN => (saddle_forbids hD hρ hN.1 hmin hN.2 y).1 hlow

end Saddle
