import Saddle.Minimax
import Mathlib.LinearAlgebra.Matrix.Kronecker

/-!
# Idle probes

Lemma `lem:probes` of the appendix on probes. Attach to each place `x` an idle probe with its own space and the
same energy `K`, counted in the place's density: `D'ₓ = Dₓ ⊗ I + Kₓ`, where `Kₓ` acts as `K` on probe `x` and as
the identity on the other probes. For every weighting `N`, `λ_min(∑ₓ Nₓ D'ₓ) = λ_min(∑ₓ Nₓ Dₓ) + λ_min(K)`
(Eq. `eq:probes`), so the probes shift every floor, and so the common value, by `λ_min(K)`, and leave every optimal
weighting optimal. The proof is the paper's: the weighted sum is `A_N ⊗ I + ∑ₓ Nₓ Kₓ`; no state takes `A_N ⊗ I`
below `λ_min(A_N)` or any `Kₓ` below `λ_min(K)`, the weights sum to one, and a ground state of `A_N` beside every
probe in a ground state of `K` reaches both at once.
-/

open Matrix
open scoped ComplexOrder Kronecker

namespace Saddle

variable {n X p : Type*} [Fintype n] [DecidableEq n] [Fintype X] [DecidableEq X] [Fintype p] [DecidableEq p]

/-- Probe `x`'s energy on the probes' joint space: `K` on probe `x`, the identity on every other probe. -/
noncomputable def probe (K : Matrix p p ℂ) (x : X) : Matrix (X → p) (X → p) ℂ :=
  (K ⊗ₖ (1 : Matrix ({z // z ≠ x} → p) ({z // z ≠ x} → p) ℂ)).submatrix
    (Equiv.piSplitAt x fun _ => p) (Equiv.piSplitAt x fun _ => p)

/-- The densities with the probes attached: `D'ₓ = Dₓ ⊗ I + Kₓ`. -/
noncomputable def probed (D : X → Matrix n n ℂ) (K : Matrix p p ℂ) (x : X) :
    Matrix (n × (X → p)) (n × (X → p)) ℂ :=
  D x ⊗ₖ 1 + 1 ⊗ₖ probe K x

omit [Fintype p] in
/-- A probe's energy is Hermitian. -/
theorem probe_isHermitian {K : Matrix p p ℂ} (hK : K.IsHermitian) (x : X) : (probe K x).IsHermitian := by
  refine IsHermitian.submatrix ?_ _
  unfold Matrix.IsHermitian
  rw [conjTranspose_kronecker, hK.eq, conjTranspose_one]

omit [Fintype n] [Fintype p] in
/-- The densities with the probes attached are Hermitian. -/
theorem probed_isHermitian {D : X → Matrix n n ℂ} (hD : ∀ x, (D x).IsHermitian) {K : Matrix p p ℂ}
    (hK : K.IsHermitian) (x : X) : (probed D K x).IsHermitian := by
  refine IsHermitian.add ?_ ?_
  · unfold Matrix.IsHermitian
    rw [conjTranspose_kronecker, (hD x).eq, conjTranspose_one]
  · unfold Matrix.IsHermitian
    rw [conjTranspose_kronecker, (probe_isHermitian hK x).eq, conjTranspose_one]

/-! ### No state goes below either part -/

omit [Fintype n] [DecidableEq n] [Fintype X] [DecidableEq X] [Fintype p] [DecidableEq p] in
/-- `(A - c I) ⊗ I = A ⊗ I - c I`. -/
theorem sub_smul_kronecker_one {m q : Type*} [DecidableEq m] [DecidableEq q] (A : Matrix m m ℂ) (c : ℂ) :
    (A - c • 1) ⊗ₖ (1 : Matrix q q ℂ) = A ⊗ₖ 1 - c • ((1 : Matrix m m ℂ) ⊗ₖ (1 : Matrix q q ℂ)) := by
  ext ⟨a, r⟩ ⟨b, s⟩
  simp only [kroneckerMap_apply, Matrix.sub_apply, Matrix.smul_apply, smul_eq_mul]
  ring

omit [Fintype n] [DecidableEq n] [Fintype X] [DecidableEq X] [Fintype p] [DecidableEq p] in
/-- `I ⊗ (B - c I) = I ⊗ B - c I`. -/
theorem one_kronecker_sub_smul {m q : Type*} [DecidableEq m] [DecidableEq q] (B : Matrix q q ℂ) (c : ℂ) :
    (1 : Matrix m m ℂ) ⊗ₖ (B - c • 1) = 1 ⊗ₖ B - c • ((1 : Matrix m m ℂ) ⊗ₖ (1 : Matrix q q ℂ)) := by
  ext ⟨a, r⟩ ⟨b, s⟩
  simp only [kroneckerMap_apply, Matrix.sub_apply, Matrix.smul_apply, smul_eq_mul]
  ring

omit [Fintype n] [DecidableEq n] [Fintype p] in
/-- Shifting `K` by a constant shifts each probe's energy by it. -/
theorem probe_sub_smul_one (K : Matrix p p ℂ) (x : X) (c : ℂ) : probe K x - c • 1 = probe (K - c • 1) x := by
  rw [probe, probe, sub_smul_kronecker_one, one_kronecker_one,
    ← submatrix_one_equiv (α := ℂ) (Equiv.piSplitAt x fun _ => p)]
  rfl

omit [Fintype X] [DecidableEq X] [Fintype p] [DecidableEq p] in
/-- No state gives `A ⊗ I` less than `λ_min(A)`. -/
theorem lowest_le_expectation_kronecker_one {q : Type*} [Fintype q] [DecidableEq q] {A : Matrix n n ℂ}
    (hA : A.IsHermitian) {ρ : Matrix (n × q) (n × q) ℂ} (hρ : ρ ∈ states (n × q)) :
    lowest hA ≤ expectation ρ (A ⊗ₖ (1 : Matrix q q ℂ)) := by
  have h := (trace_mul_re_nonneg hρ.1 ((sub_lowest_posSemidef hA).kronecker (PosSemidef.one (n := q)))).1
  rw [sub_smul_kronecker_one, one_kronecker_one, expectation_sub_smul_one hρ] at h
  linarith

/-- No state gives a probe's energy less than `λ_min(K)`. -/
theorem lowest_le_expectation_probe {K : Matrix p p ℂ} (hK : K.IsHermitian) (x : X)
    {ρ : Matrix (n × (X → p)) (n × (X → p)) ℂ} (hρ : ρ ∈ states (n × (X → p))) :
    lowest hK ≤ expectation ρ (1 ⊗ₖ probe K x) := by
  have hpsd : (probe K x - (lowest hK : ℂ) • 1).PosSemidef := by
    rw [probe_sub_smul_one]
    exact ((sub_lowest_posSemidef hK).kronecker PosSemidef.one).submatrix _
  have h := (trace_mul_re_nonneg hρ.1 ((PosSemidef.one (n := n)).kronecker hpsd)).1
  rw [one_kronecker_sub_smul, one_kronecker_one, expectation_sub_smul_one hρ] at h
  linarith

/-! ### Every probe in its ground state -/

/-- Every probe in the state `b`: `B(s) = ∏ₓ b(sₓ)`. -/
def product (b : p → ℂ) : (X → p) → ℂ :=
  fun s => ∏ z, b (s z)

omit [Fintype p] [DecidableEq p] in
/-- Split off probe `x`: `B = b ⊗ ∏_{z ≠ x} b`. -/
theorem product_split (b : p → ℂ) (x : X) (c : p) (r : {z // z ≠ x} → p) :
    product b ((Equiv.piSplitAt x fun _ => p).symm (c, r)) = b c * ∏ z : {z // z ≠ x}, b (r z) := by
  rw [product, Fintype.prod_eq_mul_prod_subtype_ne _ x]
  congr 1
  · simp [Equiv.piSplitAt]
  · refine Finset.prod_congr rfl fun z _ => ?_
    simp [Equiv.piSplitAt, z.2]

/-- If `b` is an eigenvector of `K`, every probe in `b` is an eigenvector of each probe's energy. -/
theorem probe_mulVec_product {K : Matrix p p ℂ} {b : p → ℂ} {κ : ℂ} (hb : K *ᵥ b = κ • b) (x : X) :
    probe K x *ᵥ product b = κ • product b := by
  have hF : (K ⊗ₖ (1 : Matrix ({z // z ≠ x} → p) ({z // z ≠ x} → p) ℂ)) *ᵥ
      (product b ∘ (Equiv.piSplitAt x fun _ => p).symm) = κ • (product b ∘ (Equiv.piSplitAt x fun _ => p).symm) := by
    ext ⟨a, r⟩
    have hba := congrFun hb a
    simp only [mulVec, dotProduct, Pi.smul_apply, smul_eq_mul] at hba
    simp only [mulVec, dotProduct, Fintype.sum_prod_type, kroneckerMap_apply, Function.comp_apply, product_split,
      one_apply, Pi.smul_apply, smul_eq_mul, mul_ite, mul_one, mul_zero, ite_mul, zero_mul, Finset.sum_ite_eq,
      Finset.mem_univ, ite_true]
    simp only [← mul_assoc]
    rw [← Finset.sum_mul, hba]
  rw [probe, submatrix_mulVec_equiv, hF]
  ext s
  simp only [Function.comp_apply, Pi.smul_apply, Equiv.symm_apply_apply]

omit [DecidableEq p] in
/-- `B† B = ∏ₓ b† b`: every probe in a unit vector is a unit vector. -/
theorem star_product_dotProduct (b : p → ℂ) (hb : star b ⬝ᵥ b = 1) :
    star (product (X := X) b) ⬝ᵥ product b = 1 := by
  simp only [dotProduct, product, Pi.star_apply, star_prod, ← Finset.prod_mul_distrib]
  rw [← Fintype.piFinset_univ, ← Finset.prod_univ_sum (fun _ => Finset.univ) fun _ c => star (b c) * b c]
  simp only [dotProduct, Pi.star_apply] at hb
  have hb' : ∑ c, (starRingEnd ℂ) (b c) * b c = 1 := hb
  simp [hb']

/-! ### The lemma -/

omit [DecidableEq n] in
/-- The value a state gives a sum of operators is the sum of the values. -/
theorem expectation_add_right (ρ A B : Matrix n n ℂ) :
    expectation ρ (A + B) = expectation ρ A + expectation ρ B := by
  simp only [expectation, Matrix.mul_add, trace_add, Complex.add_re]

omit [Fintype n] [DecidableEq n] [DecidableEq X] [Fintype p] [DecidableEq p] in
/-- `A_N ⊗ I = ∑ₓ Nₓ (Dₓ ⊗ I)`. -/
theorem weighted_kronecker_one {q : Type*} [DecidableEq q] (D : X → Matrix n n ℂ) (N : X → ℝ) :
    weighted (fun x => D x ⊗ₖ (1 : Matrix q q ℂ)) N = weighted D N ⊗ₖ 1 := by
  ext ⟨i, r⟩ ⟨j, s⟩
  simp [weighted, Matrix.sum_apply, kroneckerMap_apply, Finset.sum_mul, mul_assoc]

/-- **Lemma `lem:probes`.** For every weighting `N`, identical idle probes raise the floor by the lowest eigenvalue
of their energy: `λ_min(∑ₓ Nₓ D'ₓ) = λ_min(∑ₓ Nₓ Dₓ) + λ_min(K)` (Eq. `eq:probes`). -/
theorem floor_probed [Nonempty n] [Nonempty p] {D : X → Matrix n n ℂ} (hD : ∀ x, (D x).IsHermitian)
    {K : Matrix p p ℂ} (hK : K.IsHermitian) {N : X → ℝ} (hN : N ∈ simplex X) :
    floor (probed_isHermitian hD hK) N = floor hD N + lowest hK := by
  refine (lowest_isLeast (weighted_isHermitian (probed_isHermitian hD hK) N)).unique ⟨?_, ?_⟩
  · -- a ground state of `A_N` beside every probe in a ground state of `K`
    obtain ⟨⟨σ, hσ, hσeq⟩, -⟩ := lowest_isLeast (weighted_isHermitian hD N)
    obtain ⟨i, hi⟩ := exists_eigenvalues_eq_lowest hK
    set b : p → ℂ := ⇑(hK.eigenvectorBasis i)
    have hb : K *ᵥ b = ((lowest hK : ℝ) : ℂ) • b := by
      rw [hK.mulVec_eigenvectorBasis i, hi]
      ext c
      simp [Complex.real_smul, b]
    have hb1 : star b ⬝ᵥ b = 1 := by
      have h := (eigenstate_mem_states hK i).2
      rwa [eigenstate, trace_vecMulVec, dotProduct_comm] at h
    have hB1 := star_product_dotProduct (X := X) b hb1
    set T := vecMulVec (product (X := X) b) (star (product (X := X) b))
    have hT : T ∈ states (X → p) := by
      refine ⟨posSemidef_vecMulVec_self_star _, ?_⟩
      rw [trace_vecMulVec, dotProduct_comm, hB1]
    have hTprobe : ∀ x, expectation T (probe K x) = lowest hK := fun x => by
      rw [expectation_vecMulVec, probe_mulVec_product hb x, dotProduct_smul, hB1, smul_eq_mul, mul_one,
        Complex.ofReal_re]
    have hx : ∀ x, expectation (σ ⊗ₖ T) (probed D K x) = expectation σ (D x) + lowest hK := fun x => by
      rw [probed, expectation_add_right, ← hTprobe x]
      simp only [expectation, ← mul_kronecker_mul, trace_kronecker, hT.2, hσ.2,
        mul_one, one_mul]
    refine ⟨σ ⊗ₖ T, ⟨hσ.1.kronecker hT.1, by rw [trace_kronecker, hσ.2, hT.2, one_mul]⟩, ?_⟩
    simp only
    rw [expectation_weighted]
    simp only [hx, mul_add, Finset.sum_add_distrib, ← Finset.sum_mul, hN.2, one_mul]
    rw [← expectation_weighted]
    exact congrArg (· + lowest hK) hσeq
  · -- no state goes below
    rintro _ ⟨ρ, hρ, rfl⟩
    simp only
    rw [expectation_weighted]
    simp only [probed, expectation_add_right, mul_add, Finset.sum_add_distrib]
    have hA : ∑ x, N x * expectation ρ (D x ⊗ₖ 1) = expectation ρ (weighted D N ⊗ₖ 1) := by
      rw [← weighted_kronecker_one, expectation_weighted]
    have hlowA := lowest_le_expectation_kronecker_one (weighted_isHermitian hD N) hρ
    have hlowK : ∑ x, N x * lowest hK ≤ ∑ x, N x * expectation ρ (1 ⊗ₖ probe K x) :=
      Finset.sum_le_sum fun x _ => mul_le_mul_of_nonneg_left (lowest_le_expectation_probe hK x hρ) (hN.1 x)
    rw [← Finset.sum_mul, hN.2, one_mul] at hlowK
    rw [hA]
    exact add_le_add hlowA hlowK

/-- **Lemma `lem:probes`, its consequence.** Idle probes leave every optimal weighting optimal and shift the common
value by `λ_min(K)`: a weighting certifies the greatest floor with the probes exactly when it does without them,
and its floor is raised by `λ_min(K)`. -/
theorem saddle_probes [Nonempty n] [Nonempty p] {D : X → Matrix n n ℂ} (hD : ∀ x, (D x).IsHermitian)
    {K : Matrix p p ℂ} (hK : K.IsHermitian) {N : X → ℝ} (hN : N ∈ simplex X) :
    (IsMaxOn (floor (probed_isHermitian hD hK)) (simplex X) N ↔ IsMaxOn (floor hD) (simplex X) N) ∧
      floor (probed_isHermitian hD hK) N = floor hD N + lowest hK := by
  refine ⟨?_, floor_probed hD hK hN⟩
  simp only [isMaxOn_iff]
  refine forall₂_congr fun M hM => ?_
  rw [floor_probed hD hK hM, floor_probed hD hK hN, add_le_add_iff_right]

end Saddle
