import Saddle.Evolution

/-!
# What the weighted evolution leaves out

Corollary `cor:crossing` of the appendix on the saddle. In the weighted sum of Eq. `eq:evolution`, a link enters at
the mean weight of its two ends and a place's own piece at its own weight. So a link whose two places both carry no
weight has no term, and a place with no weight whose links all join places with none appears in no term at all: the
weighted sum is the sum of the terms that do not involve it.
-/

open Matrix

namespace Saddle

variable {n X : Type*} [Fintype X]

/-- **Corollary `cor:crossing`.** If a place `k` carries no weight and every link at it joins places that carry
none, the weighted sum of Eq. `eq:evolution` is the sum of the terms not involving `k`: the own pieces of the
other places, at their weights, and the links not at `k`, at the mean weight of their ends. -/
theorem weighted_avoids_weightless [DecidableEq X] (own : X → Matrix n n ℂ) (links : Finset (X × X))
    (link : X × X → Matrix n n ℂ) (hlinks : ∀ l ∈ links, l.1 ≠ l.2) (N : X → ℝ) {k : X} (hk : N k = 0)
    (hnbr : ∀ l ∈ links, (l.1 = k ∨ l.2 = k) → N l.1 = 0 ∧ N l.2 = 0) :
    weighted (fun x => own x + (1 / 2 : ℂ) • ∑ l ∈ links.filter (fun l => l.1 = x ∨ l.2 = x), link l) N =
      ∑ x ∈ Finset.univ.filter (fun x => x ≠ k), (N x : ℂ) • own x +
        ∑ l ∈ links.filter (fun l => ¬(l.1 = k ∨ l.2 = k)), (((N l.1 + N l.2) / 2 : ℝ) : ℂ) • link l := by
  rw [weighted_of_pieces own links link hlinks N]
  congr 1
  · rw [← Finset.sum_filter_add_sum_filter_not Finset.univ (fun x => x ≠ k)]
    have hrest : ∑ x ∈ Finset.univ.filter (fun x => ¬x ≠ k), (N x : ℂ) • own x = 0 := by
      refine Finset.sum_eq_zero fun x hx => ?_
      have hxk : x = k := by simpa using (Finset.mem_filter.mp hx).2
      rw [hxk, hk, Complex.ofReal_zero, zero_smul]
    rw [hrest, add_zero]
  · rw [← Finset.sum_filter_add_sum_filter_not links (fun l => l.1 = k ∨ l.2 = k)]
    have hat : ∑ l ∈ links.filter (fun l => l.1 = k ∨ l.2 = k), (((N l.1 + N l.2) / 2 : ℝ) : ℂ) • link l = 0 := by
      refine Finset.sum_eq_zero fun l hl => ?_
      obtain ⟨hl₁, hl₂⟩ := hnbr l (Finset.mem_filter.mp hl).1 (Finset.mem_filter.mp hl).2
      rw [hl₁, hl₂, add_zero, zero_div, Complex.ofReal_zero, zero_smul]
    rw [hat, zero_add]

end Saddle
