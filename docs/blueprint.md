<!-- GENERATED from the paper's sources; edit those, not this file. -->
# Blueprint

Each derived claim of [CLAIMS.md](../CLAIMS.md), with the proposition of the paper that proves it and the Lean declaration that proves it again, checked by machine. A declaration is proved when it compiles with no `sorry` and rests on no axiom beyond Lean's standard three; one still being written is owed. `make lean` rebuilds the proofs and the report this table reads.

| Claim | Proved in the paper by | Lean declaration | In Lean |
| --- | --- | --- | --- |
| saddle.exchange | prop:exchange | `saddle_value` (lean/Saddle/Minimax.lean) | proved |
| saddle.conditions | prop:conditions | `saddle_conditions` (lean/Saddle/Conditions.lean) | proved |
| saddle.certificate | cor:gap | `saddle_bracket` (lean/Saddle/Certificate.lean); `saddle_certificate` (lean/Saddle/Certificate.lean) | proved |
| saddle.forbids | cor:forbids | `saddle_forbids` (lean/Saddle/Forbids.lean) | proved |
| saddle.rates | cor:forbids, prop:rates | `saddle_rates_right` (lean/Saddle/Rates.lean); `saddle_rates_left` (lean/Saddle/Rates.lean); `saddle_rates_unique` (lean/Saddle/Rates.lean); `saddle_rates_forbidden` (lean/Saddle/Rates.lean) | proved |
| uniform.exact | prop:uniform | `uniform_exact` (lean/Saddle/Uniform.lean) | proved |
| programs.real | lem:real | `real_states_attain` (lean/Saddle/Real.lean) | proved |
| saddle.evolution | cor:evolution | `saddle_evolution` (lean/Saddle/Evolution.lean); `weighted_of_pieces` (lean/Saddle/Evolution.lean) | proved |
| evolution.conditioned | prop:conditioned | `conditioned_evolution` (lean/Saddle/Conditioned.lean) | proved |
| division.any | prop:division | `every_weighting_some_division` (lean/Saddle/Division.lean) | proved |
| division.canonical | app:division | — | not formalized: The parts are defined by partial traces over the factors of a many-place tensor product, which Mathlib does not provide for matrices; stating the decomposition needs them built first, which the proving phase takes up. |
| division.unique | app:division | — | not formalized: The proof rests on Schur's lemma and on the irreducibility of the traceless operators at a place under conjugation by its unitaries, which Mathlib does not state for products of unitary groups acting on many-place tensor products; the proving phase takes it up with the partial traces of division.canonical. |
| bell.factorization | prop:factorization | `chsh_of_one_description` (lean/Saddle/Bell.lean); `chsh_of_parts` (lean/Saddle/Bell.lean) | proved |
| history.conditioning | app:bell | — | not formalized: The null space is found by conjugating with the steps' products into a path Laplacian on the step register tensored with the identity, an operator identity on a tensor product that Mathlib does not state for matrices indexed by a register and a rest; the proving phase takes it up. |
| signals.wall | cor:crossing | `weighted_avoids_weightless` (lean/Saddle/Crossing.lean) | proved |
| saddle.dropout | prop:dropout | `dropout_iff` (lean/Saddle/Dropout.lean); `dropout_value` (lean/Saddle/Dropout.lean); `dropout_depth` (lean/Saddle/Dropout.lean); `dropout_never` (lean/Saddle/Dropout.lean) | proved |
| parity.rings | app:dropout | — | not formalized: The even ring's law is first-order perturbation theory of a nondegenerate ground state in the field, which Mathlib cannot yet state; the odd ring's case is the drop-out depth, proved in Lean as saddle.dropout. |
| probes.idle | lem:probes | `floor_probed` (lean/Saddle/Probes.lean); `saddle_probes` (lean/Saddle/Probes.lean) | proved |
| equilibrium.coldest | app:coldest | — | not formalized: The state of greatest entropy is found through maximum-entropy duality for density matrices, which Mathlib cannot yet state. |
| saddle.gain | cor:gain | `saddle_gain` (lean/Saddle/Uniform.lean) | proved |
| single.weighting | app:saddle | — | not formalized: The proof uses second-order perturbation theory of a single eigenvalue, analytic in the weighting, which Mathlib does not yet state for matrices. |
| coldest.symmetric | app:coldest | — | not formalized: The proof uses maximum-entropy duality for density matrices, which Mathlib does not yet state, as Proposition prop:coldest does. |
