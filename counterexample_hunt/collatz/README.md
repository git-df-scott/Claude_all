# Collatz — cycle lower bound + frontier probe

A Collatz counterexample is one of two things: a **divergent trajectory** (not
finitely certifiable — no computation can ever confirm divergence) or a
**nontrivial cycle** (a perfectly finite, instantly checkable certificate).
This module quantifies exactly how large that second object must be.

## `cycle_bound.py` — rigorous cycle-size lower bound

Self-contained derivation (see the module docstring for the mathematics):
combines the cycle product identity `2^K/3^L = Π(1 + 1/(3n_i))` with the law
of best approximation for the continued fraction of log₂3, computed in exact
fixed-point integer interval arithmetic (scale 10⁸⁰) — no floating point
enters any inequality. Convergents with denominator ≤ 10⁶ are independently
checked by exact big-integer comparison of 2^p vs 3^q, and the computed
continued fraction `[1; 1, 1, 2, 2, 3, 1, 5, 2, 23, ...]` matches the known
expansion of log₂3.

Result (actual output in `cycle_bound_output.log`):

| assumed verified frontier | odd members L | total standard-map length |
|---|---|---|
| 2⁷¹ (Barina 2025) | **≥ 65,470,613,321** | **≥ 169,239,080,335** |
| 2⁶⁸ (Barina 2020) | ≥ 6,586,818,670 | ≥ 17,026,679,261 |

For comparison, the best published bound via finer m-cycle theory is
L ≥ 1.375×10¹¹ odd members (Hercher, *J. Integer Seq.* 26, 2023, combined
with the 2⁷¹ verification); our self-contained bound lands within a factor
~2 of it. Either way the conclusion is the same: **the smallest possible
counterexample cycle has tens of billions of elements**, every one of them
exceeding 2⁷¹. No conceivable search stumbles onto that object; it would
have to be constructed from structure nobody has found.

## `frontier_sweep.c` — probe above the verified frontier

Verifies convergence (drop below 2⁷¹, i.e. into exhaustively verified
territory) for the first 10⁸ odd numbers above 2⁷¹, in ~3 s
(`frontier_sweep_output.log`). Worst case in the window: 623 steps, with a
maximum excursion of 2.0×10²⁹ — nine orders of magnitude above the start,
which is why the drift heuristic (E[log-shrink] = log 3 − 2 log 2 < 0) is a
statement about averages, not a safety proof. This probe is engine
demonstration, not progress: extending the actual frontier by one binary
order takes GPU-years (Barina's project), and the cycle bound above shows
what the frontier buys.
