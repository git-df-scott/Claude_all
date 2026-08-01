# Andrews–Curtis search engine

`ac_search.c` — bounded exhaustive/greedy search over Andrews–Curtis moves
for balanced 2-generator presentations of the trivial group.

## Design

- State: pair of freely reduced words over {x, X, y, Y} (X = x⁻¹),
  canonicalized by sorting the pair (the move set is r1/r2-symmetric).
- 14 elementary moves, all classical AC moves: r_i → r_i·r_j^{±1} (i ≠ j),
  r_i → r_i⁻¹, and r_i → g·r_i·g⁻¹ for each generator letter g.
- Exact visited set (open-addressing hash table with full state comparison —
  no lossy fingerprinting, so "exhausted" claims are rigorous).
- Greedy expansion: always expand the state of smallest total relator length.
- Per-relator length cap; a run either **finds** a trivialization (and
  replays the full move path from the initial presentation, verifying it
  lands on the trivial pair — printed in the log), is **EXHAUSTED**
  (rigorous: no move sequence trivializes the input while every intermediate
  relator stays within the cap; note the claim is for single-letter
  conjugations — a classical conjugation by a long word is a chain of these,
  and its intermediate stages must also respect the cap), or is
  **TRUNCATED** at the state budget (inconclusive).

Build and reproduce:

```
gcc -O2 -o ac_search ac_search.c
./ac_search "xxYYY"   "xyxYXY" --cap 11 --max-states 60000000    # AK(2)
./ac_search "xxxYYYY" "xyxYXY" --cap 12 --max-states 100000000   # AK(3), exhausts
./ac_search "xxxYYYY" "xyxYXY" --cap 20 --max-states 150000000   # AK(3), deep (~11 GB RAM)
```

## Results on the Akbulut–Kirby family AK(n) = ⟨x, y | xⁿ = yⁿ⁺¹, xyx = yxy⟩

| Run | Outcome | States | Log |
|---|---|---|---|
| AK(2), cap 11 | **Trivialized: 21 moves, replay-verified, 0.5 s** (peak intermediate relator length 7) | 24,103 stored | `runs_ak2_cap11.log` |
| AK(3), cap 12 | **EXHAUSTED — rigorously no trivialization with all relators ≤ 12 letters**; best total length reached: 13 (= start) | 99,344,336 stored = expanded | `runs_ak3_cap12_exhaustive.log` |
| AK(3), cap 13 | Truncated at budget, no trivialization | 100M stored / 32.4M expanded | (numbers from run; deterministic) |
| AK(3), cap 20 | Truncated at budget, no trivialization, best total length still 13 | 150M stored / 18.2M expanded | `runs_ak3_cap20_greedy.log` |

Context: AK(2) is known AC-trivializable (Miasnikov, *IJAC* 1999). AK(3) is
the unique minimal potential counterexample — every balanced 2-generator
presentation of total length ≤ 13 is AC-trivializable except possibly AK(3)
(Havas–Ramsay, *IJAC* 2003; certificate frontier since pushed to length 14,
arXiv:2607.23611). The RL attack of Shehper et al. (arXiv:2408.15332,
NeurIPS 2025) resolved many Miller–Schupp candidates but no AK(n), n ≥ 3.
AK(3) is now known *stably* AC-trivial (same line of work; also Lisitsa,
arXiv:2501.18601), so it can only witness failure of the standard (unstable)
conjecture.

**The disproof gap, precisely:** a counterexample claim is the Π₁ statement
"no finite AC sequence trivializes this presentation." Search can only remove
candidates or fail; certification needs an invariant constant on AC-classes
and nontrivial on some balanced presentation of the trivial group. None is
known, and Borovik–Lubotzky–Myasnikov's finite-quotient no-go rules out the
most natural source of one. That — not search-space size — is why nobody has
disproved Andrews–Curtis despite most experts believing it false.
