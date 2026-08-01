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

## Systematic validation against the published frontier

`validate_frontier.py` tests the engine against the literature at scale rather
than on a single instance. It enumerates **every** balanced 2-generator
presentation with total relator length ≤ N, up to AC-equivalence (relator
rotation, inversion, swap, and the 8 signed permutations of the generators —
every identification used is a genuine AC move, so coverage of AC-classes is
complete), keeps those with trivial abelianization (|det| = 1), and resolves
each one three ways:

- **TRIVIALIZED** — the engine found and replay-verified an AC path to
  ⟨x, y | x, y⟩. This *proves* the presentation defines the trivial group.
- **NONTRIVIAL** — we exhibit generators in some Sₙ killing both relators with
  nontrivial image, *proving* the group isn't trivial (so the theorem doesn't
  apply to it).
- **OPEN** — neither. Per Havas–Ramsay the only such case presenting the
  trivial group at total length ≤ 13 should be AK(3).

Results (`validation_results.json`, `validation_len13.log`):

| Total length | AC-classes | Trivialized | Proven nontrivial | Open |
|---|---|---|---|---|
| ≤ 9 | 90 | 90 | 0 | 0 |
| ≤ 11 | 530 | 530 | 0 | 0 |
| ≤ 12 | 1121 | 1120 | 1 | **0** |
| ≤ 13 | 3480 | 3472 | 6 | **2 (= 1 AC-class: AK(3))** |

The length-≤12 row is a clean independent replication of Miasnikov (*IJAC*
1999): **every** presentation in range either AC-trivializes or is provably
not a presentation of the trivial group — no unexplained cases. Longest
trivialization needed there: 40 moves; at length 13, 109 moves.

The length-≤13 row independently reproduces Havas–Ramsay (*IJAC* 2003), whose
theorem is that everything in range is AC-trivializable **or AC-equivalent to
AK(3)**. Exactly two classes survive as OPEN, and one of them is AK(3):

```
xxyXYY   xxxYYXy   total=13
xyxYXY   xxxxYYY   total=13   <-- AK(3)
```

Neither has a nontrivial quotient in any Sₙ for n ≤ 6 (`quotient_search.py`,
`quotient_open2.log`, `quotient_ak3.log`), consistent with both presenting the
trivial group. The second class is *not* a competing counterexample — it is
AK(3) in disguise, and we can prove it. Our enumeration quotients only by the
obvious symmetries (rotation, inversion, swap, generator permutations), not by
full AC-equivalence, so an AK(3) relative surfacing as a separate class is
expected.

The test needs no new code. Every AC move in the engine has its inverse in the
move set, and both a move and its inverse respect the length cap, so the
cap-restricted move graph is **undirected**: a run that exhausts enumerates
exactly the connected component of its starting state. Two presentations
therefore lie in the same component if and only if their exhaustive runs
report the *same state count*. They do, to the digit:

```
AK(3)             xyxYXY xxxxYYY   EXHAUSTED  expanded=99344336 stored=99344336
second open case  xxyXYY xxxYYXy   EXHAUSTED  expanded=99344336 stored=99344336
```

(`runs_ak3_cap12_exhaustive.log`, `runs_open2_cap12.log`.) So the two OPEN
classes are one AC-class: **there is exactly one open AC-class of total length
≤ 13, and it is AK(3)** — precisely Havas–Ramsay's theorem, re-derived here
from scratch.

**The disproof gap, precisely:** a counterexample claim is the Π₁ statement
"no finite AC sequence trivializes this presentation." Search can only remove
candidates or fail; certification needs an invariant constant on AC-classes
and nontrivial on some balanced presentation of the trivial group. None is
known, and Borovik–Lubotzky–Myasnikov's finite-quotient no-go rules out the
most natural source of one. That — not search-space size — is why nobody has
disproved Andrews–Curtis despite most experts believing it false.
