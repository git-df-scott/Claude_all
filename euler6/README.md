# Hunting a sixth-power counterexample to Euler's sum of powers conjecture

Target equation ((6,1,5) in Lander–Parkin–Selfridge notation):

```
a^6 + b^6 + c^6 + d^6 + e^6 = f^6,   positive integers.
```

Euler's conjecture (1769) says no such solution exists. The k=5 case fell in
1966 (Lander–Parkin, 27^5+84^5+110^5+133^5 = 144^5) and k=4 in 1988
(Elkies, plus Frye's minimal 95800^4+217519^4+414560^4 = 422481^4). For
k ≥ 6 the question is open: no solution is known, and prior exhaustive
searches (Meyrignac's EulerNet ecosystem, ca. 2002) established there is
**no solution with f ≤ 730,000**.

## Structure theorem (elementary, but it carries the whole project)

Sixth powers are extremely constrained modulo 7, 9 and 8:

- mod 7: x^6 ≡ 0 (7|x) or 1 (Fermat)
- mod 9: x^6 ≡ 0 (3|x) or 1 (Euler, φ(9)=6)
- mod 8: x^6 ≡ 0 (x even) or 1

So mod each of 7, 9, 8, the left side counts its terms coprime to that
prime, and the count (0..5) must equal f's own indicator (0 or 1). A count
of 0 forces every term and f to share the prime, contradicting primitivity.
Hence in a **primitive** solution:

> exactly one of a..e is odd, exactly one is coprime to 3, exactly one is
> coprime to 7 — and f is coprime to 42.

(The three "exempt" slots may coincide on one term or spread over up to
three terms.) Since any solution is a scalar multiple of a primitive one
with smaller f, searching primitives is fully general.

## Programs

### `search.c` — exhaustive DFS searcher

Depth-first over a ≥ b ≥ c ≥ d ≥ e with:

- exact analytic bounds per level (largest remaining term satisfies
  x^6 ≥ R/j, x^6 ≤ R);
- the class constraints above tracked as a bitmask; once a class bit is
  spent, all deeper terms are divisible by that prime, so the loops stride
  by lcm(spent primes) up to 42;
- residue feasibility masks: achievable residues of sums of j sixth powers
  mod 13 and 43 (sixth powers hit 3 of 13 and 8 of 43 residues), with
  remainders tracked incrementally (no 128-bit division in the hot path);
- O(1) perfect-sixth-power finish.

**Validation:** run in fifth-power mode (`./search 5 2 150`) the identical
machinery re-discovers Lander–Parkin's counterexample in ~30 ms.
Sixth-power mode confirms the literature (no solutions) in every range it
has been run on. Cost grows like ~f^3.3, so this program alone cannot pass
the 730k frontier on small hardware — that motivates the next program.

### `caseA2.c` — deep sweep of the concentrated-exemption case

Suppose all three exemptions land on a single term t (necessarily coprime
to 42). Then the other four terms are divisible by 42, hence

```
f^6 ≡ t^6 (mod 42^6),   42^6 = 5,489,031,744.
```

Because f is a unit mod 42^6, t/f must be one of the **144 sixth roots of
unity mod 42^6** (4 mod 2^6 × 6 mod 3^6 × 6 mod 7^6, via CRT). For each f
there are at most 144 admissible t, of which only those with t < f matter:
the expected number of candidate pairs (f,t) with f ≤ F is

```
~ (2/7) · 144 · F² / (2 · 42^6)    — negligible until F is in the millions.
```

This sub-case can therefore be swept to heights that no generic
enumeration can touch. Each candidate needs `m = (f^6 - t^6)/42^6` to be a
sum of four sixth powers with bases ≤ (f-1)/42, tested by:

1. **Exact class budgets.** In the reduced equation the number of odd bases
   is exactly m mod 8, the number coprime to 3 exactly m mod 9, coprime to
   7 exactly m mod 7 (each must be ≤ 4). Budgets of 0 force divisibility
   (strides up to 42); budgets equal to the remaining count force
   coprimality.
2. **Residue masks** mod 64, 27, 49, 13, 43 on every DFS node
   (incrementally maintained).
3. **Blocked Bloom filter** over *all* pair sums x^6+y^6 (~16 bits/pair,
   8 hash bits confined to one 64-byte cache line, built once in parallel):
   the two largest parts are enumerated and the remainder is answered in
   one memory access. Positives get exact two-pointer verification, so
   false positives cost time, never correctness.
4. 128-bit arithmetic throughout; `m` is extracted from the factorization
   f^6−t^6 = (f−t)(f+t)(f²−ft+t²)(f²+ft+t²) with gcd-cascade division by
   42^6 to stay within 128 bits (f ≤ 10^8).

Sub-case coverage note: the complementary cases spread the exemptions over
two or three terms, leaving only three or two terms divisible by 42; their
candidate density scales like F³/42^6 or worse, so the same modular
leverage does not thin them below the generic cost — they remain covered
only up to the classical exhaustive bound. The concentrated case is
precisely the one where the 42-structure concentrates enough to break past
the frontier. Under a naive equidistribution model it carries ~1/25 of the
expected (logarithmically sparse) solution mass.

**Validation:** the 144 roots are verified by direct powering; DFS
self-tests recover planted decompositions (mixed-class and all-42 cases)
and reject a residue-infeasible target; the [700k, 730k] band — inside the
published exhaustively-searched region — runs clean (124 candidates, no
solutions), consistent with the literature.

## Results

- `search.c`, sixth powers, f ≤ 5,000: no solutions (consistency check).
- `caseA2.c`, concentrated case (running tally; see HANDOFF.md for the
  authoritative table):

| f range | candidates | solutions |
|---|---|---|
| 700,000–730,000 | 124 | 0 (inside known bound — consistency) |
| 730,000–1,000,000 | 1,314 | 0 |
| 1,000,000–1,500,000 | 3,523 | 0 |
| 1,500,000–2,000,000 | 5,129 | 0 |
| 2,000,000–2,500,000 | 7,222 | 0 |
| 2,500,000–3,200,000 | 12,443 | 0 |
| 3,200,000–4,000,000 | in progress | — |

## Where a real breakthrough would have to come from

The variety X: x⁶+y⁶+z⁶+u⁶+v⁶ = w⁶ in P⁵ is a smooth sextic fourfold with
trivial canonical bundle — a **Calabi–Yau fourfold**. That is structurally
encouraging: it is the same reason Elkies' quartic surface (a K3, the CY
surface case) was not obstructed by Bombieri–Lang and ultimately carried
rational points found via its elliptic fibrations. No analogous fibration
argument over ℚ is currently known for the sextic fourfold; finding an
elliptic fibration on X defined over ℚ with a section of infinite order —
or a rational curve on X over ℚ with all coordinates nonzero — would be
the genuine Elkies-style path. The Shioda inductive structure (Fermat
varieties dominated by products of Fermat curves) exists only over
cyclotomic extensions, and the Fermat sextic *curve* has no nontrivial
rational points (FLT, n=6), so no points descend that way.

## Build & run

```
gcc -O3 -march=native -fopenmp -o search search.c -lm
./search 5 2 150          # validation: finds Lander-Parkin in ~30ms
./search 6 2 20000        # exhaustive sixth-power search

gcc -O3 -march=native -fopenmp -o caseA2 caseA2.c -lm
./caseA2 730000 3200000   # concentrated-case sweep past the frontier
```
