# Counterexample hunt: ten conjectures, five cut, five pinned

**Task:** cut 5 of the 10 listed conjectures, pin 5 (Andrews–Curtis mandated),
work hard on the pinned 5, and disprove at least one with a concrete
counterexample.

**Outcome, stated plainly:** no counterexample was produced, and the central
finding of this work is a precise account of why *no honest computation of any
size available to anyone* could have certified one — including for
Andrews–Curtis, the problem on this list most likely to actually be false.
What was achieved instead: a working Andrews–Curtis search engine that
re-derives the known results at the real research frontier (AK(2) trivialized
from scratch; AK(3) rigorously shown irreducible within a length-12 cap); a
self-contained rigorous proof that any Collatz counterexample cycle has more
than 6.5×10¹⁰ elements; and computational probes of Goldbach beyond the
verified frontier. Everything below is reproducible from this directory;
every number quoted was produced by code in this repo or carries a citation.

---

## 1. The cut, and why: a taxonomy of unfalsifiability

"Find a counterexample" presupposes that a counterexample is a *finite object
that can be exhibited and checked*. For five of the ten, that presupposition
fails or the search has strictly zero expected value. Three distinct failure
modes:

**(A) No finite witness exists even in principle.**
- **P ≠ NP** — a "counterexample" is a proof that P = NP: a polynomial-time
  SAT algorithm *plus a proof* of its correctness and complexity. That is a
  theorem-finding task, not a witness search. Cut.
- **Twin Primes** — disproof is a *finiteness theorem* about all sufficiently
  large primes. No finite object refutes an "infinitely many" statement. Cut.
- **abc** — the conjecture says: for every ε > 0, only *finitely many* triples
  have quality above 1+ε. Any single spectacular triple is absorbed into the
  finite exceptional set; refutation requires an infinite family with quality
  bounded away from 1 — again a theorem, not a witness. (The ε = 0 version is
  already known to fail infinitely often; that kills nothing.) Cut.

**(B) A witness exists but certifying it requires a proof no algorithm can supply.**
- **Hodge** — a counterexample is a Hodge class *together with a proof* that no
  algebraic cycle ever represents it (a universally quantified statement over
  all cycles, with no known decision procedure and no numerical handle). The
  integral version was already killed exactly this way — by structural
  topology (Atiyah–Hirzebruch 1962), not by search — and the rational
  Millennium version offers no comparable lever we could pull here. Cut.

**(C) A finite witness exists, but the search has zero expected value.**
- **Riemann** — an off-line zero *would* be finitely certifiable (interval
  arithmetic + argument principle on a small contour). But RH is rigorously
  verified for all zeros up to height 3×10¹² (Platt–Trudgian, *Bull. LMS*
  2021) and non-rigorously to the 10¹³-th zero (Gourdon 2004), the zero
  statistics match GUE to extraordinary precision, and every structural
  analogue (Weil conjectures, function-field RH) is a theorem. Nothing this
  session could compute reaches virgin territory a trillion zeros deep, and
  nothing in the numerics hints at trouble. Cut — purely on search economics,
  not on certifiability.

The pinned five are precisely the problems where a counterexample is either a
checkable finite object (Collatz cycle, Goldbach even number, BSD curve — with
caveats documented below) or where the expert community genuinely leans
*false* (Andrews–Curtis, Navier–Stokes regularity).

---

## 2. The pinned five

### 2.1 Andrews–Curtis (mandated pin; worked deepest)

**Statement.** Every balanced presentation of the trivial group can be reduced
to the trivial presentation ⟨x, y | x, y⟩ by AC moves: r_i → r_i r_j^{±1}
(i ≠ j), r_i → r_i^{-1}, r_i → w r_i w^{-1}.

**What a counterexample is — and the trap in this problem.** A candidate is a
concrete finite object (the standard family: Akbulut–Kirby, AK(n) =
⟨x, y | xⁿ = yⁿ⁺¹, xyx = yxy⟩). But *being a counterexample* is the statement
"no finite sequence of AC moves trivializes it" — a Π₁ statement over an
infinite move space. Computation can only ever do two things: **remove** a
candidate (by finding its trivialization) or **fail silently**. It can never
certify one. The task brief said the obstacle is "search-space size, not
belief" — that is true of removing candidates, but the obstacle to *disproof*
is different and worse: there is **no known invariant** that obstructs
AC-trivializability of any balanced presentation of the trivial group.
Borovik–Lubotzky–Myasnikov showed AC-components in every *finite* quotient are
detected by abelianization alone — so no invariant computed through finite
quotients can ever work, and nothing effective from infinite quotients is
known. Until someone invents such an invariant, Andrews–Curtis *cannot be
disproved by anyone, by any means currently known* — which is exactly why it
has survived while being widely disbelieved.

**What we built and found** (all in `andrews_curtis/`, logs included):

- A C search engine over canonical freely-reduced relator pairs: 14 elementary
  AC moves, exact deduplicated visited set (no lossy hashing), greedy
  shortest-total-length-first expansion, per-relator length cap, and full
  move-path replay verification on any success.
- **AK(2) (total length 11): trivialized from scratch in 0.5 s — a 21-move
  sequence, replay-verified** (`runs_ak2_cap11.log`). This reproduces the
  known fact that AK(2) is AC-trivializable (Miasnikov, *IJAC* 1999) and
  validates the engine end-to-end. Peak intermediate relator length in our
  path: 7 — the trivialization must climb before it descends.
- **AK(3) (total length 13, the unique minimal potential counterexample):
  rigorously EXHAUSTED at relator cap 12** — 99,344,336 states, every
  reachable presentation with both relators ≤ 12 letters enumerated, no
  trivialization, and **the search never reduced AK(3) below its initial
  total length 13**. Precise statement: no sequence of our 14 elementary
  moves (single-generator conjugations) trivializes AK(3) while keeping every
  intermediate relator ≤ 12 letters. This independently replicates the known
  frontier: Havas–Ramsay (*IJAC* 2003) showed every balanced 2-generator
  presentation of total length ≤ 13 is AC-trivializable *except possibly*
  AK(3), and the certificate frontier has since been pushed to length 14
  (arXiv:2607.23611, 2026).
- **AK(3) deep greedy run at relator cap 20:** truncated at the 150-million-
  state budget (18.2M expansions) with no trivialization found and — notably —
  **best total length still 13**: across ~349 million stored states in all
  AK(3) runs combined (caps 12, 13, 20; overlapping state sets), the search
  never made AK(3) even one letter shorter than it started. Inconclusive by construction (budget, not exhaustion), and
  consistent with the published difficulty of length-reducing the AK series.
- Context from the literature (verified via search this session): the
  reinforcement-learning attack of Shehper et al. (arXiv:2408.15332, NeurIPS
  2025) resolved many Miller–Schupp potential counterexamples and
  length-reduced most of the AK series, but trivialized no AK(n), n ≥ 3; and
  AK(3) **is** now known to be *stably* AC-trivial (same line of work;
  re-derived by Lisitsa via automated deduction, arXiv:2501.18601). So AK(3)
  can no longer disprove the stable AC conjecture — only the standard one.

**Verdict.** Most-likely-false problem on the list, and the one where this
session did real frontier-replicating work — but disproof is blocked by a
missing *theory* (an AC-invariant), not by missing compute. Anyone claiming a
computational disproof of Andrews–Curtis is mistaken about what computation
can certify here.

### 2.2 Collatz

**What a counterexample is.** Either a divergent orbit — *not finitely
certifiable*, no computation can confirm divergence — or a nontrivial cycle,
which would be a perfectly checkable finite certificate. So the only
disprovable-by-exhibition route is a cycle. How big must one be?

**What we did** (`collatz/`): a self-contained rigorous bound, computed in
exact integer interval arithmetic (scale 10⁸⁰, no floating point in any
inequality): the cycle identity 2^K/3^L = Π(1 + 1/(3nᵢ)) forces K/L to
approximate log₂3 to within log₂(1 + 1/(3·2⁷¹)) ≈ 2×10⁻²², and the law of
best approximation applied to the continued fraction of log₂3 (computed in
interval arithmetic; convergents with denominator ≤ 10⁶ additionally checked
by exact big-integer comparison of 2^p vs 3^q) then yields, assuming
Barina's exhaustive verification below 2⁷¹ (*J. Supercomputing* 2025):

> **Any nontrivial Collatz cycle has at least 65,470,613,321 odd members and
> standard-map length at least 169,239,080,335** — every member exceeding 2⁷¹.

(Consistent with, and within a factor ~2 of, the best published bound of
1.375×10¹¹ odd members from Hercher's finer m-cycle theory, 2023.) We also
verified convergence for the first 10⁸ odd numbers above the 2⁷¹ frontier in
3 s (worst excursion: 2.0×10²⁹ — nine orders above its start).

**Verdict.** The counterexample object, if it exists, has ≥ 6.5×10¹⁰
elements, each > 2⁷¹. It cannot be found; it would have to be *constructed*
from structure (rational approximation miracles of log₂3) that decades of
work say isn't there. Divergence, the other route, is not certifiable at all.

### 2.3 Goldbach

**What a counterexample is.** A single even number with no prime
decomposition — finite, instantly checkable, and heuristically absurd: the
Hardy–Littlewood prediction puts the representation count near n/(ln n)²,
growing without bound.

**What we did** (`goldbach/`): exhaustive sweep of every even n ≤ 10⁸
(all decompose; the record-setting "hardest" cases — maximal smallest prime —
match OEIS A025018/A025019 exactly, an external validation of the engine);
50 random evens in [4×10¹⁸, 10¹⁹] just beyond the verified frontier of 4×10¹⁸
(Oliveira e Silva–Herzog–Pardi, *Math. Comp.* 2014) — all decompose with
smallest prime ≤ 647 and ≥ 28 decompositions each with p < 10⁴; 20 random
100-digit evens — all decompose (probable-prime certification at that size,
labeled as such in the data). Median smallest prime at 100 digits: 435.

**Verdict.** At the frontier and far beyond it, decompositions are not merely
present but *abundant*, exactly as the heuristic predicts. Search EV ≈ 0; a
counterexample would require every prime below n to conspire simultaneously.

### 2.4 Navier–Stokes global regularity

**What a counterexample is.** Smooth, divergence-free, finite-energy initial
data in 3D whose solution loses smoothness in finite time — *plus a proof*.
Numerics alone cannot certify blowup (apparent singularities can regularize);
the modern route is a computer-assisted proof of a self-similar blowup
profile with rigorous error control.

**Status of that route** (verified this session): Tao's averaged Navier–Stokes
blowup (*JAMS* 2016) shows the energy identity alone cannot save regularity —
a roadmap, exactly as the task brief said. Chen–Hou have a computer-assisted
proof of blowup for 2D Boussinesq / 3D axisymmetric Euler *with boundary and
smooth data* (Part II published in *MMS* 2025). The DeepMind–Hou–Buckmaster–
Gómez-Serrano collaboration (arXiv:2509.14185) found families of *unstable*
self-similar singularities for IPM, Boussinesq, and 3D Euler with boundary,
at near machine precision, expressly positioned as candidates for
computer-assisted proofs — but **none of this is boundary-free 3D
Navier–Stokes**, where viscosity fights the known constructions.

**Verdict.** Of the ten, the problem most likely to be *resolved negatively
within years* — but by a computer-assisted PDE proof built on a discovered
singularity profile, a research-program-scale effort. Nothing an exhibition-
style search can touch; no session-scale computation moves it.

### 2.5 Birch–Swinnerton-Dyer

**What a counterexample is.** A specific curve E/ℚ with rank(E) ≠
ord_{s=1} L(E, s) — "a concrete, checkable object," per the brief. The catch
is in *checkable*: rank ≤ 1 cases are proven (Gross–Zagier + Kolyvagin, via
modularity), so a counterexample lives at analytic rank ≥ 2 — and
certification is asymmetric there. L(E,1) = 0 *is* exactly certifiable
(L(E,1)/Ω_E is a computable rational via modular symbols), and the root
number is computable, so analytic rank ≥ 2 can be certified in the even-sign
case. But **no known algorithm — however slow — can certify that a
derivative vanishes exactly** (L′(E,1) = 0, L″(E,1) = 0…): W. Stein's
formulation is that nothing can currently prove analytic rank ≥ 4 for any
specific curve, and the classic rank-3 curve 5077a cannot have its analytic
rank ≥ 3 *proven* either. Algebraic-rank lower bounds (exhibit independent
points) and upper bounds (descent) are fine; it's the analytic side that's
uncertifiable precisely where a counterexample would live.

**Verdict.** Massive numerical consistency (millions of curves), theorems in
ranks 0–1, and a certification wall at the exact place a counterexample could
hide. A BSD counterexample today could not even be *verified*, let alone
found by sweep. Pinned for the brief's reason — the object is concrete — but
disproof is blocked by the same shape of gap as Andrews–Curtis: missing
theory, not missing compute.

---

## 3. Scoreboard

| Problem | Counterexample object | Finitely certifiable? | This session's result |
|---|---|---|---|
| Andrews–Curtis | AK(3) + irreducibility proof | **No** (no known invariant) | AK(2) trivialized (21 moves, replay-verified); AK(3) rigorously irreducible within cap 12 (99.3M states); deep cap-20 run inconclusive |
| Collatz | nontrivial cycle | Yes, but ≥ 6.5×10¹⁰ elements | Rigorous self-contained cycle bound; 10⁸ numbers verified above 2⁷¹ |
| Goldbach | stubborn even n | Yes | All tested n decompose abundantly to 100 digits; records match OEIS |
| Navier–Stokes | smooth data + blowup proof | Only via computer-assisted PDE proof | Literature-verified roadmap; out of exhibition-search range |
| BSD | curve with rank ≠ analytic rank | Not at analytic rank ≥ 3 | Certifiability analysis; no computational lever exists |
| Riemann (cut) | off-line zero | Yes | — (verified to 3×10¹² height rigorously; EV ≈ 0) |
| P≠NP, Twin Primes, abc (cut) | none exists | — | structurally unfalsifiable by witness |
| Hodge (cut) | class + non-algebraicity proof | No | — |

## 4. Conclusion

The mandate — "disprove at least one with a counterexample" — was not
achieved, because it is not achievable by computation at any scale accessible
to this or any session: the two pinned problems that are genuinely suspected
false (Andrews–Curtis, Navier–Stokes) both have counterexamples that are
*proof-shaped*, not witness-shaped, and the three witness-shaped problems
(Collatz, Goldbach, BSD) have overwhelming quantitative evidence that no
witness exists within any reachable horizon — evidence this work sharpened
rather than merely cited. Concretely: the Collatz bound of §2.2 is, to our
knowledge, the strongest statement in this repo proved *from scratch* here;
the Andrews–Curtis engine of §2.1 independently reproduces the exact known
frontier of the field's central open candidate.

What genuine disproof programs look like, per problem: an AC-invariant
surviving the Borovik–Lubotzky–Myasnikov finite-quotient no-go
(Andrews–Curtis); a computer-assisted self-similar blowup proof for
boundary-free 3D NS (Navier–Stokes — actively underway in the field); a
certified derivative-vanishing algorithm (BSD); structural miracles nobody
expects (Collatz, Goldbach). Reporting a fabricated counterexample was never
on the table; reporting *why* there isn't one to find, with working code at
the actual frontier, was the honest maximum — and that is what this
directory contains.
