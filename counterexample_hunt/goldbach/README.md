# Goldbach module

Small computational module for the counterexample-hunt report. Everything below
is from an actual run (`python3 goldbach.py`, seed 20260801, Python 3.11.15 +
numpy 2.4.6, total runtime 15.4 s; console transcript in `run.log`, full
machine-readable data in `results.json`).

## What was run

1. **Exhaustive sweep, even n in [4, 10^8]** (numpy sieve, 49,999,999 evens,
   12.2 s). Every even n decomposes as p + q with both prime. The tracked
   statistic is the maximum over n of the *smallest* prime p usable in a
   decomposition:

   - **max smallest-p = 1093, first attained at n = 60,119,912** — i.e. every
     even n ≤ 10^8 is a sum of two primes with the smaller prime ≤ 1093.
   - The full record sequence (n, smallest-p records: 6→3, 12→5, 30→7, 98→19,
     …, 37,998,938→1039, 60,119,912→1093) is in `results.json` and matches
     OEIS A025018/A025019, which independently validates the sweep.

2. **50 random even n in [4×10^18, 10^19]** — just above the verified
   frontier. For each: smallest-p decomposition and r(n) = number of
   decompositions with p < 10^4. Primality here is *proven* (deterministic
   Miller–Rabin, valid below 3.317×10^24).

   - smallest p: min 3, median 37, max 647
   - r(n), p < 10^4: min 28, median 45, max 118

3. **20 random 100-digit even n.** Same statistics; primality is
   **probable prime only** (strong Miller–Rabin, base 2 + 40 random bases —
   not a proof, no Lucas/BPSW step).

   - smallest p: min 23, median 435, max 4651
   - r(n), p < 10^4: min 3, median 9.5, max 23

## Status of the conjecture

Goldbach has been verified exhaustively for all even n ≤ 4×10^18
(T. Oliveira e Silva, S. Herzog, S. Pardi, *Empirical verification of the even
Goldbach conjecture and computation of prime gaps up to 4·10^18*, Math. Comp.
83 (2014), 2033–2060). Nothing this module can reach is new territory in the
exhaustive range, and the samples above the frontier are spot checks, not a
verification.

## Why the search EV is ~0

The Hardy–Littlewood heuristic predicts the number of representations of an
even n as a sum of two primes grows like n·C₂·∏((p−1)/(p−2))/(ln n)², i.e.
without bound — and the data above shows exactly that: even restricting the
smaller prime to p < 10^4, every sampled 19-digit n already has dozens of
representations, and 100-digit n still have several. A counterexample would
therefore require n − p to be composite *simultaneously* for every prime
p < n/2 — a conspiracy across all primes below n whose heuristic probability
decays faster than exponentially in n/(ln n)². That is why the expected value
of hunting for a Goldbach counterexample is essentially zero. The conjecture
stays on counterexample-hunt lists only because its counterexample would at
least be a finite, instantly checkable object: a single even number that
fails, verifiable by anyone in seconds.
