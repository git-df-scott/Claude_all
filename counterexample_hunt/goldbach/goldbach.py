#!/usr/bin/env python3
"""Goldbach conjecture computational module.

Part of an honest counterexample hunt: Goldbach is verified exhaustively to
4e18 (Oliveira e Silva, Herzog, Pardi, 2014), so nothing here can plausibly
find a counterexample. This module (a) re-verifies a small window [4, 1e8]
exhaustively and tracks the maximum over n of the *smallest* prime p in any
decomposition n = p + q, and (b) samples even n just above the verified
frontier and at 100 digits, finding smallest-p decompositions and counting
representations with p < 1e4, to show representation counts grow rather than
thin out (Hardy-Littlewood: ~ n / (ln n)^2 scale).

Primality: deterministic Miller-Rabin below 3.317e24 (proven, known base set);
above that, strong-probable-prime Miller-Rabin with base 2 plus 40 random
bases (reported as probable prime, not proven).
"""

import argparse
import json
import math
import random
import statistics
import sys
import time

import numpy as np

# Deterministic Miller-Rabin: these bases prove primality for all
# n < 3,317,044,064,679,887,385,961,981 (Sorenson & Webster 2015).
DET_BASES = (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37)
DET_LIMIT = 3317044064679887385961981
N_RANDOM_BASES = 40

_RNG_BASES = random.Random(0xB0B5)  # only for random MR bases above DET_LIMIT


def _miller_rabin(n, bases):
    """Strong-probable-prime test of odd n > 2 against each base."""
    d = n - 1
    s = (d & -d).bit_length() - 1
    d >>= s
    for a in bases:
        a %= n
        if a < 2:
            continue
        x = pow(a, d, n)
        if x == 1 or x == n - 1:
            continue
        for _ in range(s - 1):
            x = x * x % n
            if x == n - 1:
                break
        else:
            return False
    return True


# Product of odd primes in (37, 5000): one gcd rejects most composites cheaply.
def _primorial(lo, hi):
    r = 1
    for m in range(lo, hi):
        if all(m % p for p in range(2, int(m**0.5) + 1)):
            r *= m
    return r


_TRIAL_PRIMORIAL = _primorial(38, 5000)


def is_prime(n):
    """True if n is prime (proven below DET_LIMIT, probable above)."""
    if n < 2:
        return False
    for p in DET_BASES:  # doubles as a tiny trial-division list
        if n % p == 0:
            return n == p
    if n < DET_LIMIT:
        return _miller_rabin(n, DET_BASES)
    if math.gcd(n, _TRIAL_PRIMORIAL) > 1:
        return False
    bases = [2] + [_RNG_BASES.randrange(3, n - 1) for _ in range(N_RANDOM_BASES)]
    return _miller_rabin(n, bases)


def build_sieve(limit):
    """Boolean primality array for 0..limit."""
    s = np.ones(limit + 1, dtype=bool)
    s[:2] = False
    for i in range(2, math.isqrt(limit) + 1):
        if s[i]:
            s[i * i :: i] = False
    return s


def exhaustive_sweep(limit, chunk_evens=2**21):
    """Verify Goldbach for every even n in [4, limit].

    For each even n, finds the smallest prime p with n - p prime (which
    automatically satisfies p <= n/2: if the first hit had p > n - p, the
    prime n - p would have hit earlier). Returns max smallest-p, the first n
    attaining it, and the full record sequence.
    """
    t0 = time.time()
    sieve = build_sieve(limit)
    table_top = min(limit, 50000)
    odd_primes = [int(p) for p in np.flatnonzero(sieve[: table_top + 1])[1:]]

    max_p = 2  # n = 4 = 2 + 2 is the only case using p = 2
    records = [[4, 2]]  # running records of smallest-p, in increasing n
    evens_verified = 1

    for start in range(6, limit + 1, chunk_evens * 2):
        stop = min(start + chunk_evens * 2, limit + 1)
        nvals = np.arange(start, stop, 2, dtype=np.int64)
        evens_verified += nvals.size
        small = np.zeros(nvals.size, dtype=np.int32)  # smallest p per n
        idx = np.arange(nvals.size)  # still-unresolved positions (sorted)
        for p in odd_primes:
            if idx.size == 0:
                break
            if int(nvals[idx[0]]) < 2 * p:
                # All primes p' < p failed for this n, and any q > n/2 hit
                # would have been found from the small side first.
                raise AssertionError(
                    f"Goldbach counterexample candidate: n = {int(nvals[idx[0]])}"
                )
            hit = sieve[nvals[idx] - p]
            small[idx[hit]] = p
            idx = idx[~hit]
        if idx.size:
            raise RuntimeError(
                f"prime table (<= {table_top}) exhausted at n = {int(nvals[idx[0]])}"
            )
        # True running records: smallest_p(n) > smallest_p(m) for all m < n,
        # including all m in previous chunks (threshold seeded with max_p).
        runmax = np.maximum.accumulate(small)
        threshold = np.empty_like(runmax)
        threshold[0] = max_p
        np.maximum(runmax[:-1], max_p, out=threshold[1:])
        for i in np.flatnonzero(small > threshold):
            records.append([int(nvals[i]), int(small[i])])
        max_p = records[-1][1]
        print(
            f"  sweep: verified up to {stop - 1:,} "
            f"(current max smallest-p = {max_p} at n = {records[-1][0]:,})",
            file=sys.stderr,
        )
    argmax_n = records[-1][0]
    assert all(a[0] < b[0] and a[1] < b[1] for a, b in zip(records, records[1:])), \
        "record sequence must strictly increase in n and p"

    return {
        "range": [4, limit],
        "evens_verified": evens_verified,
        "max_smallest_p": max_p,
        "attained_at_n": argmax_n,
        "smallest_p_records": records,
        "elapsed_seconds": round(time.time() - t0, 2),
    }


def sample_group(name, numbers, odd_primes, restricted_bound=10**4):
    """For each even n: smallest-p decomposition and count of p < bound."""
    t0 = time.time()
    restricted_primes = [p for p in odd_primes if p < restricted_bound]
    samples = []
    for n in numbers:
        smallest_p = None
        count = 0
        for p in restricted_primes:
            if is_prime(n - p):
                count += 1
                if smallest_p is None:
                    smallest_p = p
        if smallest_p is None:  # would be a counterexample candidate; dig deeper
            for p in odd_primes:
                if is_prime(n - p):
                    smallest_p = p
                    break
            else:
                raise AssertionError(f"no decomposition found for n = {n}")
        q = n - smallest_p
        assert is_prime(smallest_p) and smallest_p + q == n
        samples.append(
            {"n": str(n), "smallest_p": smallest_p, "q": str(q),
             "restricted_count": count}
        )
        print(f"  {name}: n has {len(str(n))} digits, smallest p = {smallest_p}, "
              f"r(n; p<1e4) = {count}", file=sys.stderr)

    ps = [s["smallest_p"] for s in samples]
    rs = [s["restricted_count"] for s in samples]
    return {
        "count": len(samples),
        "restricted_bound": restricted_bound,
        "smallest_p_stats": {"min": min(ps), "median": statistics.median(ps),
                             "max": max(ps)},
        "restricted_count_stats": {"min": min(rs), "median": statistics.median(rs),
                                   "max": max(rs)},
        "samples": samples,
        "elapsed_seconds": round(time.time() - t0, 2),
    }


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--limit", type=int, default=10**8,
                    help="exhaustive sweep upper bound (default 1e8)")
    ap.add_argument("--frontier-samples", type=int, default=50)
    ap.add_argument("--huge-samples", type=int, default=20)
    ap.add_argument("--seed", type=int, default=20260801)
    ap.add_argument("--out", default="results.json")
    args = ap.parse_args()

    t0 = time.time()
    rng = random.Random(args.seed)

    print(f"[1/3] Exhaustive sweep of even n in [4, {args.limit:,}] ...",
          file=sys.stderr)
    exhaustive = exhaustive_sweep(args.limit)

    odd_primes = [int(p) for p in np.flatnonzero(build_sieve(10**6))[1:]]

    # Even n uniform in [4e18, 1e19]: just above the verified frontier 4e18.
    frontier_ns = [2 * rng.randrange(2 * 10**18, 5 * 10**18)
                   for _ in range(args.frontier_samples)]
    print(f"[2/3] {len(frontier_ns)} samples in [4e18, 1e19] "
          "(primality proven, deterministic MR) ...", file=sys.stderr)
    frontier = sample_group("frontier", frontier_ns, odd_primes)
    frontier["range"] = ["4e18", "1e19"]
    frontier["primality"] = "proven (deterministic Miller-Rabin, n < 3.317e24)"

    huge_ns = []
    while len(huge_ns) < args.huge_samples:
        n = rng.randrange(10**99, 10**100) & ~1
        huge_ns.append(n)
    print(f"[3/3] {len(huge_ns)} random 100-digit even samples "
          "(probable-prime only) ...", file=sys.stderr)
    huge = sample_group("100-digit", huge_ns, odd_primes)
    huge["range"] = ["1e99", "1e100"]
    huge["primality"] = ("probable prime (strong Miller-Rabin, base 2 + "
                         f"{N_RANDOM_BASES} random bases; not a proof)")

    results = {
        "meta": {
            "generated_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
            "seed": args.seed,
            "python": sys.version.split()[0],
            "numpy": np.__version__,
            "verified_frontier_note": (
                "Goldbach verified exhaustively to 4e18 by Oliveira e Silva, "
                "Herzog & Pardi (Math. Comp. 83, 2014)."),
            "total_runtime_seconds": None,  # filled below
        },
        "exhaustive": exhaustive,
        "frontier_samples": frontier,
        "huge_samples": huge,
    }
    results["meta"]["total_runtime_seconds"] = round(time.time() - t0, 2)

    with open(args.out, "w") as f:
        json.dump(results, f, indent=1)

    fp, fr = frontier["smallest_p_stats"], frontier["restricted_count_stats"]
    hp, hr = huge["smallest_p_stats"], huge["restricted_count_stats"]
    print(f"""
== Goldbach module summary ==
Exhaustive sweep: every even n in [4, {args.limit:,}] has a prime
  decomposition ({exhaustive['evens_verified']:,} evens verified,
  {exhaustive['elapsed_seconds']} s).
  Max over n of the SMALLEST prime in a decomposition:
    p = {exhaustive['max_smallest_p']} first attained at n = {exhaustive['attained_at_n']:,}
  (i.e. every even n <= {args.limit:.0e} is prime + prime with the small
   prime <= {exhaustive['max_smallest_p']}).

Frontier samples ({frontier['count']} even n in [4e18, 1e19], primality PROVEN):
  smallest p:            min {fp['min']}, median {fp['median']}, max {fp['max']}
  r(n) with p < 1e4:     min {fr['min']}, median {fr['median']}, max {fr['max']}

100-digit samples ({huge['count']} even n, PROBABLE primes only):
  smallest p:            min {hp['min']}, median {hp['median']}, max {hp['max']}
  r(n) with p < 1e4:     min {hr['min']}, median {hr['median']}, max {hr['max']}

All sampled n decompose easily; restricted representation counts stay far
from 0, consistent with Hardy-Littlewood growth ~ n/(ln n)^2. Full data:
{args.out}
Total runtime: {results['meta']['total_runtime_seconds']} s""")


if __name__ == "__main__":
    main()
