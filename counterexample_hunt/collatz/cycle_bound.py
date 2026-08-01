"""Rigorous lower bound on the size of any nontrivial Collatz cycle.

Derivation (all classical; self-contained here):

Work with the odd-to-odd ("Syracuse") map: for odd n, the next odd value is
(3n+1)/2^a. Suppose a nontrivial cycle has L odd members n_1,...,n_L and
K = sum of the a_i (total halvings). Multiplying n_{i+1} = 3n_i(1+1/(3n_i))/2^{a_i}
around the cycle gives

    2^K / 3^L = prod_i (1 + 1/(3 n_i)),   hence   0 < K - L*log2(3) <= L*eps,

where eps = log2(1 + 1/(3*n_min)) and n_min is the smallest cycle member.
Exhaustive verification (D. Barina, J. Supercomputing 2025: every n < 2^71
converges) forces n_min > M = 2^71, so eps <= log2(1 + 1/(3M)) ~ 2e-22:
K/L approximates alpha = log2(3) from above extraordinarily well.

The law of best approximation (Khinchin, Continued Fractions, Thm 17): if
p_k/q_k are the continued-fraction convergents of an irrational alpha, then
for all integers K, L with 1 <= L < q_{k+1}:  |L*alpha - K| >= |q_k*alpha - p_k|.

So whenever delta_k := |q_k*alpha - p_k| >= q_{k+1}*eps, a cycle with
L < q_{k+1} would give delta_k <= |L*alpha - K| = K - L*alpha <= L*eps
< q_{k+1}*eps <= delta_k, a contradiction; hence L >= q_{k+1}.

Everything below is computed in fixed-point INTEGER interval arithmetic
(scale 10^80); no floating point enters any inequality. Convergents with
small denominators are additionally checked exactly via big-integer
comparison of 2^p vs 3^q.
"""

import math

S = 10 ** 80  # fixed-point scale

M_2_71 = 2 ** 71  # Barina 2025 verification frontier
M_2_68 = 2 ** 68  # Barina 2020 frontier, for a conservative variant


def atanh_inv_interval(k):
    """Interval [lo, hi] (scale S) containing atanh(1/k), k >= 2."""
    total = 0
    j = 0
    while True:
        term = S // ((2 * j + 1) * k ** (2 * j + 1))
        if term == 0:
            break
        total += term
        j += 1
    # each floored term undercounts by < 1 ulp; tail of the series is bounded
    # by its first omitted term times k^2/(k^2-1)
    tail = (S * k * k) // ((2 * j + 1) * k ** (2 * j + 1) * (k * k - 1)) + 1
    return total, total + j + tail


def interval_div(a_lo, a_hi, b_lo, b_hi):
    """[a]/[b] in scale S, for positive intervals."""
    return (a_lo * S) // b_hi, (a_hi * S) // b_lo + 1


def cf_expansion(lo, hi, max_terms=60):
    """Continued-fraction terms of the number in [lo, hi]/S, while the
    interval determines them unambiguously."""
    terms = []
    while len(terms) < max_terms:
        a_lo, a_hi = lo // S, hi // S
        if a_lo != a_hi:
            break
        terms.append(a_lo)
        f_lo, f_hi = lo - a_lo * S, hi - a_lo * S
        if f_lo <= 0:
            break
        lo, hi = (S * S) // f_hi, (S * S) // f_lo + 1
    return terms


def convergents(terms):
    ps, qs = [], []
    p2, p1, q2, q1 = 0, 1, 1, 0
    for a in terms:
        p2, p1 = p1, a * p1 + p2
        q2, q1 = q1, a * q1 + q2
        ps.append(p1)
        qs.append(q1)
    return ps, qs


def main():
    # ln2 = 2*atanh(1/3),  ln3 = ln2 + 2*atanh(1/5)
    t3_lo, t3_hi = atanh_inv_interval(3)
    t5_lo, t5_hi = atanh_inv_interval(5)
    ln2 = (2 * t3_lo, 2 * t3_hi)
    ln3 = (ln2[0] + 2 * t5_lo, ln2[1] + 2 * t5_hi)
    a_lo, a_hi = interval_div(ln3[0], ln3[1], ln2[0], ln2[1])
    width = a_hi - a_lo
    print(f"alpha = log2(3) in [{a_lo}, {a_hi}] / 10^80  (interval width {width} ulp)")
    assert abs(a_lo / S - math.log2(3)) < 1e-14, "float cross-check failed"

    terms = cf_expansion(a_lo, a_hi)
    print(f"continued fraction of log2(3): {terms[:25]} ...  ({len(terms)} terms)")
    ps, qs = convergents(terms)

    # exact big-integer side check for small convergents:
    # p/q > log2(3)  <=>  2^p > 3^q  — and convergents must alternate sides
    prev_side = None
    for p, q in zip(ps, qs):
        if q > 10 ** 6:
            break
        side = pow(2, p) > pow(3, q)
        assert side != prev_side, "convergents failed to alternate"
        prev_side = side
    print("exact 2^p vs 3^q side check passed for all convergents with q <= 10^6")

    for label, M in (("2^71 (Barina 2025)", M_2_71), ("2^68 (Barina 2020)", M_2_68)):
        # eps = log2(1 + 1/(3M)) <= (1/(3M)) / ln2  (since ln(1+u) < u)
        eps_hi = ((S // (3 * M) + 1) * S) // ln2[0] + 1
        best_L = None
        for k in range(1, len(ps) - 1):
            d_lo = qs[k] * a_lo - ps[k] * S
            d_hi = qs[k] * a_hi - ps[k] * S
            assert d_lo > 0 or d_hi < 0, f"delta interval straddles 0 at k={k}"
            delta_lo = d_lo if d_lo > 0 else -d_hi
            if delta_lo >= qs[k + 1] * eps_hi:
                best_L = max(best_L or 0, qs[k + 1])
        assert best_L is not None
        # K > L*alpha, so standard-map cycle length L + K > L*(1 + alpha)
        k_min = (best_L * a_lo) // S + 1
        print(f"\nverified frontier M = {label}: any nontrivial cycle has")
        print(f"  odd members  L >= {best_L:,}")
        print(f"  halvings     K >= {k_min:,}")
        print(f"  standard-map cycle length L + K >= {best_L + k_min:,}")


if __name__ == "__main__":
    main()
