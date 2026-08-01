"""Todd-Coxeter coset enumeration over the trivial subgroup H = 1.

For a balanced 2-generator presentation <x, y | r1, r2> the enumeration of
cosets of H = 1 is an enumeration of the group elements themselves, so:

  * completes with index 1   -> the group is TRIVIAL (a proof)
  * completes with index n>1 -> the group is finite of order n (so nontrivial)
  * exceeds the coset budget -> INCONCLUSIVE (may be infinite, may just be hard)

Only the third outcome is a non-answer; the first two are rigorous, because a
completed coset table is a concrete permutation representation of the group on
the cosets and the enumeration is a decision procedure whenever it halts.

Generators are the letters x, X, y, Y with X = x^-1 and Y = y^-1, matching the
rest of this directory.  Words are plain strings, e.g. "xxyXYY".

Algorithm: HLT (Haselgrove-Leech-Trotter) with scan-and-fill, plus the standard
union-find / queue based COINCIDENCE routine (Holt, *Handbook of Computational
Group Theory*, ch. 5).  The main loop is repeated until a full pass over all
live cosets produces neither a new coset nor a coincidence; that guarantees the
final table is closed under every generator and every relator traces back to its
own coset even if a late coincidence invalidated an early scan.  The result is
then independently re-checked by ``verify_table``.

Storage is flat ``array('i')`` blocks (4 bytes/entry), not dicts of objects, so
millions of cosets cost tens of megabytes.

Usage:
    python3 todd_coxeter.py R1 R2 [--max-cosets N]
    python3 todd_coxeter.py --selftest
"""

import sys
import time
from array import array

LETTERS = "xXyY"
GEN_INDEX = {"x": 0, "X": 1, "y": 2, "Y": 3}
INV = (1, 0, 3, 2)

DEFAULT_MAX_COSETS = 5_000_000

TRIVIAL = "TRIVIAL"
FINITE = "FINITE"
INCONCLUSIVE = "INCONCLUSIVE"


class _Overflow(Exception):
    """Raised internally when the coset budget is exhausted."""


def _zeros(count):
    a = array("i")
    if count > 0:
        a.frombytes(bytes(a.itemsize * count))
    return a


def parse_word(word):
    """'xxYY' -> [0, 0, 3, 3].  Raises ValueError on a bad letter."""
    out = []
    for ch in word:
        if ch not in GEN_INDEX:
            raise ValueError("bad generator letter %r (expected one of %s)"
                             % (ch, LETTERS))
        out.append(GEN_INDEX[ch])
    return out


def unparse_word(word):
    return "".join(LETTERS[g] for g in word)


def free_reduce(word):
    """Cancel adjacent g g^-1 pairs."""
    out = []
    for g in word:
        if out and out[-1] == INV[g]:
            out.pop()
        else:
            out.append(g)
    return out


def cyclic_reduce(word):
    """Free-reduce, then cancel a matching first/last pair.

    Replacing a relator by a cyclic conjugate of it does not change the normal
    closure, so this is safe and it shortens the scans.
    """
    w = free_reduce(word)
    while len(w) > 1 and w[0] == INV[w[-1]]:
        w = w[1:-1]
    return w


def _prepare(words, cyclic):
    """Normalise words into (forward, inverse-letter) index tuples.

    ``cyclic`` must be True only for relators.  Cyclic reduction replaces
    g w g^-1 by w, which preserves the normal closure (so it is safe for a
    relator) but *not* the subgroup generated (so it is wrong for a subgroup
    generator: <g w g^-1> != <w> in general).
    """
    out = []
    for r in words:
        w = cyclic_reduce(parse_word(r)) if cyclic else free_reduce(parse_word(r))
        if not w:
            continue  # freely trivial, carries no information
        out.append((tuple(w), tuple(INV[g] for g in w)))
    return out


def _enumerate(relators, max_cosets, subgroup_gens=()):
    """Core enumeration.  Returns a dict with the raw table and statistics."""
    rels = _prepare(relators, cyclic=True)
    subs = _prepare(subgroup_gens, cyclic=False)

    if max_cosets < 1:
        raise ValueError("max_cosets must be >= 1")
    if max_cosets > 2 ** 31 - 2:
        raise ValueError("max_cosets must fit a signed 32-bit array('i') entry")

    cap = min(max_cosets, 4096)
    table = _zeros(4 * (cap + 1))   # table[4*c + g], 0 == undefined
    p = _zeros(cap + 1)             # union-find: p[c] == c iff c is alive
    n = 0                           # cosets ever defined
    n_coinc = 0                     # coincidences processed

    def grow(need):
        nonlocal cap, table, p
        new = cap
        while new < need:
            new = min(max_cosets, new * 2)
        table.extend(_zeros(4 * (new - cap)))
        p.extend(_zeros(new - cap))
        cap = new

    def define(a, x):
        nonlocal n
        if n >= max_cosets:
            raise _Overflow
        n += 1
        c = n
        if c > cap:
            grow(c)
        p[c] = c
        table[4 * a + x] = c
        table[4 * c + INV[x]] = a
        return c

    def rep(k):
        """Union-find find, with path compression."""
        l = k
        while True:
            m = p[l]
            if m == l:
                break
            l = m
        m = k
        while True:
            nxt = p[m]
            if nxt == m:
                break
            p[m] = l
            m = nxt
        return l

    def coincidence(a, b):
        """Identify cosets a and b, propagating all induced identifications."""
        nonlocal n_coinc
        q = []

        def merge(k, l):
            k = rep(k)
            l = rep(l)
            if k != l:
                if k > l:
                    k, l = l, k
                p[l] = k
                q.append(l)

        merge(a, b)
        i = 0
        while i < len(q):
            e = q[i]
            i += 1
            n_coinc += 1
            for x in range(4):
                f = table[4 * e + x]
                if f == 0:
                    continue
                ix = INV[x]
                # Drop the dead edge in both directions before reinstalling it
                # on the surviving representatives.
                table[4 * e + x] = 0
                if table[4 * f + ix] == e:
                    table[4 * f + ix] = 0
                e1 = rep(e)
                f1 = rep(f)
                g = table[4 * e1 + x]
                if g != 0:
                    merge(f1, g)
                else:
                    h = table[4 * f1 + ix]
                    if h != 0:
                        merge(e1, h)
                    else:
                        table[4 * e1 + x] = f1
                        table[4 * f1 + ix] = e1

    def scan_and_fill(alpha, w, wi):
        """Scan relator w at coset alpha, defining cosets to close the gap."""
        f = alpha
        i = 0
        b = alpha
        j = len(w) - 1
        while True:
            while i <= j:
                t = table[4 * f + w[i]]
                if t == 0:
                    break
                f = t
                i += 1
            if i > j:
                if f != b:
                    coincidence(f, b)
                return
            while j >= i:
                t = table[4 * b + wi[j]]
                if t == 0:
                    break
                b = t
                j -= 1
            if j < i:
                coincidence(f, b)
                return
            if j == i:
                # the scan closes with a single deduction
                table[4 * f + w[i]] = b
                table[4 * b + wi[i]] = f
                return
            define(f, w[i])

    status = None
    t0 = time.time()
    try:
        n = 1
        grow(1)
        p[1] = 1
        for w, wi in subs:
            scan_and_fill(1, w, wi)

        # HLT main loop, repeated until a whole pass changes nothing.  A pass
        # that defines no coset and finds no coincidence certifies that every
        # live coset has a complete row and every relator closes at it.
        while True:
            n_before = n
            c_before = n_coinc
            alpha = 1
            while alpha <= n:
                if p[alpha] == alpha:
                    for w, wi in rels:
                        scan_and_fill(alpha, w, wi)
                        if p[alpha] != alpha:
                            break
                    if p[alpha] == alpha:
                        for x in range(4):
                            if table[4 * alpha + x] == 0:
                                define(alpha, x)
                alpha += 1
            if n == n_before and n_coinc == c_before:
                break
        status = FINITE
    except _Overflow:
        status = INCONCLUSIVE

    elapsed = time.time() - t0

    index = None
    if status == FINITE:
        index = 0
        for c in range(1, n + 1):
            if p[c] == c:
                index += 1
        if index == 1:
            status = TRIVIAL

    return {
        "status": status,
        "index": index,
        "peak_cosets": n,
        "coincidences": n_coinc,
        "elapsed": elapsed,
        "table": table,
        "p": p,
        "n": n,
        "relators": rels,
        "subgroup": subs,
        "max_cosets": max_cosets,
    }


def verify_table(result):
    """Independent check that a completed table really is a coset table.

    Checks, for every live coset a:
      1. all four generator images are defined and live,
      2. table[a][g] = b implies table[b][g^-1] = a,
      3. every relator, traced from a letter by letter, returns to a.
    Returns (ok, message).
    """
    if result["status"] not in (TRIVIAL, FINITE):
        return False, "table is incomplete (enumeration did not finish)"
    table = result["table"]
    p = result["p"]
    n = result["n"]
    live = [c for c in range(1, n + 1) if p[c] == c]
    if len(live) != result["index"]:
        return False, "live coset count %d != reported index %d" % (
            len(live), result["index"])
    for a in live:
        for g in range(4):
            b = table[4 * a + g]
            if b == 0:
                return False, "coset %d has no image under %s" % (a, LETTERS[g])
            if p[b] != b:
                return False, "coset %d maps to dead coset %d under %s" % (
                    a, b, LETTERS[g])
            if table[4 * b + INV[g]] != a:
                return False, "edge %d -%s-> %d has no inverse edge" % (
                    a, LETTERS[g], b)
    for a in live:
        for w, _wi in result["relators"]:
            c = a
            for g in w:
                c = table[4 * c + g]
            if c != a:
                return False, "relator %s does not close at coset %d (ends %d)" % (
                    unparse_word(w), a, c)
    # every subgroup generator must fix the coset of H itself (coset 1)
    one = 1
    while p[one] != one:
        one = p[one]
    for w, _wi in result.get("subgroup", ()):
        c = one
        for g in w:
            c = table[4 * c + g]
        if c != one:
            return False, "subgroup generator %s does not fix coset H" % (
                unparse_word(w))
    # the action must be transitive from coset 1, else the table is not a
    # coset table at all
    seen = {one}
    stack = [one]
    while stack:
        c = stack.pop()
        for g in range(4):
            d = table[4 * c + g]
            if d not in seen:
                seen.add(d)
                stack.append(d)
    if len(seen) != len(live):
        return False, "action is not transitive (%d of %d cosets reachable)" % (
            len(seen), len(live))
    return True, "table closed and complete on %d cosets" % len(live)


def enumerate_cosets(r1, r2, max_cosets=DEFAULT_MAX_COSETS):
    """Enumerate cosets of H = 1 in <x, y | r1, r2>.

    Returns (status, index_or_None, peak_cosets) where status is one of
    "TRIVIAL", "FINITE", "INCONCLUSIVE".
    """
    res = _enumerate([r1, r2], max_cosets)
    return res["status"], res["index"], res["peak_cosets"]


def enumerate_cosets_multi(relators, max_cosets=DEFAULT_MAX_COSETS,
                           subgroup_gens=()):
    """Same, for any number of relators (and optional subgroup generators)."""
    res = _enumerate(relators, max_cosets, subgroup_gens)
    return res["status"], res["index"], res["peak_cosets"]


def abelianized_order(r1, r2):
    """Order of the abelianization of <x, y | r1, r2>, or None if infinite.

    G^ab = Z^2 / im(M) where M has rows the exponent-sum vectors of the
    relators, so |G^ab| = |det M| (infinite when det M = 0).  This is an
    *independent* check on the enumeration: |G^ab| divides |G| when G is
    finite, and det M = 0 forces G to be infinite (so the enumeration must
    never complete).  In particular index 1 requires |det M| = 1.
    """
    rows = []
    for r in (r1, r2):
        w = parse_word(r)
        ex = sum(1 for g in w if g == 0) - sum(1 for g in w if g == 1)
        ey = sum(1 for g in w if g == 2) - sum(1 for g in w if g == 3)
        rows.append((ex, ey))
    det = rows[0][0] * rows[1][1] - rows[0][1] * rows[1][0]
    return None if det == 0 else abs(det)


def _abelian_check(relators, status, index):
    """Cross-check a result against the abelianization.  (ok, message)."""
    if len(relators) != 2:
        return True, ""
    try:
        d = abelianized_order(relators[0], relators[1])
    except ValueError:
        return True, ""
    if d is None:
        if status != INCONCLUSIVE:
            return False, "det = 0 so G^ab is infinite, but got %s" % status
        return True, "ab: infinite"
    if status == INCONCLUSIVE:
        return True, "ab: |G^ab| = %d" % d
    if index % d != 0:
        return False, "|G^ab| = %d does not divide index %d" % (d, index)
    return True, "ab: |G^ab| = %d divides %d" % (d, index)


# --------------------------------------------------------------------------
# self-test


_KNOWN = [
    # (relators, expected status, expected index or None, note)
    # -- balanced 2-relator cases (the shape the CLI takes) -----------------
    (["x", "y"], TRIVIAL, 1, "trivial presentation"),
    (["x", "yy"], FINITE, 2, "Z2"),
    (["xxxxx", "y"], FINITE, 5, "Z5"),
    (["xxxYY", "xyXY"], INCONCLUSIVE, None,
     "abelian Z (det = 0 because [x,y] has zero exponent sums)"),
    (["xxYY", "xxYXYX"], FINITE, 8, "Q8 = <a,b | a^2 = b^2 = (ab)^2>"),
    (["xxYYY", "xyxyxyXXXX"], FINITE, 24,
     "binary tetrahedral <2,3,3> = SL(2,3)"),
    (["xxYYY", "xyxyxyxyXXXXXX"], FINITE, 48, "binary octahedral <2,3,4>"),
    (["xxYYY", "xyxyxyxyxyXXXXXXXX"], FINITE, 120,
     "binary icosahedral <2,3,5> = SL(2,5)"),
    (["xxYYY", "xyxYXY"], TRIVIAL, 1, "AK(2)"),
    (["xxxYYYY", "xyxYXY"], TRIVIAL, 1, "AK(3)"),
    (["xxyXy", "xyyyyxY"], None, None,
     "known NONTRIVIAL (quotient in S5) - must NOT be index 1"),
    (["xxx", "yy"], INCONCLUSIVE, None, "Z3 * Z2 = PSL(2,Z), infinite"),
    (["xx", "yy"], INCONCLUSIVE, None, "infinite dihedral"),
    (["xyXY", "xxxyyy"], INCONCLUSIVE, None, "Z x Z3, infinite (det = 0)"),
    (["xxYY", "xyxY"], None, None, "small case, no claim"),
    # -- more than 2 relators (internal API only) ---------------------------
    (["xyXY"], INCONCLUSIVE, None, "free abelian Z^2, infinite"),
    (["xx", "yy", "xyxy"], FINITE, 4, "Klein four group"),
    (["xxx", "yy", "xyxy"], FINITE, 6, "S3"),
    (["xxxx", "yy", "xyxy"], FINITE, 8, "dihedral of order 8"),
    (["xxxxx", "yy", "xyxy"], FINITE, 10, "dihedral of order 10"),
    (["xxx", "yyy", "xyXY"], FINITE, 9, "Z3 x Z3"),
    (["xxx", "yy", "xyXY"], FINITE, 6, "Z6"),
    (["xxxxxxxxxxxx", "yyyyyyyy", "xyXY"], FINITE, 96, "Z12 x Z8"),
    (["xxx", "yy", "xyxyxy"], FINITE, 12, "A4 = triangle group (2,3,3)"),
    (["xxxx", "yy", "xyxyxy"], FINITE, 24, "S4 = triangle group (2,3,4)"),
    (["xxxxx", "yy", "xyxyxy"], FINITE, 60, "A5 = triangle group (2,3,5)"),
    (["xxxxxx", "yy", "xyxyxy"], INCONCLUSIVE, None,
     "triangle group (2,3,6), infinite (Euclidean)"),
    (["xxx", "yyy", "xyxyxy"], INCONCLUSIVE, None,
     "triangle group (3,3,3), infinite (Euclidean)"),
    (["yy", "xxxxx", "yxyxyxyx", "yXXyxxyXXyxx"], FINITE, 120,
     "S5, Moore presentation"),
    (["yy", "xxxxxx", "yxyxyxyxyx", "yXXyxxyXXyxx", "yXXXyxxxyXXXyxxx"],
     FINITE, 720, "S6, Moore presentation"),
    (["yy", "xxx", "yxyxyxyxyxyxyx", "YXyxYXyxYXyxYXyx"], FINITE, 168,
     "PSL(2,7) = <a,b | a^2, b^3, (ab)^7, [a,b]^4>"),
]


def selftest(max_cosets=2_000_000, verbose=True):
    failures = 0
    for relators, exp_status, exp_index, note in _KNOWN:
        res = _enumerate(relators, max_cosets)
        if exp_status is None:
            ok = True  # informational only, no known answer asserted
        else:
            ok = res["status"] == exp_status and (
                exp_index is None or res["index"] == exp_index)
        vmsg = ""
        if res["status"] in (TRIVIAL, FINITE):
            vok, vtext = verify_table(res)
            if not vok:
                ok = False
            vmsg = "  [verify: %s]" % vtext
        aok, atext = _abelian_check(relators, res["status"], res["index"])
        if not aok:
            ok = False
            vmsg += "  [ABELIAN CHECK FAILED: %s]" % atext
        if not ok:
            failures += 1
        if verbose:
            got = res["status"]
            if res["index"] is not None:
                got += " index %d" % res["index"]
            if exp_status is None:
                want = "(no claim)"
            else:
                want = exp_status + (" index %d" % exp_index if exp_index else "")
            print("%-4s <x,y | %-40s>  ->  %-18s (want %-18s peak %8d, %6.2fs)  %s%s"
                  % ("ok" if ok else "FAIL", ", ".join(relators), got, want,
                     res["peak_cosets"], res["elapsed"], note, vmsg))
    if verbose:
        print("\n%d/%d self-tests passed" % (len(_KNOWN) - failures, len(_KNOWN)))
    return failures == 0


# --------------------------------------------------------------------------
# randomised differential stress test of the coincidence handling
#
# Every transformation below leaves the presented group unchanged up to
# isomorphism (relator inversion / rotation / conjugation and r1 -> r1*r2 leave
# the normal closure alone; swapping the generators is an automorphism of the
# free group).  They do however drive the enumeration down completely different
# paths, with different coincidence patterns.  So: any two variants that both
# complete must report the same index.  A missed or spurious coincidence would
# almost surely break that agreement.


def _invert(w):
    swap = {"x": "X", "X": "x", "y": "Y", "Y": "y"}
    return "".join(swap[c] for c in reversed(w))


def _rotate(w, k):
    if not w:
        return w
    k %= len(w)
    return w[k:] + w[:k]


def _swap_gens(w):
    swap = {"x": "y", "y": "x", "X": "Y", "Y": "X"}
    return "".join(swap[c] for c in w)


def _flip_x(w):
    swap = {"x": "X", "X": "x", "y": "y", "Y": "Y"}
    return "".join(swap[c] for c in w)


def _mul_reduce(a, b):
    return unparse_word(free_reduce(parse_word(a) + parse_word(b)))


def _random_variant(r1, r2, rng):
    if rng.random() < 0.5:
        r1, r2 = r2, r1
    if rng.random() < 0.5:
        r1 = _invert(r1)
    if rng.random() < 0.5:
        r2 = _invert(r2)
    r1 = _rotate(r1, rng.randrange(max(1, len(r1))))
    r2 = _rotate(r2, rng.randrange(max(1, len(r2))))
    if rng.random() < 0.4:
        g = rng.choice("xXyY")
        r1 = _mul_reduce(_mul_reduce(g, r1), {"x": "X", "X": "x",
                                              "y": "Y", "Y": "y"}[g])
    if rng.random() < 0.4:
        # r1 -> r1 * r2^{+-1}: same normal closure, same group
        r1 = _mul_reduce(r1, r2 if rng.random() < 0.5 else _invert(r2))
        if not r1:
            r1 = "x" if rng.random() < 0.5 else "y"
    if rng.random() < 0.5:
        r1, r2 = _swap_gens(r1), _swap_gens(r2)
    if rng.random() < 0.5:
        r1, r2 = _flip_x(r1), _flip_x(r2)
    return r1, r2


def stresstest(trials=400, variants=6, max_cosets=30000, seed=20260801,
               verbose=True):
    """Random differential test: invariance of the index under re-presentation."""
    import random
    rng = random.Random(seed)
    letters = "xXyY"
    mismatches = 0
    verify_fail = 0
    abel_fail = 0
    completed = 0
    for _ in range(trials):
        r1 = "".join(rng.choice(letters) for _ in range(rng.randrange(1, 8)))
        r2 = "".join(rng.choice(letters) for _ in range(rng.randrange(1, 8)))
        if not free_reduce(parse_word(r1)) or not free_reduce(parse_word(r2)):
            continue
        seen = {}
        for k in range(variants):
            v1, v2 = (r1, r2) if k == 0 else _random_variant(r1, r2, rng)
            res = _enumerate([v1, v2], max_cosets)
            aok, atext = _abelian_check([v1, v2], res["status"], res["index"])
            if not aok:
                abel_fail += 1
                if verbose:
                    print("ABELIAN FAIL %s %s: %s" % (v1, v2, atext))
            if res["status"] in (TRIVIAL, FINITE):
                completed += 1
                vok, vtext = verify_table(res)
                if not vok:
                    verify_fail += 1
                    if verbose:
                        print("VERIFY FAIL %s %s: %s" % (v1, v2, vtext))
                seen.setdefault(res["index"], []).append((v1, v2))
        if len(seen) > 1:
            mismatches += 1
            if verbose:
                print("INDEX MISMATCH for <%s, %s>: %s" % (r1, r2, seen))
    if verbose:
        print("stress: %d presentations x %d variants, %d completed runs; "
              "%d index mismatches, %d verify failures, %d abelian failures"
              % (trials, variants, completed, mismatches, verify_fail,
                 abel_fail))
    return mismatches == 0 and verify_fail == 0 and abel_fail == 0


# --------------------------------------------------------------------------
# CLI


def main(argv):
    args = []
    max_cosets = DEFAULT_MAX_COSETS
    i = 1
    do_selftest = False
    do_stress = False
    while i < len(argv):
        a = argv[i]
        if a == "--stresstest":
            do_stress = True
        elif a == "--max-cosets":
            i += 1
            max_cosets = int(argv[i].replace("_", ""))
        elif a.startswith("--max-cosets="):
            max_cosets = int(a.split("=", 1)[1].replace("_", ""))
        elif a == "--selftest":
            do_selftest = True
        elif a in ("-h", "--help"):
            print(__doc__)
            return 0
        else:
            args.append(a)
        i += 1

    if do_stress:
        return 0 if stresstest() else 1

    if do_selftest:
        return 0 if selftest(max_cosets=max_cosets) else 1

    if not args:
        print("usage: python3 todd_coxeter.py R1 R2 [--max-cosets N]",
              file=sys.stderr)
        return 2

    try:
        res = _enumerate(args, max_cosets)
    except ValueError as exc:
        print("error: %s" % exc, file=sys.stderr)
        return 2
    if res["status"] == INCONCLUSIVE:
        print("INCONCLUSIVE: exceeded %d cosets" % max_cosets)
    elif res["status"] == TRIVIAL:
        print("TRIVIAL: index 1")
    else:
        print("FINITE: index %d" % res["index"])
    if res["status"] in (TRIVIAL, FINITE):
        ok, msg = verify_table(res)
        print("verification: %s (%s)" % ("PASS" if ok else "FAIL", msg))
    if len(args) == 2:
        aok, atext = _abelian_check(args, res["status"], res["index"])
        print("abelianization check: %s (%s)"
              % ("PASS" if aok else "FAIL", atext))
    print("peak cosets: %d" % res["peak_cosets"])
    print("coincidences: %d" % res["coincidences"])
    print("elapsed: %.3f s" % res["elapsed"])
    return 0


if __name__ == "__main__":
    try:
        rc = main(sys.argv)
    except BrokenPipeError:
        rc = 0
    sys.exit(rc)
