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


def _prepare(relators):
    """Normalise relator strings into (forward, inverse-letter) index tuples."""
    rels = []
    for r in relators:
        w = cyclic_reduce(parse_word(r))
        if not w:
            continue  # freely trivial relator, no information
        rels.append((tuple(w), tuple(INV[g] for g in w)))
    return rels


def _enumerate(relators, max_cosets, subgroup_gens=()):
    """Core enumeration.  Returns a dict with the raw table and statistics."""
    rels = _prepare(relators)
    subs = _prepare(subgroup_gens)

    if max_cosets < 1:
        raise ValueError("max_cosets must be >= 1")

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

    def scan_only(alpha, w, wi):
        """Scan without defining.  True if the scan closed on alpha."""
        f = alpha
        i = 0
        b = alpha
        j = len(w) - 1
        while i <= j:
            t = table[4 * f + w[i]]
            if t == 0:
                break
            f = t
            i += 1
        if i > j:
            return f == b
        while j >= i:
            t = table[4 * b + wi[j]]
            if t == 0:
                break
            b = t
            j -= 1
        return False

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


# --------------------------------------------------------------------------
# self-test


_KNOWN = [
    # (relators, expected status, expected index or None, note)
    # -- 2-relator cases (the shape the CLI takes) --------------------------
    (["x", "y"], TRIVIAL, 1, "trivial presentation"),
    (["x", "yy"], FINITE, 2, "Z2"),
    (["xxxxx", "y"], FINITE, 5, "Z5"),
    (["xyXY", "xxxyyy"], FINITE, 3, "Z3 (abelian, det = 3)"),
    (["xxYYY", "xxYXYXYX"], FINITE, 24, "binary tetrahedral <2,3,3> = SL(2,3)"),
    (["xxYYY", "xxYXYXYXYX"], FINITE, 48, "binary octahedral <2,3,4>"),
    (["xxYYY", "xxYXYXYXYXYX"], FINITE, 120, "binary icosahedral <2,3,5> = SL(2,5)"),
    (["xxYYY", "xyxYXY"], TRIVIAL, 1, "AK(2)"),
    (["xxxYYYY", "xyxYXY"], TRIVIAL, 1, "AK(3)"),
    (["xxx", "yy"], INCONCLUSIVE, None, "Z3 * Z2 = PSL(2,Z), infinite"),
    (["xx", "yy"], INCONCLUSIVE, None, "infinite dihedral"),
    (["xxYY", "xyxY"], None, None, "small case, no claim"),
    # -- more than 2 relators (internal API only) ---------------------------
    (["xyXY"], INCONCLUSIVE, None, "free abelian Z^2, infinite"),
    (["xx", "yy", "xyxy"], FINITE, 4, "Klein four group"),
    (["xxx", "yy", "xyxy"], FINITE, 6, "S3"),
    (["xxxx", "yy", "xyxy"], FINITE, 8, "dihedral of order 8"),
    (["xxxxx", "yy", "xyxy"], FINITE, 10, "dihedral of order 10"),
    (["xxx", "yyy", "xyXY"], FINITE, 9, "Z3 x Z3"),
    (["xxx", "yy", "xyxyxy"], FINITE, 12, "A4 = triangle group (2,3,3)"),
    (["xxxx", "yy", "xyxyxy"], FINITE, 24, "S4 = triangle group (2,3,4)"),
    (["xxxxx", "yy", "xyxyxy"], FINITE, 60, "A5 = triangle group (2,3,5)"),
    (["xxxxxx", "yy", "xyxyxy"], INCONCLUSIVE, None,
     "triangle group (2,3,6), infinite (Euclidean)"),
    (["xxx", "yyy", "xyxyxy"], INCONCLUSIVE, None,
     "triangle group (3,3,3), infinite (Euclidean)"),
    (["xxxxx", "yy", "xyxyxyxy", "xyxyXYXY"], FINITE, 120, "S5"),
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
# CLI


def main(argv):
    args = []
    max_cosets = DEFAULT_MAX_COSETS
    i = 1
    do_selftest = False
    while i < len(argv):
        a = argv[i]
        if a == "--max-cosets":
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

    if do_selftest:
        return 0 if selftest(max_cosets=max_cosets) else 1

    if not args:
        print("usage: python3 todd_coxeter.py R1 R2 [--max-cosets N]",
              file=sys.stderr)
        return 2

    res = _enumerate(args, max_cosets)
    if res["status"] == INCONCLUSIVE:
        print("INCONCLUSIVE: exceeded %d cosets" % max_cosets)
    elif res["status"] == TRIVIAL:
        print("TRIVIAL: index 1")
    else:
        print("FINITE: index %d" % res["index"])
    if res["status"] in (TRIVIAL, FINITE):
        ok, msg = verify_table(res)
        print("verification: %s (%s)" % ("PASS" if ok else "FAIL", msg))
    print("peak cosets: %d" % res["peak_cosets"])
    print("coincidences: %d" % res["coincidences"])
    print("elapsed: %.3f s" % res["elapsed"])
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
