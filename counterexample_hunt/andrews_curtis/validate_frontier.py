"""Validate the AC engine against the published length-13 frontier.

Havas & Ramsay (IJAC 2003) proved: every balanced 2-generator presentation of
the TRIVIAL group with total relator length <= 13 is AC-trivializable, except
possibly AK(3) = <x,y | x^3 = y^4, xyx = yxy>, which is AC-equivalent to the
unique remaining class. Miasnikov (IJAC 1999) covers total length <= 12 with
no exceptions.

We cannot decide triviality of a presentation directly, so each enumerated
presentation is resolved one of three ways:

  TRIVIALIZED  the engine found an AC path to <x,y | x,y>. This *proves* the
               presentation defines the trivial group and is AC-trivial.
  NONTRIVIAL   we exhibit a homomorphism onto a nontrivial subgroup of some
               S_n killing both relators. This *proves* the group is not
               trivial, so the presentation is outside the theorem's scope.
  OPEN         neither. Per the literature the only OPEN case with total
               length <= 13 that presents the trivial group should be AK(3);
               other OPEN cases must be nontrivial groups with no small
               permutation quotient.

Enumeration is over cyclically reduced, rotation-canonical relator pairs with
trivial abelianization (|det| = 1), deduplicated modulo the AC-preserving
symmetries: relator inversion, relator swap, and the 8 signed permutations of
the generators. Every identification used is a genuine AC-equivalence, so
coverage of AC-classes is complete.
"""

import itertools
import json
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor

LETTERS = "xXyY"
MAX_TOTAL = 13
ENGINE = "./ac_search"

# escalating (cap, max_states) attempts for each presentation
ATTEMPTS = [(14, 3_000_000), (17, 20_000_000)]


def inv_letter(c):
    return c ^ 1


def freely_reduced_words(length):
    """All freely reduced words of the given length, as int tuples."""
    if length == 0:
        yield ()
        return
    stack = [(c,) for c in range(4)]
    while stack:
        w = stack.pop()
        if len(w) == length:
            yield w
            continue
        for c in range(4):
            if c != inv_letter(w[-1]):
                stack.append(w + (c,))


def is_cyclically_reduced(w):
    return len(w) == 1 or w[-1] != inv_letter(w[0])


def min_rotation(w):
    return min(w[i:] + w[:i] for i in range(len(w)))


def word_inverse(w):
    return tuple(inv_letter(c) for c in reversed(w))


def generator_symmetries():
    """The 8 permutations of {x,X,y,Y} commuting with inversion."""
    syms = []
    for swap in (False, True):
        for fx in (False, True):
            for fy in (False, True):
                p = [0] * 4
                for g, flip in ((0, fx), (1, fy)):
                    tgt = (g ^ 1) if swap else g
                    p[2 * g] = 2 * tgt + (1 if flip else 0)
                    p[2 * g + 1] = 2 * tgt + (0 if flip else 1)
                syms.append(tuple(p))
    return syms


SYMS = generator_symmetries()


def exponent_sums(w):
    a = b = 0
    for c in w:
        if c < 2:
            a += 1 if c == 0 else -1
        else:
            b += 1 if c == 2 else -1
    return a, b


def canonical_pair(r1, r2):
    """Least representative of the pair's class under all AC symmetries."""
    best = None
    for sym in SYMS:
        s1 = tuple(sym[c] for c in r1)
        s2 = tuple(sym[c] for c in r2)
        for a in (s1, word_inverse(s1)):
            ca = min_rotation(a)
            for b in (s2, word_inverse(s2)):
                cb = min_rotation(b)
                key = tuple(sorted((ca, cb), key=lambda w: (len(w), w)))
                if best is None or key < best:
                    best = key
    return best


def enumerate_presentations(max_total):
    """Deduplicated presentations with |det| = 1 and total length <= max_total."""
    by_len = {}
    for length in range(1, max_total):
        words = []
        for w in freely_reduced_words(length):
            if is_cyclically_reduced(w) and min_rotation(w) == w:
                words.append((w, exponent_sums(w)))
        by_len[length] = words
    seen = set()
    for l1 in range(1, max_total):
        for l2 in range(l1, max_total - l1 + 1):
            if l1 + l2 > max_total:
                break
            for w1, (a11, a12) in by_len[l1]:
                for w2, (a21, a22) in by_len[l2]:
                    if abs(a11 * a22 - a12 * a21) != 1:
                        continue
                    seen.add(canonical_pair(w1, w2))
    return sorted(seen, key=lambda p: (len(p[0]) + len(p[1]), p))


def to_str(w):
    return "".join(LETTERS[c] for c in w)


def run_engine(pres):
    """Return (status, moves, detail). status in {TRIVIALIZED, UNRESOLVED}."""
    r1, r2 = to_str(pres[0]), to_str(pres[1])
    last = ""
    for cap, budget in ATTEMPTS:
        out = subprocess.run(
            [ENGINE, r1, r2, "--cap", str(cap), "--max-states", str(budget)],
            stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, text=True,
        ).stdout
        if "TRIVIALIZATION FOUND" in out:
            moves = int(out.split("TRIVIALIZATION FOUND — ")[1].split(" ")[0])
            assert "Replay verified" in out, f"replay failed for {r1} {r2}"
            return "TRIVIALIZED", moves, f"cap={cap}"
        last = "EXHAUSTED" if "EXHAUSTED" in out else "TRUNCATED"
    return "UNRESOLVED", None, last


# --- nontriviality certificates via permutation quotients ---------------------

def compose(p, q):
    return tuple(p[i] for i in q)


def evaluate(word, a, b, n):
    ident = tuple(range(n))
    ainv = tuple(sorted(range(n), key=lambda i: a[i]))
    binv = tuple(sorted(range(n), key=lambda i: b[i]))
    imgs = (a, ainv, b, binv)
    r = ident
    for c in word:
        r = compose(r, imgs[c])
    return r


def nontrivial_quotient(pres, max_n=5):
    """Find (n, a, b) with both relators killed and <a,b> nontrivial."""
    r1, r2 = pres
    for n in range(2, max_n + 1):
        ident = tuple(range(n))
        perms = list(itertools.permutations(range(n)))
        for a in perms:
            for b in perms:
                if a == ident and b == ident:
                    continue
                if evaluate(r1, a, b, n) == ident and evaluate(r2, a, b, n) == ident:
                    return n, a, b
    return None


def main():
    max_total = int(sys.argv[1]) if len(sys.argv) > 1 else MAX_TOTAL
    print(f"enumerating balanced presentations, total length <= {max_total} ...")
    pres = enumerate_presentations(max_total)
    print(f"{len(pres)} distinct AC-classes with trivial abelianization")

    ak3 = canonical_pair(tuple(LETTERS.index(c) for c in "xxxYYYY"),
                         tuple(LETTERS.index(c) for c in "xyxYXY"))
    if max_total >= 13:
        assert ak3 in pres, "AK(3) missing from enumeration — enumeration is wrong"
        print(f"AK(3) present in enumeration as {to_str(ak3[0])} {to_str(ak3[1])}")

    with ThreadPoolExecutor(max_workers=4) as pool:
        results = list(pool.map(run_engine, pres))

    trivialized = [(p, r) for p, r in zip(pres, results) if r[0] == "TRIVIALIZED"]
    unresolved = [p for p, r in zip(pres, results) if r[0] == "UNRESOLVED"]
    print(f"\nTRIVIALIZED: {len(trivialized)}   UNRESOLVED by search: {len(unresolved)}")

    nontrivial, still_open = [], []
    for p in unresolved:
        q = nontrivial_quotient(p)
        (nontrivial if q else still_open).append((p, q))
    print(f"of the unresolved: {len(nontrivial)} proven NONTRIVIAL groups, "
          f"{len(still_open)} OPEN")
    for p, (n, a, b) in nontrivial:
        print(f"  nontrivial: {to_str(p[0]):<14} {to_str(p[1]):<14} "
              f"via S_{n}, x->{a}, y->{b}")

    print("\nOPEN cases (should be AK(3) and nontrivial groups without small quotients):")
    for p, _ in still_open:
        mark = "  <-- AK(3)" if p == ak3 else ""
        print(f"  {to_str(p[0]):<14} {to_str(p[1]):<14} "
              f"total={len(p[0]) + len(p[1])}{mark}")

    max_moves = max((r[1] for _, r in trivialized), default=0)
    print(f"\nlongest trivialization found: {max_moves} moves")

    with open("validation_results.json", "w") as f:
        json.dump({
            "max_total_length": max_total,
            "classes": len(pres),
            "trivialized": len(trivialized),
            "proven_nontrivial": len(nontrivial),
            "open": [[to_str(p[0]), to_str(p[1])] for p, _ in still_open],
            "nontrivial_certificates": [
                {"r1": to_str(p[0]), "r2": to_str(p[1]),
                 "n": n, "x_image": list(a), "y_image": list(b)}
                for p, (n, a, b) in nontrivial
            ],
            "ak3_canonical": [to_str(ak3[0]), to_str(ak3[1])],
            "longest_trivialization_moves": max_moves,
            "attempts": ATTEMPTS,
        }, f, indent=2)


if __name__ == "__main__":
    main()
