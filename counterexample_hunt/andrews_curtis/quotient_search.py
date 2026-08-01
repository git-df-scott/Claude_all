"""Search for a nontrivial permutation quotient of a given presentation.

Extends the S_5 check in validate_frontier.py to larger symmetric groups.
Finding generators a, b in S_n that kill both relators, with <a,b> nontrivial,
proves the presented group is not the trivial group.

Usage: python3 quotient_search.py R1 R2 [max_n]
"""

import itertools
import sys

LETTERS = "xXyY"


def evaluate(word, imgs, ident):
    r = ident
    for c in word:
        p = imgs[c]
        r = tuple(p[i] for i in r)
    return r


def search(r1, r2, n):
    ident = tuple(range(n))
    perms = list(itertools.permutations(range(n)))
    invs = {p: tuple(sorted(range(n), key=lambda i: p[i])) for p in perms}
    for a in perms:
        ai = invs[a]
        for b in perms:
            if a == ident and b == ident:
                continue
            imgs = (a, ai, b, invs[b])
            if evaluate(r1, imgs, ident) != ident:
                continue
            if evaluate(r2, imgs, ident) == ident:
                return a, b
    return None


def main():
    r1 = tuple(LETTERS.index(c) for c in sys.argv[1])
    r2 = tuple(LETTERS.index(c) for c in sys.argv[2])
    max_n = int(sys.argv[3]) if len(sys.argv) > 3 else 7
    for n in range(2, max_n + 1):
        hit = search(r1, r2, n)
        print(f"S_{n}: {'FOUND ' + str(hit) if hit else 'none'}", flush=True)
        if hit:
            print(f"NONTRIVIAL: {sys.argv[1]} {sys.argv[2]} surjects onto a "
                  f"nontrivial subgroup of S_{n}")
            return
    print(f"no nontrivial quotient in S_n for n <= {max_n}")


if __name__ == "__main__":
    main()
