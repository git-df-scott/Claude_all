# HANDOFF — Euler sixth-power counterexample hunt

**Read this first. It contains everything needed to continue the project with
zero context loss.** Written 2026-07-31 by the session that built this
directory. The working branch is `claude/euler-sixth-power-d4ywmn`; the open
draft PR is https://github.com/git-df-scott/Claude_all/pull/3.

## 1. Mission

Find positive integers with `a^6+b^6+c^6+d^6+e^6 = f^6` — a counterexample to
Euler's sum of powers conjecture (1769) for k=6, open for 257 years. The k=5
case fell in 1966 (Lander–Parkin, 144^5), k=4 in 1988 (Elkies/Frye). For k=6
nothing is known; published exhaustive searches (ca. 2002) cleared f ≤ 730,000.

## 2. What has been established (do not re-derive)

### Structure theorem (proved, load-bearing)
Sixth powers are ≡ 0 or 1 mod 7, mod 9, and mod 8. Counting terms mod each
prime shows: in any **primitive** solution, exactly one of a..e is odd,
exactly one is coprime to 3, exactly one is coprime to 7, and **f is coprime
to 42**. Searching primitives is fully general (any solution is a multiple of
a primitive one with smaller f).

### The concentrated ("case A") reduction
If all three exempt slots land on ONE term t (t coprime to 42), the other four
terms are divisible by 42, forcing `f^6 ≡ t^6 (mod 42^6)`, `42^6 =
5,489,031,744`. The valid t are exactly `u·f mod 42^6` where u ranges over the
**144 sixth roots of unity mod 42^6** (4 mod 2^6 × 6 mod 3^6 × 6 mod 7^6, by
CRT — verified by direct powering in code). Expected candidates (f,t) with
f ≤ F: `(2/7)·144·F²/(2·42^6)` — sparse enough to sweep millions deep.
Each candidate needs `m = (f^6−t^6)/42^6` written as a sum of four sixth
powers with bases ≤ (f−1)/42; in that reduced equation the number of odd bases
is EXACTLY m mod 8, coprime-to-3 bases EXACTLY m mod 9, coprime-to-7 EXACTLY
m mod 7 (each ≤ 4 or infeasible).

The complementary cases (exemptions spread over 2–3 terms) have candidate
density ~F³/42^6 or worse — NOT thinnable this way; they remain covered only
to 730k. Case A carries ~1/25 of heuristic solution mass. Expected hits for
the whole 730k→4M sweep: ~1–5%. That is the honest number.

### Results so far (recorded; /tmp logs die with the container — trust THIS table)

| f range | case | candidates | solutions |
|---|---|---|---|
| ≤ 5,000 | all (exhaustive, search.c) | — | 0 |
| 700,000–730,000 | concentrated | 124 | 0 (consistency check) |
| 730,000–1,000,000 | concentrated | 1,314 | 0 |
| 1,000,000–1,500,000 | concentrated | 3,523 | 0 |
| 1,500,000–2,000,000 | concentrated | 5,129 | 0 |
| 2,000,000–2,500,000 | concentrated | 7,222 | 0 |
| 2,500,000–3,200,000 | concentrated | 12,443 | 0 |
| 3,200,000–4,000,000 | concentrated | 18,003 | 0 |
| 4,000,000–4,300,000 | concentrated | 8,050 | 0 |

## 3. The code (all in euler6/, all committed)

- `search.c` — exhaustive (6,1,5) DFS with class strides + mod 13/43 sieves.
  VALIDATED: `./search 5 2 150` re-finds Lander–Parkin's 144^5 in ~30 ms.
  Scales ~f^3.3; useless past ~20k on small hardware; kept for verification.
- `caseA2.c` — THE workhorse. 144-root candidate enumeration; exact class
  budgets with forced strides (up to 42); residue masks mod 64/27/49/13/43
  (incremental); blocked Bloom filter over all pair sums x^6+y^6 (8 hash bits
  in one 64-byte line, ~16 bits/pair default, optional 3rd argv = bits/pair);
  Bloom positives verified exactly (false positives cost time, never
  correctness). Usage: `./caseA2 <fmin> <fmax> [bits_per_pair]`.
- `caseA.c` — slower predecessor, provenance only. `audit.c` — completeness
  audit; brute force matches the 144-root enumeration exactly on sampled f
  (run it after any change to root generation).
- Build everything: `gcc -O3 -march=native -fopenmp -o X X.c -lm`.

### Hard limits to respect
- `fmax ≤ 1e8` (128-bit overflow guard, enforced in code).
- Bloom RAM ≈ `(fmax/42)² / 2 × bits_per_pair / 8` bytes. 15 GB box → fmax
  ≈ 4M at 16 bpp. 32 GB → ~5.5M. 64 GB → ~8M (12 bpp → ~10M, fp ~10⁻³, fine).
- Don't drop below ~12 bpp; exact-verify cost explodes.

## 4. Immediate next actions for the successor session

1. `git clone` / checkout branch `claude/euler-sixth-power-d4ywmn`, build
   caseA2, run the two sanity checks (§3 search.c validation; `./caseA2
   700000 730000` → 124 candidates, 0 found, <1 min).
2. **Re-run `./caseA2 2500000 3200000`** unless the PR/README already records
   its completion (the predecessor session may have died mid-chunk; ranges are
   idempotent, re-running is safe and takes ~1–2 h on 4 cores).
3. Run `./caseA2 3200000 4000000 12` (12 bpp fits 15 GB).
4. After each chunk: update the results tables in README.md AND HANDOFF.md,
   commit, push, and edit the PR #3 body. THE REPO IS THE ONLY DURABLE MEMORY.
5. Background long runs with the framework's run_in_background (they survive
   past tool timeouts); schedule an hourly send_later check-in to harvest
   results, then re-arm. pkill kills framework background tasks — don't.

## 5. The user's PC campaign (agreed division of labor)

Cloud territory **730k → 4.3M is COMPLETE** (55,684 candidates, 0 solutions). The user's PC owns **4.3M upward** (starting at 4M is harmless overlap) (more cores
+ RAM). Recipe already given to the user: clone branch, build, sanity-check,
then `./caseA2 4000000 5000000` etc. Coordinate by claiming ranges in the PR.
If the user asks, build these two upgrades (a few hours work, high value for
a multi-week campaign):
- checkpointing (periodically write max-f-completed; resume flag);
- disk-backed sorted pair-sum table as Bloom alternative (removes RAM
  ceiling, ~30% slower).

## 6. If a SOLUTION line ever prints

1. DO NOT announce. Verify with independent exact arithmetic (python:
   `a**6+b**6+c**6+d**6+e**6 == f**6`), from the printed parts (parts42 are
   already multiplied by 42; t is the fifth term; f is f).
2. Check gcd/primitivity, sanity-check each part < f.
3. If it verifies: commit the tuple to the repo immediately (durable
   timestamped priority), update PR, tell the user with the verification
   transcript. It would be the first (6,1,5) solution in history — treat
   with according care.

## 7. Known heuristics & honest framing (for continuity of judgment)

- Expected solutions grow ~C·log F; C is small (nothing to 730k while k=5
  solved at f=144). Each e-fold of height costs ~F⁴ compute. Brute force
  alone will not win; this sweep buys real but small probability.
- The genuinely-new-math path: the variety is a **Calabi–Yau fourfold** (sextic
  in P⁵, trivial canonical bundle) — the same structural class as Elkies' K3.
  An elliptic fibration over ℚ with infinite-order section, or a ℚ-rational
  curve with all coordinates nonzero, would beat any amount of CPU. Shioda's
  inductive structure only works over cyclotomic fields (and FLT n=6 blocks
  descent from the Fermat sextic curve). This is the research direction if
  the user wants "more math, less compute."
- Precedent (July 2026): Fable-5-assisted counterexample to the Jacobian
  conjecture (n ≥ 3, Alpoge announcement, not yet peer-reviewed). Construction
  beat enumeration there; same lesson.

## 8. Environment gotchas (cloud session specifics)

- Default branch of the repo is `claude/access-skills-md-file-iv9090` (PR #3
  targets it; there is no main/master).
- Commit as `-c user.email=gerrykelvin22@yahoo.com -c user.name="Claude"`;
  end commit messages with the Claude Code co-author footer.
- Binaries/logs are gitignored (euler6/.gitignore, paths relative to euler6/).
- No `gh` CLI — use GitHub MCP tools for PR edits.
- /tmp task outputs are LOST on container death — flush results to the repo
  early and often.
