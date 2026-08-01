"""Triage the OPEN classes from a validate_frontier.py run.

Each open class is one of three things, and we try to decide which:

  AK(3)-equivalent  a path to AK(3) exists within the cap, so it is the known
                    candidate in disguise, not a new one.
  NONTRIVIAL        a nontrivial permutation quotient exists (searched further
                    than validate_frontier.py does), so it is not a
                    presentation of the trivial group at all.
  NEW CANDIDATE     neither -- a genuinely new minimal candidate counterexample
                    at this length, or a nontrivial group with no small quotient.

Usage: python3 triage_open.py validation_results_len14.json [cap] [max_states]
"""

import json
import subprocess
import sys

AK3 = ("xyxYXY", "xxxxYYY")


def target_test(r1, r2, cap, budget):
    out = subprocess.run(
        ["./ac_search", r1, r2, "--cap", str(cap), "--max-states", str(budget),
         "--target", AK3[0], AK3[1]],
        stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, text=True,
    ).stdout
    if "TARGET REACHED" in out:
        return "AK3_EQUIVALENT"
    if "EXHAUSTED" in out:
        return "NOT_AK3_within_cap"
    return "INCONCLUSIVE"


def quotient_test(r1, r2, max_n=6):
    out = subprocess.run(
        ["python3", "quotient_search.py", r1, r2, str(max_n)],
        stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, text=True,
    ).stdout
    return "NONTRIVIAL" if "NONTRIVIAL:" in out else "no_small_quotient"


def main():
    path = sys.argv[1]
    cap = int(sys.argv[2]) if len(sys.argv) > 2 else 14
    budget = int(sys.argv[3]) if len(sys.argv) > 3 else 60_000_000
    data = json.load(open(path))
    openc = [tuple(p) for p in data["open"]]
    print(f"{len(openc)} open classes from {path}; cap={cap} budget={budget}\n")

    verdicts = []
    for r1, r2 in openc:
        if (r1, r2) == AK3:
            verdict, q, t = "AK(3) ITSELF", "-", "-"
        else:
            q = quotient_test(r1, r2)
            t = target_test(r1, r2, cap, budget)
            if q == "NONTRIVIAL":
                verdict = "NOT A TRIVIAL-GROUP PRESENTATION"
            elif t == "AK3_EQUIVALENT":
                verdict = "AK(3) in disguise"
            elif t == "NOT_AK3_within_cap":
                verdict = "NEW CANDIDATE (not AK(3) within cap)"
            else:
                verdict = "UNDECIDED (search budget)"
        verdicts.append({"r1": r1, "r2": r2, "quotient": q,
                         "ak3_test": t, "verdict": verdict})
        print(f"  {r1:<16} {r2:<16} {verdict}")

    with open("triage_open_results.json", "w") as f:
        json.dump({"source": path, "cap": cap, "max_states": budget,
                   "results": verdicts}, f, indent=2)


if __name__ == "__main__":
    main()
