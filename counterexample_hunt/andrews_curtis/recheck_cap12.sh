#!/bin/sh
# The escalate-the-cap strategy was backwards: a LOWER cap prunes the state
# space and concentrates the greedy search, so it can find trivializations a
# higher cap misses. Re-attack at cap 12 with a large state budget.
cd "$(dirname "$0")" || exit 1
while read -r r1 r2; do
  printf "%-10s %-12s " "$r1" "$r2"
  ./ac_search "$r1" "$r2" --cap 12 --max-states 120000000 2>/dev/null \
    | grep -oE "TRIVIALIZATION FOUND — [0-9]+ moves|EXHAUSTED[^:]*|TRUNCATED[^:]*" | head -1
done < "$1"
