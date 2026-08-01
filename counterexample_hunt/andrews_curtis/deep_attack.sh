#!/bin/sh
# Adversarial pass: try much harder to trivialize each surviving candidate
# than the validation sweep did (which escalated only to cap 17 / 20M states).
cd "$(dirname "$0")" || exit 1
while read -r r1 r2; do
  printf "%-10s %-12s " "$r1" "$r2"
  ./ac_search "$r1" "$r2" --cap 24 --max-states 40000000 2>/dev/null \
    | grep -oE "TRIVIALIZATION FOUND — [0-9]+ moves|EXHAUSTED[^:]*|TRUNCATED[^:]*" | head -1
done < "${1:-targets_len14_remaining.txt}"
