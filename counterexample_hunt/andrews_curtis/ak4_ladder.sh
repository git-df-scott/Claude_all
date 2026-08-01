#!/bin/sh
# AK(4)'s cap-12 component is small (5.4M), so exhaustion can be pushed much
# further in cap than AK(3) allows. Climb until the budget stops us.
cd "$(dirname "$0")" || exit 1
for cap in 13 14 15 16 17; do
  printf "cap %2d: " "$cap"
  ./ac_search "xxxxYYYYY" "xyxYXY" --cap "$cap" --max-states 120000000 2>/dev/null \
    | grep -oE "TRIVIALIZATION FOUND — [0-9]+ moves|(EXHAUSTED|TRUNCATED)[^:]*: [^.]*\. expanded=[0-9]+ stored=[0-9]+ best_total_length=[0-9]+"
done
