#!/bin/sh
# Partition the unclassified candidates into genuine AC-classes: repeatedly
# take the first one still unclassified, exhaustively enumerate its cap-12
# component with every other as a simultaneous target, and remove the whole
# class it captures. Exhaustion makes each absence rigorous.
cd "$(dirname "$0")" || exit 1
pool=partition_pool.txt
for i in $(seq 1 12); do
  n=$(grep -c . "$pool") || break
  [ "$n" -le 1 ] && { echo "pool exhausted ($n left)"; break; }
  rep=$(head -1 "$pool"); tail -n +2 "$pool" > .rest.txt
  printf "class %d: rep=%s  pool=%d -> " "$i" "$rep" "$n"
  ./ac_search $rep --cap 12 --max-states 150000000 --targets .rest.txt > .part.log 2>&1
  hits=$(grep -c "TARGET IN COMPONENT" .part.log)
  status=$(grep -oE "EXHAUSTED|TRUNCATED" .part.log | head -1)
  echo "captured $hits siblings ($status)"
  grep "TARGET IN COMPONENT" .part.log | sed 's/TARGET IN COMPONENT: //; s/ (.*//' > .hits.txt
  grep -vxF -f .hits.txt .rest.txt > .next.txt 2>/dev/null || cp .rest.txt .next.txt
  mv .next.txt "$pool"
done
echo "final unclassified: $(grep -c . $pool)"
