/*
 * frontier_sweep — Collatz convergence check for odd numbers immediately
 * above the exhaustively verified frontier (2^71, Barina 2025).
 *
 * For each odd n in [2^71 + 1, 2^71 + 2*COUNT), iterate the map until the
 * value drops below 2^71; every value below the frontier is already verified
 * to reach 1, so dropping below it proves convergence of n.
 *
 * A divergent orbit would hit the step or overflow guard; none is expected
 * (this is a demonstration probe, not meaningful progress on the conjecture).
 *
 * Usage: ./frontier_sweep COUNT
 */
#include <inttypes.h>
#include <stdio.h>
#include <stdlib.h>

typedef unsigned __int128 u128;

int main(int argc, char **argv) {
    uint64_t count = argc > 1 ? strtoull(argv[1], NULL, 10) : 100000000ULL;
    const u128 frontier = (u128)1 << 71;
    const u128 guard = (u128)1 << 120;
    const uint64_t step_cap = 1000000;

    uint64_t max_steps = 0, max_steps_off = 0;
    u128 max_excursion = 0;
    uint64_t max_exc_off = 0;

    for (uint64_t i = 0; i < count; i++) {
        u128 v = frontier + 1 + 2 * (u128)i;
        uint64_t steps = 0;
        while (v >= frontier) {
            if (v & 1) {
                if (v > guard) { printf("OVERFLOW GUARD HIT at offset %" PRIu64 "\n", i); return 1; }
                v = 3 * v + 1;
            } else {
                v >>= 1;
            }
            if (++steps > step_cap) { printf("STEP CAP HIT at offset %" PRIu64 "\n", i); return 1; }
            if (v > max_excursion) { max_excursion = v; max_exc_off = i; }
        }
        if (steps > max_steps) { max_steps = steps; max_steps_off = i; }
    }
    printf("verified convergence (drop below 2^71) for %" PRIu64 " odd numbers above 2^71\n", count);
    printf("max steps to drop below frontier: %" PRIu64 " (n = 2^71 + %" PRIu64 ")\n",
           max_steps, 1 + 2 * max_steps_off);
    printf("max excursion: %.6g (n = 2^71 + %" PRIu64 ")\n",
           (double)max_excursion, 1 + 2 * max_exc_off);
    printf("no divergence, no cycle candidates found (as expected)\n");
    return 0;
}
