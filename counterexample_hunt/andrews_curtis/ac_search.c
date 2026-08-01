/*
 * ac_search — bounded exhaustive search over Andrews–Curtis moves for
 * balanced 2-generator presentations.
 *
 * State: a pair of freely reduced words (r1, r2) over {x, X, y, Y}
 * (X = x^-1, Y = y^-1), canonicalized by sorting the pair.
 *
 * Elementary moves (all are classical AC moves, symmetric in r1/r2):
 *   r_i -> r_i * r_j          r_i -> r_i * r_j^-1     (i != j)
 *   r_i -> r_i^-1
 *   r_i -> g * r_i * g^-1     for each generator letter g in {x, X, y, Y}
 *
 * A state is expanded only if both relators stay within --cap letters.
 * If the search exhausts the queue without reaching total length 2, then
 * NO sequence of these elementary moves trivializes the input while every
 * intermediate relator has length <= cap. (Note: this is stated for
 * single-letter conjugations; a classical conjugation by a long word is a
 * chain of these and its intermediate stages must also respect the cap.)
 *
 * Success = reaching total length 2. AC moves preserve the normal closure,
 * so from a presentation of the trivial group the only length-2 states are
 * {x^±1, y^±1}, i.e. the trivial presentation up to inversion.
 *
 * With --target T1 T2 the search instead reports whether the target
 * presentation is reachable, which decides AC-equivalence of the two inputs
 * within the cap: every move's inverse is also a move and respects the cap,
 * so the cap-restricted move graph is undirected and the reachable set is
 * exactly the connected component of the start.
 *
 * Usage: ./ac_search R1 R2 --cap N [--max-states N] [--target T1 T2]
 */
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAXCAP 60
#define WBUF (2 * MAXCAP + 4)

static const char LET[4] = {'x', 'X', 'y', 'Y'};
static int inv(int c) { return c ^ 1; }

typedef struct { uint8_t s[WBUF]; int n; } Word;

/* stack-based free reduction of the concatenation of two buffers */
static void reduce_cat(const uint8_t *a, int na, const uint8_t *b, int nb, Word *out) {
    int n = 0;
    for (int i = 0; i < na + nb; i++) {
        int c = i < na ? a[i] : b[i - na];
        if (n > 0 && out->s[n - 1] == (uint8_t)inv(c)) n--;
        else out->s[n++] = (uint8_t)c;
    }
    out->n = n;
}

static void word_inv(const Word *w, Word *out) {
    for (int i = 0; i < w->n; i++) out->s[i] = (uint8_t)inv(w->s[w->n - 1 - i]);
    out->n = w->n;
}

static void conj_gen(const Word *w, int g, Word *out) {
    uint8_t buf[WBUF];
    buf[0] = (uint8_t)g;
    memcpy(buf + 1, w->s, (size_t)w->n);
    buf[w->n + 1] = (uint8_t)inv(g);
    Word tmp;
    reduce_cat(buf, w->n + 2, NULL, 0, &tmp);
    *out = tmp;
}

static int word_cmp(const Word *a, const Word *b) {
    if (a->n != b->n) return a->n - b->n;
    return memcmp(a->s, b->s, (size_t)a->n);
}

/* packed state: [len1][len2][r1 bytes][r2 bytes], r1 <= r2 */
static int pack(const Word *r1, const Word *r2, uint8_t *out) {
    const Word *a = r1, *b = r2;
    if (word_cmp(a, b) > 0) { a = r2; b = r1; }
    out[0] = (uint8_t)a->n;
    out[1] = (uint8_t)b->n;
    memcpy(out + 2, a->s, (size_t)a->n);
    memcpy(out + 2 + a->n, b->s, (size_t)b->n);
    return 2 + a->n + b->n;
}

static void unpack(const uint8_t *p, Word *r1, Word *r2) {
    r1->n = p[0];
    r2->n = p[1];
    memcpy(r1->s, p + 2, (size_t)r1->n);
    memcpy(r2->s, p + 2 + r1->n, (size_t)r2->n);
}

static uint64_t fnv1a(const uint8_t *p, int n) {
    uint64_t h = 1469598103934665603ULL;
    for (int i = 0; i < n; i++) { h ^= p[i]; h *= 1099511628211ULL; }
    return h;
}

/* node storage (SoA) + arena of packed states */
static uint64_t *node_off;
static uint32_t *node_parent;
static uint8_t  *node_move;
static uint8_t  *arena;
static uint64_t arena_used, arena_capacity;
static uint32_t n_nodes;
static uint64_t max_states;

/* hash table: open addressing, stores node index+1 (0 = empty) */
static uint32_t *table;
static uint64_t table_mask;

static uint32_t lookup_or_insert(const uint8_t *key, int klen, uint32_t parent,
                                 uint8_t move, int *is_new) {
    uint64_t h = fnv1a(key, klen) & table_mask;
    for (;;) {
        uint32_t slot = table[h];
        if (slot == 0) {
            if (n_nodes >= max_states) { *is_new = -1; return 0; }
            if (arena_used + (uint64_t)klen > arena_capacity) {
                arena_capacity *= 2;
                arena = realloc(arena, arena_capacity);
                if (!arena) { fprintf(stderr, "arena OOM\n"); exit(1); }
            }
            memcpy(arena + arena_used, key, (size_t)klen);
            node_off[n_nodes] = arena_used;
            node_parent[n_nodes] = parent;
            node_move[n_nodes] = move;
            arena_used += (uint64_t)klen;
            table[h] = n_nodes + 1;
            *is_new = 1;
            return n_nodes++;
        }
        const uint8_t *q = arena + node_off[slot - 1];
        int qlen = 2 + q[0] + q[1];
        if (qlen == klen && memcmp(q, key, (size_t)klen) == 0) { *is_new = 0; return slot - 1; }
        h = (h + 1) & table_mask;
    }
}

/* bucket queue by total relator length */
typedef struct { uint32_t *v; uint64_t n, cap, head; } Bucket;
static Bucket buckets[2 * MAXCAP + 3];

static void bucket_push(int tlen, uint32_t node) {
    Bucket *b = &buckets[tlen];
    if (b->n == b->cap) {
        b->cap = b->cap ? b->cap * 2 : 1024;
        b->v = realloc(b->v, b->cap * sizeof(uint32_t));
        if (!b->v) { fprintf(stderr, "bucket OOM\n"); exit(1); }
    }
    b->v[b->n++] = node;
}

static const char *MOVE_NAMES[] = {
    "r1 <- r1*r2", "r1 <- r1*r2^-1", "r2 <- r2*r1", "r2 <- r2*r1^-1",
    "r1 <- r1^-1", "r2 <- r2^-1",
    "r1 <- x r1 x^-1", "r1 <- X r1 X^-1", "r1 <- y r1 y^-1", "r1 <- Y r1 Y^-1",
    "r2 <- x r2 x^-1", "r2 <- X r2 X^-1", "r2 <- y r2 y^-1", "r2 <- Y r2 Y^-1",
};

/* apply move m to (r1, r2) in place; returns 0 if result is invalid (empty relator) */
static int apply_move(int m, Word *r1, Word *r2, int cap) {
    Word t;
    switch (m) {
        case 0: reduce_cat(r1->s, r1->n, r2->s, r2->n, &t); *r1 = t; break;
        case 1: word_inv(r2, &t);
                { Word u; reduce_cat(r1->s, r1->n, t.s, t.n, &u); *r1 = u; } break;
        case 2: reduce_cat(r2->s, r2->n, r1->s, r1->n, &t); *r2 = t; break;
        case 3: word_inv(r1, &t);
                { Word u; reduce_cat(r2->s, r2->n, t.s, t.n, &u); *r2 = u; } break;
        case 4: word_inv(r1, &t); *r1 = t; break;
        case 5: word_inv(r2, &t); *r2 = t; break;
        default:
            if (m < 10) { conj_gen(r1, m - 6, &t); *r1 = t; }
            else        { conj_gen(r2, m - 10, &t); *r2 = t; }
    }
    if (r1->n == 0 || r2->n == 0) return 0;
    if (r1->n > cap || r2->n > cap) return 0;
    return 1;
}

static void print_word(const Word *w) {
    for (int i = 0; i < w->n; i++) putchar(LET[w->s[i]]);
}

static int parse_word(const char *s, Word *w) {
    w->n = 0;
    for (; *s; s++) {
        const char *p = memchr(LET, *s, 4);
        if (!p) return 0;
        w->s[w->n++] = (uint8_t)(p - LET);
    }
    Word t;
    reduce_cat(w->s, w->n, NULL, 0, &t);
    *w = t;
    return 1;
}

/* replay the move path from the root and verify it reaches the trivial pair */
static void replay_and_print(uint32_t node, const Word *init1, const Word *init2, int cap) {
    uint8_t *path;
    int plen = 0;
    for (uint32_t n = node; node_parent[n] != UINT32_MAX; n = node_parent[n]) plen++;
    path = malloc((size_t)plen);
    int i = plen;
    for (uint32_t n = node; node_parent[n] != UINT32_MAX; n = node_parent[n]) path[--i] = node_move[n];

    Word r1 = *init1, r2 = *init2;
    printf("TRIVIALIZATION FOUND — %d moves. Replay:\n  start: ", plen);
    print_word(&r1); printf("  "); print_word(&r2); printf("\n");
    for (i = 0; i < plen; i++) {
        /* replay must mirror the search exactly: moves were applied to the
           canonically ordered pair */
        if (word_cmp(&r1, &r2) > 0) { Word t = r1; r1 = r2; r2 = t; }
        if (!apply_move(path[i], &r1, &r2, cap)) {
            printf("REPLAY FAILED at step %d — bug in search\n", i);
            exit(1);
        }
        printf("  %-18s ", MOVE_NAMES[path[i]]);
        print_word(&r1); printf("  "); print_word(&r2); printf("\n");
    }
    if (r1.n + r2.n != 2 || r1.s[0] / 2 == r2.s[0] / 2) {
        printf("REPLAY FAILED: final state not trivial pair — bug in search\n");
        exit(1);
    }
    printf("Replay verified: reached trivial presentation.\n");
    free(path);
}

int main(int argc, char **argv) {
    if (argc < 4) {
        fprintf(stderr, "usage: %s R1 R2 --cap N [--max-states N]\n", argv[0]);
        return 2;
    }
    Word init1, init2;
    if (!parse_word(argv[1], &init1) || !parse_word(argv[2], &init2)) {
        fprintf(stderr, "bad word (letters must be x X y Y)\n");
        return 2;
    }
    int cap = 0;
    max_states = 50000000ULL;
    uint8_t target_key[WBUF * 2];
    int target_len = 0;
    for (int i = 3; i < argc - 1; i++) {
        if (!strcmp(argv[i], "--cap")) cap = atoi(argv[i + 1]);
        if (!strcmp(argv[i], "--max-states")) max_states = strtoull(argv[i + 1], NULL, 10);
        if (!strcmp(argv[i], "--target") && i + 2 < argc) {
            Word t1, t2;
            if (!parse_word(argv[i + 1], &t1) || !parse_word(argv[i + 2], &t2)) {
                fprintf(stderr, "bad target word\n");
                return 2;
            }
            target_len = pack(&t1, &t2, target_key);
        }
    }
    if (cap < 1 || cap > MAXCAP) { fprintf(stderr, "cap must be 1..%d\n", MAXCAP); return 2; }

    /* abelianization sanity check: exponent-sum matrix must have det ±1 */
    long ex[2][2] = {{0, 0}, {0, 0}};
    for (int i = 0; i < init1.n; i++) ex[0][init1.s[i] / 2] += init1.s[i] & 1 ? -1 : 1;
    for (int i = 0; i < init2.n; i++) ex[1][init2.s[i] / 2] += init2.s[i] & 1 ? -1 : 1;
    long det = ex[0][0] * ex[1][1] - ex[0][1] * ex[1][0];
    printf("input: r1=%s r2=%s  cap=%d  max-states=%llu  abelianization det=%ld%s\n",
           argv[1], argv[2], cap, (unsigned long long)max_states, det,
           det == 1 || det == -1 ? " (trivial H1, ok)" : " (NOT trivial group!)");

    uint64_t tsize = 1;
    while (tsize < max_states * 2) tsize <<= 1;
    table = calloc(tsize, sizeof(uint32_t));
    node_off = malloc(max_states * sizeof(uint64_t));
    node_parent = malloc(max_states * sizeof(uint32_t));
    node_move = malloc(max_states);
    arena_capacity = max_states * (uint64_t)(2 * cap + 4); /* worst case: no reallocs */
    arena = malloc(arena_capacity);
    if (!table || !node_off || !node_parent || !node_move || !arena) {
        fprintf(stderr, "allocation failed — lower --max-states\n");
        return 1;
    }
    table_mask = tsize - 1;

    uint8_t key[WBUF * 2];
    int klen = pack(&init1, &init2, key);
    int is_new;
    uint32_t root = lookup_or_insert(key, klen, UINT32_MAX, 0, &is_new);
    node_parent[root] = UINT32_MAX;
    bucket_push(init1.n + init2.n, root);
    if (target_len && target_len == klen && !memcmp(key, target_key, (size_t)klen)) {
        printf("TARGET REACHED at the start state (inputs are identical)\n");
        return 0;
    }
    if (target_len && (target_key[0] > cap || target_key[1] > cap)) {
        fprintf(stderr, "target exceeds cap — raise --cap\n");
        return 2;
    }

    uint64_t expanded = 0;
    int truncated = 0;
    int best = init1.n + init2.n;

    for (;;) {
        /* greedy: always expand from the lowest nonempty bucket */
        int tlen;
        for (tlen = 2; tlen <= 2 * cap; tlen++)
            if (buckets[tlen].head < buckets[tlen].n) break;
        if (tlen > 2 * cap) break;
        uint32_t cur = buckets[tlen].v[buckets[tlen].head++];
        Word r1, r2;
        unpack(arena + node_off[cur], &r1, &r2);
        if (r1.n + r2.n < best) {
            best = r1.n + r2.n;
            fprintf(stderr, "  new best total length %d after %llu expansions\n",
                    best, (unsigned long long)expanded);
        }
        expanded++;
        if ((expanded & 0xFFFFFF) == 0)
            fprintf(stderr, "  ... expanded %lluM states, stored %uM, frontier len %d\n",
                    (unsigned long long)(expanded >> 20), n_nodes >> 20, tlen);

        for (int m = 0; m < 14; m++) {
            Word n1 = r1, n2 = r2;
            if (!apply_move(m, &n1, &n2, cap)) continue;
            klen = pack(&n1, &n2, key);
            uint32_t child = lookup_or_insert(key, klen, cur, (uint8_t)m, &is_new);
            if (is_new == -1) { truncated = 1; goto done; }
            if (!is_new) continue;
            if (target_len && target_len == klen && !memcmp(key, target_key, (size_t)klen)) {
                printf("TARGET REACHED: AC-equivalent within cap %d "
                       "(expanded=%llu stored=%u)\n",
                       cap, (unsigned long long)expanded, n_nodes);
                return 0;
            }
            if (!target_len && n1.n + n2.n == 2) {
                printf("expanded=%llu stored=%u\n", (unsigned long long)expanded, n_nodes);
                replay_and_print(child, &init1, &init2, cap);
                return 0;
            }
            bucket_push(n1.n + n2.n, child);
        }
    }
done:
    printf("%s: %s. expanded=%llu stored=%u best_total_length=%d cap=%d\n",
           truncated ? "TRUNCATED (state budget hit — inconclusive)"
                     : "EXHAUSTED (rigorous: no path within cap)",
           target_len ? "target NOT reachable" : "no trivialization",
           (unsigned long long)expanded, n_nodes, best, cap);
    return truncated ? 3 : 4;
}
