/*
 * search.c — search for counterexamples to Euler's sum of powers conjecture.
 *
 * Primary target: (6,1,5)  a^6+b^6+c^6+d^6+e^6 = f^6  in positive integers.
 * Validation mode: (5,1,4) a^5+b^5+c^5+d^5 = f^5 (Lander–Parkin 1966:
 *   27^5+84^5+110^5+133^5 = 144^5), which exercises the same DFS machinery.
 *
 * Structure exploited for k=6 (primitive solutions; any solution is a
 * multiple of a primitive one with smaller f, so searching primitives up
 * to N is exhaustive up to N):
 *
 *   x^6 mod 7 is 0 (7|x) or 1;  same shape mod 9 (prime 3) and mod 8 (2).
 * With 5 terms the coprime-count on the left must equal f's (0 or 1) and
 * count 0 forces a common factor, so a primitive solution has EXACTLY ONE
 * term coprime to each of 2, 3, 7 (bits may share a term) and f coprime
 * to 42. Consequence used here: once an exemption bit is spent, every
 * deeper term must be divisible by that prime — the DFS strides by
 * lcm(spent primes), up to 42.
 *
 * Residue sieves: sums of j sixth powers hit few residues mod 13 and 43;
 * remainders are tracked incrementally (no 128-bit division in the hot path).
 *
 * Build:  gcc -O3 -march=native -fopenmp -o search search.c -lm
 * Run:    ./search <k:5|6> <fmin> <fmax>
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>

typedef unsigned __int128 u128;
typedef uint64_t u64;
typedef uint32_t u32;

static int K, NTERMS;
static u64 FMAX;
static u128 *P;               /* P[x] = x^K, x <= FMAX */

static u128 ipow(u64 x, int k){ u128 r=1, b=x; while(k){ if(k&1) r*=b; b*=b; k>>=1; } return r; }

/* ---- residue sieves ---- */
#define NM 2
static const u32 MODS[NM] = {13, 43};
static u64 sumres[NM][6];     /* residues achievable by a sum of j k-th powers */
static u64 fpres[NM];         /* residues of a single k-th power */
static u32 powmod[NM][64];    /* x^K mod m, indexed by x mod m */

static void init_sieves(void){
    for(int m=0;m<NM;m++){
        u32 md = MODS[m]; u64 single=0;
        for(u32 x=0;x<md;x++){ u64 r=1; for(int i=0;i<K;i++) r=r*x%md; powmod[m][x]=(u32)r; single|=1ULL<<r; }
        fpres[m]=single; sumres[m][0]=1ULL;
        for(int j=1;j<6;j++){
            u64 prev=sumres[m][j-1], cur=0;
            for(u32 r=0;r<md;r++) if(prev>>r&1)
                for(u32 s=0;s<md;s++) if(single>>s&1)
                    cur |= 1ULL<<((r+s)%md);
            sumres[m][j]=cur;
        }
    }
}

/* ---- integer k-th root ---- */
static u64 iroot(u128 R){
    if(R==0) return 0;
    u64 r = (u64)powl((long double)R, 1.0L/K);
    while(ipow(r+1,K) <= R) r++;
    while(r>0 && ipow(r,K) > R) r--;
    return r;
}

/* class(x): bit0 x odd, bit1 3∤x, bit2 7∤x (k=6 exemption bits) */
static inline int cls(u64 x){ return (int)((x&1) | ((x%3!=0)<<1) | ((x%7!=0)<<2)); }
static const int STRIDE[8] = {1,2,3,6,7,14,21,42};   /* lcm of spent primes */

static long long nfound = 0;
static void report(u64 f, u64 *t){
    printf("SOLUTION k=%d f=%llu :", K, (unsigned long long)f);
    for(int i=0;i<NTERMS;i++) printf(" %llu", (unsigned long long)t[i]);
    printf("\n"); fflush(stdout);
    #pragma omp atomic
    nfound++;
}

/* DFS: t[0..lvl-1] chosen (non-increasing), R remaining, rm[]=R mod MODS,
 * usedmask = exemption bits spent. */
static void dfs(u64 f, u64 *t, int lvl, u128 R, u32 rm0, u32 rm1, u64 maxv, int usedmask){
    int j = NTERMS - lvl;
    if(!(sumres[0][j]>>rm0 & 1)) return;
    if(!(sumres[1][j]>>rm1 & 1)) return;

    if(j==1){
        if(!(fpres[0]>>rm0&1) || !(fpres[1]>>rm1&1)) return;
        u64 e = iroot(R);
        if(e==0 || e>maxv || ipow(e,K)!=R) return;
        if(K==6){
            int c=cls(e);
            if((usedmask&c) || (usedmask|c)!=7) return;
        }
        t[lvl]=e; report(f,t);
        return;
    }
    u64 hi = iroot(R); if(hi>maxv) hi=maxv;
    u64 lo = iroot((R + j - 1)/j);
    if(ipow(lo,K)*(u128)j < R) lo++;

    int stride = (K==6) ? STRIDE[usedmask] : 1;
    u64 x = hi;
    if(stride>1) x -= x % stride;               /* x must be divisible by every spent prime */
    for(; x>=lo && x>0; x-=stride){
        int c=0;
        if(K==6){
            c = cls(x);
            if(c & usedmask) continue;          /* spent bit reintroduced (only when stride==1 partially) */
        }
        u128 xp = P[x];
        u32 m0 = powmod[0][x%13], m1 = powmod[1][x%43];
        t[lvl]=x;
        dfs(f, t, lvl+1, R-xp,
            rm0 + 13 - m0 >= 13 ? rm0 - m0 : rm0 + 13 - m0,
            rm1 + 43 - m1 >= 43 ? rm1 - m1 : rm1 + 43 - m1,
            x, usedmask|c);
        if(x < (u64)stride) break;
    }
}

int main(int argc, char **argv){
    if(argc<4){ fprintf(stderr,"usage: %s <k:5|6> <fmin> <fmax>\n", argv[0]); return 2; }
    K = atoi(argv[1]); NTERMS = K-1;
    u64 fmin = strtoull(argv[2],0,10), fmax = strtoull(argv[3],0,10);
    FMAX = fmax;
    init_sieves();
    P = malloc(sizeof(u128)*(FMAX+1));
    for(u64 x=0;x<=FMAX;x++) P[x]=ipow(x,K);

    if(K==6){  /* arithmetic self-test on the known (6,3,3) identity */
        if(ipow(3,6)+ipow(19,6)+ipow(22,6) != ipow(10,6)+ipow(15,6)+ipow(23,6)){
            fprintf(stderr,"SELF-TEST FAILED\n"); return 1; }
    }

    #pragma omp parallel for schedule(dynamic,1)
    for(u64 f=fmin; f<=fmax; f++){
        if(K==6 && (f%2==0 || f%3==0 || f%7==0)) continue;
        u128 T = P[f];
        u32 r0=(u32)(T%13), r1=(u32)(T%43);
        u64 t[8];
        dfs(f, t, 0, T, r0, r1, f-1, 0);
    }
    fprintf(stderr,"done: k=%d f in [%llu,%llu], solutions found: %lld\n",
            K,(unsigned long long)fmin,(unsigned long long)fmax,nfound);
    return 0;
}
