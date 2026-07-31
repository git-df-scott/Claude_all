/*
 * caseA.c — deep search of the "concentrated exemption" case of
 * a^6+b^6+c^6+d^6+e^6 = f^6.
 *
 * In a primitive solution exactly one term is odd, one is coprime to 3,
 * one is coprime to 7 (see search.c). Case A: all three exemptions sit on
 * the SAME term t, so the other four terms are divisible by 42 and
 *
 *      f^6 ≡ t^6 (mod 42^6),   f, t coprime to 42.
 *
 * The units t/f are then sixth roots of unity mod M = 42^6 = 5489031744;
 * there are exactly 4*6*6 = 144 of them (CRT over 2^6, 3^6, 7^6). For each
 * f we enumerate t = u*f mod M over the 144 roots u and keep t < f: the
 * expected number of candidate pairs below F is ~ 144 F^2 / (2M) * (2/7),
 * negligible until f is in the millions — so this case can be swept to
 * heights generic enumeration cannot reach.
 *
 * A surviving pair (f,t) needs m = (f^6 - t^6)/42^6 to be a sum of four
 * sixth powers (the quotients of the terms by 42). Key second-stage
 * structure: sixth powers are 0/1 mod 8, 9, 7, so the number of odd parts
 * is exactly m mod 8, the number of parts coprime to 3 is exactly m mod 9,
 * and coprime to 7 exactly m mod 7 (each must be <= 4, else infeasible).
 * The DFS tracks these budgets: when a budget hits 0 the remaining parts
 * are forced divisible by that prime (stride up to 42); when it equals the
 * number of remaining parts, ALL of them are forced coprime. Residue masks
 * mod 64, 27, 49, 13, 43 (tracked incrementally) prune further.
 *
 * Build:  gcc -O3 -march=native -fopenmp -o caseA caseA.c -lm
 * Run:    ./caseA <fmin> <fmax>        (fmax <= 1e8: u128 limit)
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>

typedef unsigned __int128 u128;
typedef uint64_t u64;
typedef uint32_t u32;

static const u64 M = 5489031744ULL;            /* 42^6 */
static u64 roots[160]; static int nroots = 0;

static u128 ipow6(u64 x){ u128 s=(u128)x*x; return s*s*s; }
static u64 pow6m(u64 x, u64 m){ u64 r=x%m; u128 s=(u128)r*r%m; s=s*s%m; return (u64)((u128)s*((u128)r*r%m)%m); }

static void build_roots(void){
    u64 r2[8]; int n2=0; for(u64 x=1;x<64;x+=2) if(pow6m(x,64)==1) r2[n2++]=x;
    u64 r3[8]; int n3=0; for(u64 x=1;x<729;x++) if(x%3 && pow6m(x,729)==1) r3[n3++]=x;
    u64 r7[8]; int n7=0; for(u64 x=1;x<117649;x++) if(x%7 && pow6m(x,117649)==1) r7[n7++]=x;
    for(int i=0;i<n2;i++)for(int j=0;j<n3;j++)for(int k=0;k<n7;k++){
        u64 m12=64*729ULL, x12=0, x=0;
        for(u64 t=r2[i];;t+=64) if(t%729==r3[j]){ x12=t; break; }
        for(u64 t=x12;;t+=m12) if(t%117649==r7[k]){ x=t; break; }
        roots[nroots++]=x;
    }
}

/* residue masks mod 27, 49, 13, 43, 64: sums of j sixth powers */
#define NM 4
static const u32 MODS[NM] = {27,49,13,43};
static u64 sumres[NM][5];
static u64 sumres64[5];
static u32 powmod6[NM][64];      /* x^6 mod MODS[m], indexed by x mod MODS[m] */
static void init_sieves(void){
    for(int m=0;m<NM;m++){
        u32 md=MODS[m]; u64 single=0;
        for(u32 x=0;x<md;x++){ u64 r=pow6m(x,md); powmod6[m][x]=(u32)r; single|=1ULL<<r; }
        sumres[m][0]=1;
        for(int j=1;j<5;j++){
            u64 prev=sumres[m][j-1], cur=0;
            for(u32 r=0;r<md;r++) if(prev>>r&1)
                for(u32 s=0;s<md;s++) if(single>>s&1)
                    cur|=1ULL<<((r+s)%md);
            sumres[m][j]=cur;
        }
    }
    u64 single=0;
    for(u32 x=0;x<64;x++) single|=1ULL<<(u32)pow6m(x,64);
    sumres64[0]=1;
    for(int j=1;j<5;j++){
        u64 prev=sumres64[j-1], cur=0;
        for(u32 r=0;r<64;r++) if(prev>>r&1)
            for(u32 s=0;s<64;s++) if(single>>s&1)
                cur|=1ULL<<((r+s)&63);
        sumres64[j]=cur;
    }
}

static u64 iroot6(u128 R){
    if(!R) return 0;
    u64 r=(u64)powl((long double)R,1.0L/6);
    while(ipow6(r+1)<=R) r++;
    while(r && ipow6(r)>R) r--;
    return r;
}

static long long nfound=0;
static void report(u64 f,u64 t,u64 *b,int nb){
    printf("SOLUTION f=%llu t=%llu parts42:",(unsigned long long)f,(unsigned long long)t);
    for(int i=0;i<nb;i++) printf(" %llu",(unsigned long long)(b[i]*42));
    printf("\n"); fflush(stdout);
    #pragma omp atomic
    nfound++;
}

/* DFS: R as a sum of j sixth powers, bases in [1,maxv], with exact budgets:
 * o2/o3/o7 = number of remaining bases that must be odd / coprime-3 / coprime-7.
 * r27,r49,r13,r43 = R mod 27,49,13,43 (incremental). */
static int dfs(u128 R,int j,u64 maxv,int o2,int o3,int o7,
               u32 r27,u32 r49,u32 r13,u32 r43,u64 *out){
    if(o2>j||o3>j||o7>j) return 0;
    if(!(sumres64[j]>>((u32)R&63)&1)) return 0;
    if(!(sumres[0][j]>>r27&1)) return 0;
    if(!(sumres[1][j]>>r49&1)) return 0;
    if(!(sumres[2][j]>>r13&1)) return 0;
    if(!(sumres[3][j]>>r43&1)) return 0;

    if(j==1){
        u64 e=iroot6(R);
        if(e==0||e>maxv||ipow6(e)!=R) return 0;
        if((int)(e&1)!=o2 || (int)(e%3!=0)!=o3 || (int)(e%7!=0)!=o7) return 0;
        out[0]=e; return 1;
    }
    u64 hi=iroot6(R); if(hi>maxv) hi=maxv;
    u64 lo=iroot6((R+j-1)/j);
    if(ipow6(lo)*(u128)j<R) lo++;
    /* forced divisibility: budget 0 => remaining bases divisible by that prime */
    int d2 = (o2==0)?2:1, d3 = (o3==0)?3:1, d7 = (o7==0)?7:1;
    int stride = d2*d3*d7;
    u64 x=hi; if(stride>1) x-=x%stride;
    for(; x>=lo && x>0; x-=(u64)stride){
        int x2=(int)(x&1), x3=(x%3!=0), x7=(x%7!=0);
        /* forced coprimality: budget == j => every remaining base coprime */
        if(o2==j && !x2) continue;   if(o2==0 && x2) continue;
        if(o3==j && !x3) continue;   if(o3==0 && x3) continue;
        if(o7==j && !x7) continue;   if(o7==0 && x7) continue;
        out[0]=x;
        u32 p27=powmod6[0][x%27], p49=powmod6[1][x%49], p13=powmod6[2][x%13], p43=powmod6[3][x%43];
        if(dfs(R-ipow6(x),j-1,x,o2-x2,o3-x3,o7-x7,
               r27+((r27<p27)?27:0)-p27, r49+((r49<p49)?49:0)-p49,
               r13+((r13<p13)?13:0)-p13, r43+((r43<p43)?43:0)-p43, out+1)) return 1;
        if(x<(u64)stride) break;
    }
    return 0;
}

/* try to write m as a sum of nb sixth powers (bases <= maxb); budgets from m */
static int try_decompose(u128 m,int nb,u64 maxb,u64 *out){
    int k2=(int)(m%8), k3=(int)(m%9), k7=(int)(m%7);
    if(k2>nb||k3>nb||k7>nb) return 0;
    return dfs(m,nb,maxb,k2,k3,k7,(u32)(m%27),(u32)(m%49),(u32)(m%13),(u32)(m%43),out);
}

static long long ncand=0, nsurv=0;

int main(int argc,char **argv){
    if(argc<3){ fprintf(stderr,"usage: %s <fmin> <fmax>\n",argv[0]); return 2; }
    u64 fmin=strtoull(argv[1],0,10), fmax=strtoull(argv[2],0,10);
    if(fmax>100000000ULL){ fprintf(stderr,"fmax capped at 1e8 (u128 overflow)\n"); return 2; }
    build_roots(); init_sieves();
    if(nroots!=144){ fprintf(stderr,"root count %d != 144\n",nroots); return 1; }
    for(int i=0;i<nroots;i++) if(pow6m(roots[i],M)!=1){ fprintf(stderr,"bad root\n"); return 1; }
    {   /* self-tests: recover known decompositions */
        u64 out[4];
        if(!try_decompose(ipow6(5)+ipow6(9)+ipow6(11)+ipow6(14),4,100,out))
            { fprintf(stderr,"SELF-TEST 1 FAILED\n"); return 1; }
        if(!try_decompose(ipow6(42)+ipow6(84)+ipow6(126)+ipow6(168),4,200,out))
            { fprintf(stderr,"SELF-TEST 2 FAILED\n"); return 1; }
        if(try_decompose((u128)12345677,4,100,out))
            { fprintf(stderr,"SELF-TEST 3 FAILED (false positive)\n"); return 1; }
    }

    #pragma omp parallel for schedule(dynamic,4096) reduction(+:ncand,nsurv)
    for(u64 f=fmin;f<=fmax;f++){
        if(f%2==0||f%3==0||f%7==0) continue;
        for(int i=0;i<nroots;i++){
            u64 t=(u64)((u128)roots[i]*f % M);
            if(t==0||t>=f) continue;
            ncand++;
            /* m = (f^6-t^6)/M via the four factors, dividing M out by gcds */
            u64 g[2]={f-t,f+t};
            u128 qq[2]={(u128)f*f-(u128)f*t+(u128)t*t,(u128)f*f+(u128)f*t+(u128)t*t};
            u64 rem=M; u128 m=1;
            for(int s=0;s<2;s++){
                u64 a=g[s],b=rem; while(b){u64 w=a%b;a=b;b=w;}
                rem/=a; m*=g[s]/a;
            }
            for(int s=0;s<2;s++){
                u64 a=(u64)(qq[s]%rem);
                if(a){ u64 x=a,y=rem; while(y){u64 w=x%y;x=y;y=w;} a=x; } else a=rem;
                rem/=a; m*=qq[s]/a;
            }
            if(rem!=1){ fprintf(stderr,"WARN rem=%llu f=%llu\n",(unsigned long long)rem,(unsigned long long)f); continue; }
            nsurv++;
            u64 out[4], maxb=(f-1)/42;
            if(try_decompose(m,4,maxb,out)) report(f,t,out,4);
            if(try_decompose(m,3,maxb,out)) report(f,t,out,3);
            if(try_decompose(m,2,maxb,out)) report(f,t,out,2);
            if(try_decompose(m,1,maxb,out)) report(f,t,out,1);
        }
    }
    fprintf(stderr,"done: f in [%llu,%llu] candidates=%lld processed=%lld found=%lld\n",
        (unsigned long long)fmin,(unsigned long long)fmax,ncand,nsurv,nfound);
    return 0;
}
