/* audit.c — completeness check of the (f,t) candidate enumeration in caseA2.
 * For sample f, count t in [1,f), t coprime to 42, with t^6 ≡ f^6 (mod 42^6)
 * two ways: (1) brute force over all t; (2) via the 144 roots of unity.
 * Any mismatch means the fast enumeration is wrong. */
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
typedef unsigned __int128 u128;
typedef uint64_t u64;
static const u64 M = 5489031744ULL;
static u64 pow6m(u64 x,u64 m){ u64 r=x%m; u128 s=(u128)r*r%m; s=s*s%m; return (u64)((u128)s*((u128)r*r%m)%m); }
static u64 roots[160]; static int nroots=0;
static void build_roots(void){
    u64 r2[8]; int n2=0; for(u64 x=1;x<64;x+=2) if(pow6m(x,64)==1) r2[n2++]=x;
    u64 r3[8]; int n3=0; for(u64 x=1;x<729;x++) if(x%3&&pow6m(x,729)==1) r3[n3++]=x;
    u64 r7[8]; int n7=0; for(u64 x=1;x<117649;x++) if(x%7&&pow6m(x,117649)==1) r7[n7++]=x;
    for(int i=0;i<n2;i++)for(int j=0;j<n3;j++)for(int k=0;k<n7;k++){
        u64 m12=64*729ULL,x12=0,x=0;
        for(u64 t=r2[i];;t+=64) if(t%729==r3[j]){x12=t;break;}
        for(u64 t=x12;;t+=m12) if(t%117649==r7[k]){x=t;break;}
        roots[nroots++]=x;
    }
}
int main(int argc,char**argv){
    build_roots();
    printf("nroots=%d\n",nroots);
    u64 fs[]= {730001, 800003, 899999, 999983, 1299983, 54321011};
    for(int i=0;i<6;i++){
        u64 f=fs[i];
        if(f%2==0||f%3==0||f%7==0){ printf("f=%llu not coprime to 42, skip\n",(unsigned long long)f); continue; }
        u64 f6=pow6m(f,M);
        long brute=0;
        for(u64 t=1;t<f;t++){
            if(t%2==0||t%3==0||t%7==0) continue;
            if(pow6m(t,M)==f6) brute++;
        }
        long fast=0;
        for(int j=0;j<nroots;j++){
            u64 t=(u64)((u128)roots[j]*f%M);
            if(t>0&&t<f) fast++;
        }
        printf("f=%llu brute=%ld fast=%ld %s\n",(unsigned long long)f,brute,fast,brute==fast?"OK":"MISMATCH");
    }
    return 0;
}
