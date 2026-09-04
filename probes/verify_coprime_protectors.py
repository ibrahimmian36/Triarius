#!/usr/bin/env python3
"""Session 118: the coprime-protector theorem.

Groups G_j on cylinders C_j = c_j + p_j Z with p_j pairwise
coprime; group j is a coprime layer (p_j q, a), q prime, all
layer primes distinct across groups and coprime to every p.
Then membership and kill are independent across j, so
   nu = prod_j (1 - kappa_j / p_j)   EXACTLY,
kappa_j = 1 - prod_{q in Q_j}(1 - 1/q) >= 1 - e^{-S_j}.
Positive nu forces sum kappa_j/p_j < inf, hence the split
   S_j >= 1 : sum 1/p_j < inf        (summable protectors)
   S_j <  1 : sum S_j/p_j < inf      (convergent drift).

INSTANCE: protectors = odd primes to 300; big groups (S = 0.40)
on 11 and 23, small groups (S ~ 1/(40 log^2 p)) elsewhere, layer
primes distinct across groups. Distinct primes make any layer
with S >= 1 need primes to 1e9, so the instance checks the
product formula and the forced split, not divergence; the
mutation handles divergence analytically.
MUTATION: big groups (S = 1.5) on EVERY protector: the product
prod(1 - kappa/p) over the first N protectors is reported
exactly; the infinite product is 0 by the theorem.
"""
import math, random, sys


def primes_upto(n):
    s=bytearray([1])*(n+1); s[0]=s[1]=0
    for i in range(2,int(n**.5)+1):
        if s[i]: s[i*i::i]=bytearray(len(s[i*i::i]))
    return [i for i in range(n+1) if s[i]]


def instance(seed=5):
    """protectors: odd primes to 300; big groups on 11 and 23 with
    S = 0.40; small groups S = 1/(40 log^2 p) elsewhere; layer
    primes distinct, drawn in order from (300, 1e7)."""
    rng=random.Random(seed)
    prot=[p for p in primes_upto(300) if p>=3]
    bigset={11,23}
    pool=[q for q in primes_upto(10**7) if q>300]
    groups={}; used=0
    for p in prot:
        S=0.40 if p in bigset else 1.0/(40*math.log(p)**2)
        qs=[]; s=0.0
        while s<S:
            q=pool[used]; used+=1; qs.append(q); s+=1.0/q
            if used>=len(pool): sys.exit("pool exhausted")
        c=rng.randrange(p)
        groups[p]=(c,[(q,rng.randrange(q)) for q in qs])
    assert used==len({q for p in groups for q,_ in groups[p][1]})   # distinct across groups
    return prot,bigset,groups


def kappa(qs):
    pr=1.0
    for q,_ in qs: pr*=(1-1.0/q)
    return 1-pr


def main():
    prot,bigset,groups=instance()
    ncls=sum(len(g[1]) for g in groups.values())
    nu=1.0; drift=0.0; sbig=0.0; ssmall=0.0
    for p,(c,qs) in groups.items():
        k=kappa(qs); nu*=(1-k/p); S=sum(1.0/q for q,_ in qs); drift+=S/p
        if p in bigset: sbig+=1.0/p
        else: ssmall+=S/p
    print("instance: %d protectors (%d big), %d classes" % (len(prot),len(bigset),ncls))
    print("[1] nu = prod(1 - kappa_j/p_j) = %.5f ; drift of truncation %.4f" % (nu,drift))
    print("    forced split: sum_{big} 1/p_j = %.4f (sparse, summable); sum_{small} S_j/p_j = %.4f" % (sbig,ssmall))
    rng=random.Random(1); n=6000; surv=0
    allq={q for p in groups for q,_ in groups[p][1]}
    for _ in range(n):
        r={p:rng.randrange(p) for p in prot}; rq={q:rng.randrange(q) for q in allq}
        ok=True
        for p,(c,qs) in groups.items():
            if r[p]==c and any(rq[q]==a for q,a in qs): ok=False; break
        surv+=ok
    mc=surv/float(n); ok1=abs(mc-nu)<3*math.sqrt(nu*(1-nu)/n)+0.003
    print("    MC on Zhat (%d samples): %.4f   agree: %s" % (n,mc,ok1))
    print("[2] mutation: S = 1.5 on EVERY protector; exact partial products of prod(1 - kappa/p):")
    kap=1-math.exp(-1.5)
    prods=[]; pr=1.0; P=primes_upto(2*10**6)
    marks=(100,10**3,10**4,10**5,10**6,2*10**6); mi=0
    for p in P:                      # product over primes p <= N, as in the paper
        if p<3: continue
        while mi<len(marks) and p>marks[mi]:
            prods.append((marks[mi],pr)); mi+=1
        if mi==len(marks): break
        pr*=(1-kap/p)
    for N,v in prods: print("      N=%8d  product %.4f" % (N,v))
    falling=all(b<a for (_,a),(_,b) in zip(prods,prods[1:]))
    print("    falling toward 0 (like (log N)^-%.2f): %s" % (kap,falling))
    print()
    print("[verdict]", "PASS" if ok1 and falling else "FAIL")
    return {"nprot":len(prot),"nbig":len(bigset),"ncls":ncls,"nu":round(nu,5),"drift":round(drift,4),
            "sbig":round(sbig,4),"ssmall":round(ssmall,4),"mc":round(mc,4),"prods":[(N,round(v,4)) for N,v in prods],"ok":ok1 and falling}


if __name__=="__main__":
    r=main(); sys.exit(0 if r["ok"] else 1)
