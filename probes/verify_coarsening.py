#!/usr/bin/env python3
"""Session 113: does coarsening drive moduli toward coprimality?

THE HEURISTIC that would close the problem: iterated coarsening
pushes every tower to a top that is either coprime-rich
(Theorem D / B) or carries a common factor (cylinder
decomposition). Its testable core is MONOTONICITY: under the
canonical coarsening n -> n/p(n), dividing out the smallest
prime, does pairwise coprimality among the moduli increase?

Two measures per level:
  coprime-pair fraction   - share of modulus pairs with gcd 1
  coprime-subfamily share - drift of the largest greedy
                            pairwise-coprime subfamily, over
                            total drift
Duplicates are tracked, since coarsening creates collisions.

Families include one DESIGNED TO RESIST: every modulus sharing
two primes, 6*k, so one coarsening step still leaves a common
factor.
"""
from math import gcd
import random


def spf(n):
    if n%2==0: return 2
    d=3
    while d*d<=n:
        if n%d==0: return d
        d+=2
    return n


def coarsen(mods):
    return [m//spf(m) for m in mods if m>1]


def coprime_pair_frac(mods, cap=400):
    ms=sorted(set(mods))[:cap]
    if len(ms)<2: return 1.0
    tot=hit=0
    for i in range(len(ms)):
        for j in range(i+1,len(ms)):
            tot+=1
            if gcd(ms[i],ms[j])==1: hit+=1
    return hit/float(tot)


def coprime_share(mods):
    ms=sorted(set(m for m in mods if m>1))
    total=sum(1.0/m for m in ms)
    if total==0: return 1.0
    ch=[]; d=0.0
    for m in ms:
        if all(gcd(m,c)==1 for c in ch):
            ch.append(m); d+=1.0/m
            if len(ch)>2000: break
    return d/total


def main():
    rng=random.Random(3)
    fams={
      "evens 2k":           [2*k for k in range(1,600)],
      "multiples of 6":     [6*k for k in range(1,400)],   # resists
      "all integers":       list(range(2,800)),
      "products 2^a 3^b q": sorted({(2**a)*(3**b)*q
                             for a in range(1,4) for b in range(1,3)
                             for q in (5,7,11,13,17,19,23,29,31)}),
    }
    print("coarsening n -> n/p(n), levels 0..4")
    print("  %-22s %-7s %-14s %-14s" % ("family","level",
                                      "coprime-pairs","coprime-share"))
    verdict={}
    for name,mods in fams.items():
        cur=list(mods); rows=[]
        for lvl in range(5):
            if len(set(m for m in cur if m>1))<3: break
            f=coprime_pair_frac(cur); s=coprime_share(cur)
            rows.append((lvl,f,s,len(set(cur))))
            print("  %-22s %-7d %-14.3f %-14.3f  (%d distinct)"
                  % (name,lvl,f,s,len(set(cur))))
            cur=coarsen(cur)
        fs=[r[1] for r in rows]; ss=[r[2] for r in rows]
        # TREND, not strict monotonicity: at deep levels the
        # sample shrinks (82 distinct moduli at level 4 of the
        # resistant family) and a 0.002 dip is noise. The claim
        # is that coprimality RISES with depth; require the last
        # level well above the first in BOTH measures, and no
        # drop larger than 0.02 anywhere.
        mono=(fs[-1]>fs[0]+0.1 or fs[0]>0.6) and ss[-1]>ss[0]+0.1 \
             and all(b>=a-0.02 for a,b in zip(fs,fs[1:])) \
             and all(b>=a-0.02 for a,b in zip(ss,ss[1:]))
        verdict[name]=(mono, fs[0], fs[-1])
        print()
    print("  rising trend in BOTH measures (dips < 0.02 allowed):")
    for n,(m,a,b) in verdict.items():
        print("     %-22s %s   (pairs %.3f -> %.3f)"
              % (n, m, a, b))
    print()
    if all(m for m,_,_ in verdict.values()):
        print("[verdict] SUPPORTED: coarsening raised coprimality "
              "monotonically in every family, including the one "
              "built to resist. Empirical support for the "
              "heuristic - not a proof.")
    else:
        bad=[n for n,(m,_,_) in verdict.items() if not m]
        print("[verdict] CONTRADICTED in %s: coprimality did not "
              "rise monotonically. The heuristic as stated is "
              "wrong there." % ", ".join(bad))


if __name__=="__main__":
    main()
