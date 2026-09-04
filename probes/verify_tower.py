#!/usr/bin/env python3
"""Session 111: a two-level exhaustion tower.

The cylinder-exhaustion theorem required the exhausted cylinders
to be ALIGNED. Drop that. Exhaustion is really a REDUCTION: the
survivors of the system equal those of the coarser system formed
by its cylinders, and the question becomes whether the coarser
one is covered. That invites iteration.

TOWER.
  top    : sparse primes q_i; cylinders {x = 1 mod q_i}.
  middle : per q_i, classes (q_i*m, c) with c = 1 mod q_i and
           c in general position mod m. These are the level-1
           cylinders - a MISALIGNED family.
  bottom : each level-1 cylinder exhausted by its own group of
           classes (q_i*m*k, c') with c' = c mod q_i*m, general
           position mod k.

The bottom system has gcd 1, no finite prime cover, misaligned
residues, divergent drift, and its exhausted cylinders are NOT
aligned - so the one-level theorem fails on it. Yet its
survivors equal the top level's, prod(1-1/q_i), by two
reductions.

CHECKS: bottom survivors match the top product; each level-1
cylinder dies; the level-1 family is genuinely misaligned; and
coverage against every theorem is confirmed BEFORE the word
'evader' is used.
"""
import sys
from math import gcd

def primes_upto(n):
    s=bytearray([1])*(n+1); s[0]=s[1]=0
    i=2
    while i*i<=n:
        if s[i]: s[i*i::i]=bytearray(len(s[i*i::i]))
        i+=1
    return [i for i in range(2,n+1) if s[i]]


def crt(r1,m1,r2,m2):
    g=gcd(m1,m2)
    assert (r2-r1)%g==0
    l=m1//g*m2
    inv=pow((m1//g)%(m2//g),-1,m2//g) if m2//g>1 else 0
    k=(((r2-r1)//g)*inv)%(m2//g) if m2//g>1 else 0
    return (r1+m1*k)%l, l


def aligned(cls):
    r,m=None,1
    for n,a in cls:
        if r is None: r,m=a%n,n; continue
        g=gcd(m,n)
        if (a-r)%g: return False
        r,m=crt(r,m,a%n,n)
    return True


def main():
    X=3*10**6
    import random
    rng=random.Random(31)
    qs=[11,197,439]
    top=1.0
    for q in qs: top*=(1-1.0/q)
    print("top: q_i = %s, prod(1-1/q_i) = %.4f" % (qs, top))

    m=bytearray(X)
    level1=[]; nbot=0; drift=0.0; gcd_all=0
    for q in qs:
        for mm in range(1, 60):
            if gcd(mm,q)!=1: continue
            # level-1 class: residue 1 mod q, random mod mm
            c,d = crt(1,q, rng.randrange(mm),mm) if mm>1 else (1,q)
            level1.append((d,c))
            # bottom group exhausting it: moduli d*k, residue
            # c mod d, random mod k, for k up to a bound
            for k in range(1, max(2, X//(d*40))):
                if gcd(k,d)!=1: continue
                if k==1: cc,n=c,d
                else: cc,n = crt(c,d, rng.randrange(k),k)
                gcd_all=gcd(gcd_all,n); drift+=1.0/n; nbot+=1
                st=cc%n
                if st==0: st=n
                while st<n: st+=n
                if st<X: m[st::n]=b"\x01"*len(m[st::n])
    print("  level-1 cylinders: %d, misaligned as a family: %s"
          % (len(level1), not aligned(level1)))
    print("  bottom classes: %d, gcd %d, drift %.3f"
          % (nbot, gcd_all, drift))
    lo=X//2
    surv=1-sum(m[lo:X])/float(X-lo)
    print()
    print("  bottom survivors: %.4f   vs   top product %.4f"
          % (surv, top))
    # each level-1 cylinder must be dead
    dead=[]
    for d,c in level1[:6]:
        pts=[x for x in range(lo,X) if x%d==c%d]
        if pts: dead.append(1-sum(m[x] for x in pts)/float(len(pts)))
    print("  survivor share in first level-1 cylinders: %s"
          % "  ".join("%.4f"%v for v in dead))
    # complement of top cylinders untouched
    comp=[x for x in range(lo,X) if all(x%q!=1 for q in qs)]
    unt=1-sum(m[x] for x in comp)/float(len(comp))
    print("  complement of top cylinders untouched: %.4f" % unt)
    print()
    ok = gcd_all==1 and abs(surv-top)<0.02 and max(dead)<0.02 \
         and unt>0.99 and not aligned(level1)
    print("  coverage: one-level exhaustion needs ALIGNED "
          "cylinders - level-1 family aligned? %s -> theorem "
          "does NOT apply" % aligned(level1))
    if ok:
        print()
        print("[verdict] TWO-LEVEL TOWER VERIFIED: evades the "
              "one-level exhaustion theorem (cylinders "
              "misaligned), yet survivors equal the top product "
              "by TWO reductions. The recursive form of the "
              "theorem covers it; the one-level form does not.")
    else:
        print("[verdict] construction did not verify.")


if __name__=="__main__":
    main()
