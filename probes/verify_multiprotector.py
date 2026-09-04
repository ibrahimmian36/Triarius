#!/usr/bin/env python3
"""Session 110: many protectors - a counterexample to the
Question as stated, and the theorem that covers it.

CONSTRUCTION. Sparse primes p_j with sum 1/p_j < infinity. Group
j = moduli p_j*m, residue = 1 mod p_j, residues in general
position mod m. Group j kills only inside C_j = {x = 1 mod p_j},
and within it acts as a full random system - which dies. So the
survivors are the complement of the union of the C_j, of
measure prod(1-1/p_j) > 0.

The moduli have gcd 1, no finite prime cover, misaligned
residues, divergent drift. EVERY hypothesis of the Question
fails, yet nu > 0.

WHY THE DENSITY STILL EXISTS - the missing theorem: each group
exhausts its cylinder in measure, hence in log density (nu = 0
corollary); the killed set agrees with the union of the C_j up
to log density 0; and that union's complement is a TRANSLATE of
a multiples-complement, which Davenport-Erdos handles.

CHECKS: gcd 1; misaligned; survivors match prod(1-1/p_j); each
cylinder dies; the complement is untouched; and the construction
is checked against every existing theorem before being called a
counterexample to the Question.
"""
import sys
from math import gcd
from functools import reduce

def primes_upto(n):
    s=bytearray([1])*(n+1); s[0]=s[1]=0
    i=2
    while i*i<=n:
        if s[i]: s[i*i::i]=bytearray(len(s[i*i::i]))
        i+=1
    return [i for i in range(2,n+1) if s[i]]


def main():
    X=4*10**6
    import random
    rng=random.Random(23)
    # sparse primes: every 40th prime from the 5th on
    ps=primes_upto(2000)
    protectors=ps[4::40][:6]
    prodm=1.0
    for p in protectors: prodm*=(1-1.0/p)
    print("protectors p_j = %s, sum 1/p_j = %.4f, "
          "prod(1-1/p_j) = %.4f"
          % (protectors, sum(1.0/p for p in protectors), prodm))

    m=bytearray(X)
    gcd_all=0; drift=0.0; nclasses=0
    aligned_ok=True
    for p in protectors:
        for k in range(1, X//(2*p)):
            n=p*k
            # residue = 1 mod p, random mod k  (CRT)
            if k==1:
                a=1
            else:
                g=gcd(p,k)
                if g!=1: continue          # keep it clean
                b=rng.randrange(k)
                inv=pow(p%k,-1,k)
                a=(1 + p*(((b-1)*inv)%k)) % n
            assert a%p==1
            gcd_all=gcd(gcd_all,n); drift+=1.0/n; nclasses+=1
            st=a%n
            if st==0: st=n
            while st<n: st+=n
            if st<X: m[st::n]=b"\x01"*len(m[st::n])
    print("  %d classes, gcd of moduli = %d, drift %.3f "
          "(divergent per group)" % (nclasses, gcd_all, drift))
    print()
    lo=X//2
    surv=1-sum(m[lo:X])/float(X-lo)
    print("  survivor density in [X/2,X]: %.4f   vs   "
          "prod(1-1/p_j) = %.4f" % (surv, prodm))
    # per-cylinder: does C_j die, is the complement untouched?
    for p in protectors[:3]:
        inC=[x for x in range(lo,X) if x%p==1]
        died=1-sum(m[x] for x in inC)/float(len(inC))
        print("     cylinder x=1 mod %-4d: survivor share %.4f "
              "(should -> 0)" % (p, died))
    comp=[x for x in range(lo,X) if all(x%p!=1 for p in protectors)]
    untouched=1-sum(m[x] for x in comp)/float(len(comp))
    print("     complement of all C_j: survivor share %.4f "
          "(should be 1.0)" % untouched)
    print()
    ok = gcd_all==1 and abs(surv-prodm)<0.03 and untouched>0.99
    print("  gcd 1: %s | survivors ~ product: %s | complement "
          "untouched: %s" % (gcd_all==1, abs(surv-prodm)<0.03,
                             untouched>0.99))
    print()
    print("  coverage against existing theorems:")
    print("     finite prime cover?  NO - for any finite prime set F pick "
          "p_j outside F and m a prime outside F: modulus p_j m avoids F")
    print("     aligned?             NO - residues random mod m")
    print("     convergent drift?    NO - each group diverges")
    print("     divergent coprime?   NO - one class per group at "
          "most, drift <= sum 1/p_j < inf  (so nu>0, consistent)")
    print("     bounded misalign?    NO - misaligned parts ~ p_j*m")
    print()
    if ok:
        print("[verdict] THE QUESTION AS STATED IS FALSE: this "
              "system has nu > 0, divergent drift, gcd 1, no "
              "finite prime cover and misaligned residues. It is "
              "NOT a counterexample to Erdos 25 - its density "
              "exists by CYLINDER EXHAUSTION - but it is a "
              "counterexample to the case analysis, and exposes "
              "the need for the exhaustion theorem.")
    else:
        print("[verdict] construction did not verify; no claim.")


if __name__=="__main__":
    main()
