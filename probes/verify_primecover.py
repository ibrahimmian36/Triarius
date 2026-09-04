#!/usr/bin/env python3
"""Session 108: does divergent drift + nu > 0 force a finite
prime cover of the moduli?

WHY IT MIGHT. nu > 0 rules out any divergent coprime subfamily
(Theorem D), so the moduli must share primes heavily. If every
modulus is divisible by one of finitely many primes, a coprime
subfamily has at most that many members and bounded drift -
exactly the constraint. And an infinite family of two-element
prime supports that pairwise intersect must be a STAR, hence has
a common prime.

CONJECTURE: divergent drift with nu > 0 forces the moduli to be
covered, up to a convergent remainder, by finitely many primes.
Then the system is a multiples-of-d family plus a convergent
tail and the mixed corollary applies - and the open case is
empty.

TEST: build divergent-drift systems with UNBOUNDED prime support
and see whether any sustains survivors.

CONTROLS FIRST: aligned evens must SURVIVE (nu = 1/2), all-moduli
random must DIE - else the measurement discriminates nothing.
No property of an infinite system is inferred from a truncation.
"""
import sys

def primes_upto(n):
    s=bytearray([1])*(n+1); s[0]=s[1]=0
    i=2
    while i*i<=n:
        if s[i]: s[i*i::i]=bytearray(len(s[i*i::i]))
        i+=1
    return [i for i in range(2,n+1) if s[i]]


def survivor_trend(gen, X):
    m=bytearray(X)
    for n,a in gen:
        s=a%n
        if s==0: s=n
        while s<n: s+=n
        if s<X: m[s::n]=b"\x01"*len(m[s::n])
    out=[]; Y=16384
    while Y<=X:
        out.append(1-sum(m[Y//2:Y])/float(Y-Y//2)); Y*=2
    return out


def show(name, seq, cover):
    trend = "falling" if all(b<=a+1e-9 for a,b in zip(seq,seq[1:])) else "mixed"
    print("   %-30s %s   [%s]  prime cover: %s"
          % (name, "  ".join("%.4f"%v for v in seq[-5:]),
             trend, cover))


def main():
    X=4*10**6
    import random
    rng=random.Random(9)
    print("CONTROLS")
    ev=survivor_trend(((2*k,0) for k in range(1,X//2)), X)
    show("aligned evens (nu=1/2)", ev, "{2}")
    allr=survivor_trend(((n,rng.randrange(n))
                         for n in range(2,X//2)), X)
    show("all moduli, random (nu=0)", allr, "none")
    ok = ev[-1] > 0.4 and allr[-1] < 0.01
    print("   controls discriminate: %s" % ok)
    if not ok:
        print("   [ABORT] controls failed"); sys.exit(1)
    print()
    print("CANDIDATES: divergent drift, UNBOUNDED prime support")
    ps=primes_upto(int((X//2)**0.5)+1)
    pset=set(ps)

    def semiprimes():
        for p in ps:
            for q in ps:
                if q<p: continue
                n=p*q
                if 2<=n<X//2: yield n

    sp=sorted(set(semiprimes()))
    print("   semiprimes p*q: %d moduli, drift %.3f (grows like "
          "(loglog)^2)" % (len(sp), sum(1.0/n for n in sp)))
    s1=survivor_trend(((n,rng.randrange(n)) for n in sp), X)
    show("semiprimes, random residues", s1, "unbounded")

    # squarefree with >= 2 factors, odd only (no common prime)
    odd_sp=[n for n in sp if n%2]
    print("   odd semiprimes: %d moduli, drift %.3f"
          % (len(odd_sp), sum(1.0/n for n in odd_sp)))
    s2=survivor_trend(((n,rng.randrange(n)) for n in odd_sp), X)
    show("odd semiprimes, random", s2, "unbounded, no prime 2")

    # aligned versions - alignment is what saved the evens
    s3=survivor_trend(((n,0) for n in odd_sp), X)
    show("odd semiprimes, ALIGNED", s3, "unbounded")
    print()
    live=[(n,s) for n,s in (("semiprimes random",s1),
                            ("odd semiprimes random",s2),
                            ("odd semiprimes aligned",s3))
          if s[-1] > 0.02 and s[-1] >= s[-2]-1e-9]
    if live:
        print("[verdict] CANDIDATE HARD INSTANCE: %s sustains "
              "survivors with unbounded prime support and "
              "divergent drift. Needs nu computed properly "
              "before any claim." % ", ".join(n for n,_ in live))
    else:
        print("[verdict] every unbounded-prime-support system "
              "tested loses its survivors. Consistent with the "
              "conjecture that nu > 0 plus divergent drift "
              "forces a finite prime cover - which would make "
              "the open case EMPTY. Measured, not proved.")


if __name__=="__main__":
    main()
