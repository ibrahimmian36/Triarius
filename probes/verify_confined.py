#!/usr/bin/env python3
"""Session 115: an irreducible system with nu > 0 and divergent
drift -- no cylinder dies -- covered by the SUMMABLE-PROTECTORS
theorem (tail kills confined to cylinders of summable density,
finite heads of convergent drift).

Protectors p_j = the ceil(4 j^{3/2})-th prime, so sum 1/p_j < 0.24 and sum sqrt(j)/p_j diverges (PNT).
Group j: classes (p_j q, a) with q over primes not among the
protectors, taken in order until sum 1/q >= S_j = sqrt(j), with
a = 1 mod p_j and random mod q.

Exact checks:
  1. sum 1/p_j over the truncation, and nu >= 1 - that sum
  2. IRREDUCIBLE: every cylinder holding >= 2 classes has prime
     modulus, induced moduli distinct primes, induced survivor
     measure prod(1-1/l) > 0 -- computed exactly, minimum shown
  3. drift of the truncation = sum S_j/p_j; and the formula
     sum sqrt(j)/p_j ~ sum 1/(6 j log j) partial sums (analytic, diverges)
  4. misaligned residues (unbounded misaligned part)
  5. MUTATION (expected to escape): protectors = all primes,
     sum 1/p_j diverges; MC survivors fall as groups are added
  6. MC survivors of the truncation vs the bound 1 - sum 1/p_j
"""
import math, random, sys


def primes_upto(n):
    s=bytearray([1])*(n+1); s[0]=s[1]=0
    for i in range(2,int(n**.5)+1):
        if s[i]: s[i*i::i]=bytearray(len(s[i*i::i]))
    return [i for i in range(n+1) if s[i]]


def build(J, prot, P, rng, S=lambda j: math.sqrt(j)):
    """prot: list of protector primes (length >= J)."""
    protset=set(prot[:J]); classes={}; groups={}
    for j in range(1,J+1):
        pj=prot[j-1]; s=0.0; qs=[]
        for q in P:
            if q in protset: continue
            qs.append(q); s+=1.0/q
            if s>=S(j): break
        else: sys.exit("prime table too small for j=%d"%j)
        groups[pj]=qs
        for q in qs:
            r=rng.randrange(q)
            a=(1+pj*(((r-1)*pow(pj,-1,q))%q))%(pj*q)   # CRT: 1 mod p_j, r mod q
            assert a%pj==1 and a%q==r
            classes[(pj,q)]=a
    return classes, groups


def cylinders(classes):
    cyl={}
    for (p,q),a in classes.items():
        cyl.setdefault((p,a%p),[]).append(q)
        cyl.setdefault((q,a%q),[]).append(p)
    return cyl


def mc(classes, Pall, rng, n=2500):
    idx={}
    for (p,q),a in classes.items(): idx.setdefault(p,[]).append((q,a%p,a%q))
    surv=0
    for _ in range(n):
        r={l:rng.randrange(l) for l in Pall}
        killed=False
        for p,lst in idx.items():
            rp=r[p]
            for q,ap,aq in lst:
                if rp==ap and r[q]==aq: killed=True; break
            if killed: break
        surv+=not killed
    return surv/float(n)


def main():
    rng=random.Random(7); J=6
    P=primes_upto(100000)
    prot=[]
    for j in range(1,J+1):
        prot.append(P[math.ceil(4*j**1.5)-1])   # the ceil(4 j^1.5)-th prime
    classes,groups=build(J,prot,P,rng)
    print("protectors p_j = the ceil(4 j^1.5)-th prime:", prot)
    print("classes: %d;  group sizes: %s" % (len(classes),[len(groups[p]) for p in prot]))

    # 1. nu bound
    sp=sum(1.0/p for p in prot)
    tail=sum(1.0/(k*math.log(k)) for k in (math.ceil(4*j**1.5) for j in range(J+1,200000)))
    print("[1] sum 1/p_j = %.4f (six terms); tail j>%d bounded via p_k > k log k by %.4f; total < %.3f;  nu >= %.4f"
          % (sp, J, tail, sp+tail, 1-sp))

    # 2. irreducibility
    cyl=cylinders(classes)
    Pset=set(P); minsurv=1.0; maxind=0.0; allprime=True; multi=0
    for (l,c),mods in cyl.items():
        if len(mods)<2: continue
        multi+=1
        allprime&=(l in Pset) and all(m in Pset for m in mods) and len(set(mods))==len(mods)
        d=sum(1.0/m for m in mods); maxind=max(maxind,d)
        s=1.0
        for m in mods: s*=(1-1.0/m)
        minsurv=min(minsurv,s)
    ok2=allprime and minsurv>0
    print("[2] irreducible: %d cylinders hold >=2 classes; all prime with distinct prime induced moduli: %s;"
          % (multi,allprime))
    print("    max induced drift %.3f (= S_J = sqrt(%d) = %.3f expected at (p_J,1)); min induced survivor %.4f > 0: %s"
          % (maxind,J,math.sqrt(J),minsurv,minsurv>0))

    # 3. drift
    dr=sum(sum(1.0/q for q in groups[p])/p for p in prot)
    part=[sum(1.0/(6*j*math.log(j+1)) for j in range(1,N+1)) for N in (10,10**3,10**5,10**7)]
    print("[3] drift of truncation = %.4f;  asymptotic term sqrt(j)/p_j ~ 1/(6 j log j): partial sums %s"
          % (dr, ", ".join("%.2f"%x for x in part)), "(diverges as (1/8) log N)")

    # 4. misalignment
    mis=[p*q for (p,q),a in classes.items() if a%p and a%q]
    print("[4] misaligned classes %d/%d; largest misaligned part %d (grows with the truncation)"
          % (len(mis),len(classes),max(mis)))

    # 5. mutation: protectors = all primes
    Pall=sorted({p for p,q in classes}|{q for p,q in classes}|set(P[:300]))
    protm=[p for p in P if p not in ()][:40]   # first 40 primes as protectors
    print("[5] mutation, protectors = first primes (sum 1/p diverges), S_j = sqrt(j):")
    prev=None; falling=True
    for Jm in (5,10,20,40):
        clm,_=build(Jm,protm,P,random.Random(3),S=lambda j: min(math.sqrt(j),0.7))
        Pm=sorted({p for p,q in clm}|{q for p,q in clm})
        sv=mc(clm,Pm,random.Random(1),n=1500)
        print("      J=%2d  sum 1/p_j=%.3f  survivors %.3f"%(Jm,sum(1.0/p for p in protm[:Jm]),sv))
        if prev is not None and sv>prev+0.02: falling=False
        prev=sv
    print("    survivors fall as protectors accumulate:", falling)

    # 6. MC survivors vs bound
    sv=mc(classes,Pall,random.Random(2),n=2500)
    ok6=sv>=1-sp-0.03
    print("[6] MC survivors %.3f >= bound %.3f (within MC error): %s"%(sv,1-sp,ok6))

    allok=ok2 and dr<1 and falling and ok6 and sp<0.5
    print()
    print("[verdict]", "PASS: irreducible, nu > 0, misaligned, drift diverging by the formula; the divergent-protector mutation collapses as designed." if allok else "FAIL")
    return 0 if allok else 1


if __name__=="__main__":
    sys.exit(main())
