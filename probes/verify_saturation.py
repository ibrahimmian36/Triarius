#!/usr/bin/env python3
"""Session 121: the saturation reduction, exact on finite systems.

Work modulo a fixed L; every class is a cylinder (d, c) with d | L.
A = residues mod L killed by no class.  A cylinder is NULL if it
meets A in no residue; the reduction P = the maximal null
cylinders under inclusion.  Checks:
  1. every class lies inside a member of P
  2. union of P = union of the classes (exactly, mod L)
  3. reduction is idempotent: P(P) = P
  4. saturation: for every member of P and every proper divisor
     d' of its modulus, the coarsened cylinder meets A
Controls:
  MP  many-protector shape (protectors 5, 9 = 3^2; inside, all
      sub-cylinders): must collapse to the protectors
  SAT a system that is already saturated must be fixed
  RND random systems: 1-4 must hold on every one of them
"""
import random, sys
from math import gcd


def divisors(L): return [d for d in range(1,L+1) if L%d==0]


def contains(D,C,d,c):   # is (d,c) inside (D,C)?  need D | d and c = C mod D
    return d%D==0 and c%D==C%D


def reduce_system(cls, L):
    killed=set()
    for d,c in cls:
        for x in range(c%d,L,d): killed.add(x)
    A=set(range(L))-killed
    null=[(d,c) for d in divisors(L) for c in range(d)
          if all(x not in A for x in range(c,L,d))]
    maximal=[(d,c) for (d,c) in null
             if not any((D,C)!=(d,c) and contains(D,C,d,c) for (D,C) in null)]
    return A, sorted(set(maximal)), killed


def check(cls, L, label):
    A,P,killed=reduce_system(cls,L)
    c1=all(any(contains(D,C,d,c) for (D,C) in P) for (d,c) in cls)
    unionP=set()
    for d,c in P:
        for x in range(c%d,L,d): unionP.add(x)
    c2=(unionP==killed)
    A2,P2,_=reduce_system(P,L)
    c3=(P2==P and A2==A)
    c4=all(any(x in A for x in range(c%d2,L,d2))
           for (d,c) in P for d2 in divisors(d) if d2<d)
    print("  %-34s classes %3d -> reduced %3d  contains %s  union %s  idempotent %s  saturated %s"
          % (label,len(cls),len(P),c1,c2,c3,c4))
    return c1 and c2 and c3 and c4, len(cls), len(P)


def main():
    rng=random.Random(121); L=2*2*3*3*5*7   # 1260
    ok=True
    # MP: protectors (5,1) and (9,1), filled with all their sub-cylinders mod L
    mp=[]
    for (p,r) in ((5,1),(9,1)):
        for d in divisors(L):
            if d%p==0 and d>p:
                for c in range(d):
                    if c%p==r: mp.append((d,c))
    print("  ambient modulus L = %d (finite systems modulo L)" % L)
    g,n0,n1=check(mp,L,"MP many-protector (fills 5,9)"); ok&=g and n1==2
    # SAT: coprime primes with one residue each - already maximal
    sat=[(2,1),(3,2),(5,4),(7,3)]
    g,n0,n1=check(sat,L,"SAT coprime primes"); ok&=g and n1==4
    # RND
    for k in range(6):
        cls=sorted({(d,rng.randrange(d)) for d in rng.sample([x for x in divisors(L) if x>1],12)})
        g,_,_=check(cls,L,"RND %d"%k); ok&=g
    print()
    print("[verdict]", "PASS: maximal null cylinders contain every class, reproduce the union, are idempotent and saturated; many-protector collapses to its protectors." if ok else "FAIL")
    return 0 if ok else 1


if __name__=="__main__":
    sys.exit(main())
