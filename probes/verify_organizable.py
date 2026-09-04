#!/usr/bin/env python3
"""Session 117: which systems are coprime-layered at all?

(T2) forces every class modulus to be d_C * q with q prime and
d_C the modulus of a node; the tree condition forces any two
nodes to be nested or disjoint.  Two residue classes of moduli
d, d' always intersect when d, d' are coprime (CRT), and neither
contains the other unless one modulus divides the other.  So:

  CRITERION.  A choice of nodes is a tree iff every two chosen
  nodes that intersect have one modulus dividing the other.

This kills more systems than it looks.  We test:
  T  the depth-4 tree of session 116     (expected: organizable)
  R  Prop irred2, sparse protectors      (paper CLAIMS organizable)
  S  semiprimes pq                       (expected: not)
  M  multiples of a fixed d              (expected: organizable)
"""
import math, sys
from itertools import combinations


def primes_upto(n):
    s=bytearray([1])*(n+1); s[0]=s[1]=0
    for i in range(2,int(n**.5)+1):
        if s[i]: s[i*i::i]=bytearray(len(s[i*i::i]))
    return [i for i in range(n+1) if s[i]]
P=primes_upto(3000)


def factor(n):
    f={}; m=n
    for p in P:
        if p*p>m: break
        while m%p==0: f[p]=f.get(p,0)+1; m//=p
    if m>1: f[m]=f.get(m,0)+1
    return f


def options(n,a):
    """possible nodes (d, c) for the class: d = n/q, q prime, q not | d"""
    out=[]
    for q,e in factor(n).items():
        if e!=1: continue          # (T4): q must not divide d = n/q
        d=n//q
        out.append((d,a%d))
    return out


def compatible(x,y):
    (d,c),(d2,c2)=x,y
    g=math.gcd(d,d2)
    if (c-c2)%g: return True                 # disjoint
    return d%d2==0 or d2%d==0                # intersect: need nesting


def organizable(cls, limit=200):
    """search for a tree choice; exact for small systems (backtracking
    with the most-constrained class first)."""
    opts=[options(n,a) for n,a in cls[:limit]]
    if any(not o for o in opts): return False, "a class has no legal node (no prime to exactly one power)"
    order=sorted(range(len(opts)), key=lambda i: len(opts[i]))
    chosen={}
    def bt(k):
        if k==len(order): return True
        i=order[k]
        for nd in opts[i]:
            if all(compatible(nd,v) for v in chosen.values()):
                chosen[i]=nd
                if bt(k+1): return True
                del chosen[i]
        return False
    return bt(0), "backtracking over %d classes"%len(opts)


def main():
    # T: session-116 tree, depth 2 slice: nodes 1, 2, 10 ... classes d*q
    T=[(2*q,q%2*1+2*(q%3)) for q in P[2:12]]+[(10*q,4+10*(q%5)) for q in P[12:20]]
    # R: Prop irred2 shape
    R=[(pj*q,1+pj*(q%pj)) for pj in (11,23,43) for q in P[10:18] if q!=pj]
    # S: semiprimes
    S=[(p*q,1) for i,p in enumerate(P[:6]) for q in P[i+1:i+4]]
    # M: multiples of 6
    M=[(6*q,6) for q in P[3:15]]
    for cls,lab,exp in ((T,"T depth-2 tree slice",True),
                        (R,"R Prop irred2 (paper says yes)",True),
                        (S,"S semiprimes pq",False),
                        (M,"M multiples of 6",True)):
        ok,why=organizable(cls)
        flag="" if ok==exp else "   <-- CONTRADICTS EXPECTATION"
        print("  %-32s organizable: %-5s  (%d classes)%s"%(lab,ok,len(cls),flag))
    print()
    print("  reason, when it fails: nodes of coprime moduli always meet")
    print("  (CRT) and never nest, so all chosen nodes must lie on one")
    print("  divisor chain.")
    return 0


if __name__=="__main__":
    sys.exit(main())
