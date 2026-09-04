#!/usr/bin/env python3
"""Session 116: ingredients of the coprime-layered tree theorem,
checked exactly on an explicit tree.

TREE. Depth-k nodes are cylinders modulo M_k = prod_{i<=k} m_i,
m_i = i^2 + 1 (2, 5, 10, 17). A node has m_k - 1 children (one
gap residue), so the depth measure prod(1 - 1/m_k) stays away
from 0. Node primes come from pools disjoint by depth and
coprime to 2, 5, 17, so along any branch all primes are
distinct (T2-T4). Each node is a black hole (its whole pool,
drift about 1/2 to 3/4) with probability 1/k^2, else light (one
prime). Residues are random: the tree is misaligned.

Checks:
  1. nu exactly, as sum over cells of measure x prod(1-1/q)
     along the branch, against Monte Carlo on Zhat
  2. first-crossing cylinders at level N are pairwise disjoint
     and their total measure equals mu{g > N} exactly
  3. the light system at level N has drift <= N
  4. dying implies killed: a branch with g = infinity would have
     survival 0; here every g is finite and the product positive
  5. window survivors against the profinite measure of the
     visible system at several scales (ILLUSTRATION ONLY)
"""
import math, random, sys


def primes_upto(n):
    s=bytearray([1])*(n+1); s[0]=s[1]=0
    for i in range(2,int(n**.5)+1):
        if s[i]: s[i*i::i]=bytearray(len(s[i*i::i]))
    return [i for i in range(n+1) if s[i]]


def build(K, rng):
    P=primes_upto(6000); tree_primes={2,5,17}
    pools={1:[p for p in P if 3<=p<=30], 2:[p for p in P if 31<=p<=200],
           3:[p for p in P if 201<=p<=1500], 4:[p for p in P if 1501<=p<=5800]}
    for k in pools: pools[k]=[p for p in pools[k] if p not in tree_primes]
    ms={k:k*k+1 for k in range(1,K+1)}
    # node = (depth, d, c) with cylinder c + d Z; gap cells likewise
    nodes=[]; cells=[]     # cells: (d, c, branch_primes) at the leaves of the truncation
    def rec(k, d, c, branch):
        # this node's classes
        if k>0:
            m=ms[k]; black=(rng.random()<1.0/(k*k)) if k>=2 else False
            Q=list(pools[k]) if black else [rng.choice(pools[k])]
            cls=[(d*q, (c + d*rng.randrange(q))%(d*q), q) for q in Q]
            nodes.append((k,d,c,Q,cls,black)); branch=branch+list(Q)
        if k==K:
            cells.append((d,c,tuple(branch))); return
        m=ms[k+1]
        for r in range(m-1): rec(k+1, d*m, c+d*r, branch)
        cells.append((d*m, c+d*(m-1), tuple(branch)))   # gap
    rec(0,1,0,[])
    return nodes, cells, ms


def main():
    rng=random.Random(19); K=4
    nodes,cells,ms=build(K,rng)
    classes=[cl for n in nodes for cl in n[4]]
    print("tree: %d nodes, %d cells, %d classes; black holes per depth: %s"
          % (len(nodes),len(cells),len(classes),
             {k:sum(1 for n in nodes if n[0]==k and n[5]) for k in range(1,K+1)}))
    Mk=1
    for k in ms: Mk*=ms[k]
    depthm=1.0
    for k in ms: depthm*=(1-1.0/ms[k])
    print("depth-%d cylinder measure prod(1-1/m_k) = %.4f; gaps %.4f" % (K,depthm,1-depthm))

    # 1. nu exactly and by MC
    nu=0.0; drift=sum(1.0/cl[0] for cl in classes)
    for d,c,br in cells:
        s=1.0
        for q in br: s*=(1-1.0/q)
        nu+=s/d
    print("[1] nu (exact branch product) = %.5f;  drift of truncation %.3f" % (nu,drift))
    allq=sorted({q for _,_,br in cells for q in br})
    surv=0; n=20000
    for _ in range(n):
        x=rng.randrange(Mk); r={q:rng.randrange(q) for q in allq}
        killed=False
        for dd,a,q in classes:
            if x%(dd//q)==a%(dd//q) and r[q]==a%q: killed=True; break
        surv+=not killed
    mc=surv/float(n); ok1=abs(mc-nu)<3*math.sqrt(nu*(1-nu)/n)+0.005
    print("    MC on Zhat (%d samples) = %.4f   agree: %s" % (n,mc,ok1))

    # 2. first-crossing cylinders vs mu{g > N}
    ok2=True
    print("[2] first-crossing cylinders at level N:")
    gs={}   # per node: partial sum along branch to that node
    for k,d,c,Q,cls,black in nodes:
        # ancestors: nodes containing this cylinder
        anc=[n for n in nodes if n[0]<k and c%n[1]==n[2]]
        gs[(d,c)]=sum(1.0/q for n in anc for q in n[3])+sum(1.0/q for q in Q)
    for N in (0.05,0.2,0.5,1.0):
        fc=[]
        for k,d,c,Q,cls,black in nodes:
            if gs[(d,c)]>N:
                anc=[n for n in nodes if n[0]<k and c%n[1]==n[2]]
                if all(gs[(n[1],n[2])]<=N for n in anc): fc.append((d,c))
        disjoint=all(not (c1%math.gcd(d1,d2)==c2%math.gcd(d1,d2)) for i,(d1,c1) in enumerate(fc) for (d2,c2) in fc[i+1:])
        mfc=sum(1.0/d for d,c in fc)
        mg=sum(1.0/d for d,c,br in cells if sum(1.0/q for q in br)>N)
        # cells with g>N: their branch sum exceeds N at some node, i.e. they lie in a first-crossing cylinder
        light_drift=sum(1.0/cl[0] for k,d,c,Q,cls,black in nodes if gs[(d,c)]<=N for cl in cls)
        good=disjoint and abs(mfc-mg)<1e-12 and light_drift<=N+1e-12
        ok2&=good
        print("      N=%.2f  %3d cylinders, disjoint %s, measure %.5f = mu{g>N} %.5f, light drift %.4f <= N: %s"
              % (N,len(fc),disjoint,mfc,mg,light_drift,good))

    # 4. survival positive on every finite branch
    minsurv=min(math.prod(1-1.0/q for q in br) for _,_,br in cells)
    print("[4] every branch finite: min branch survival %.4f > 0; a branch with g=inf would have product 0" % minsurv)

    # 5. window survivors vs visible profinite (illustration)
    print("[5] window [X/2,X] survivors s(X) vs nu_X of the visible system (illustration only):")
    X=4*10**6; killed=bytearray(X+1)
    for dd,a,q in classes:
        if dd<=X:
            start=a if a>=dd else a+dd*((dd-a+dd-1)//dd)
            if start<=X: killed[start::dd]=b"\x01"*len(range(start,X+1,dd))
    for Xs in (10**4,10**5,10**6,4*10**6):
        vis=[cl for cl in classes if cl[0]<=Xs]
        visq={cl[2] for cl in vis}
        nuX=0.0
        for d,c,br in cells:
            s=1.0
            for q in br:
                if q in visq: s*=(1-1.0/q)
            nuX+=s/d
        lo=Xs//2; s=sum(1 for x in range(lo,Xs+1) if not killed[x])/float(Xs-lo+1)
        print("      X=%8d  s(X)=%.4f  nu_X=%.4f  ratio %.3f" % (Xs,s,nuX,s/nuX))
    print()
    ok=ok1 and ok2 and minsurv>0
    print("[verdict]", "PASS: exact nu matches MC, first-crossing cylinders disjoint with measure mu{g>N}, light drift <= N." if ok else "FAIL")
    return 0 if ok else 1


if __name__=="__main__":
    sys.exit(main())
