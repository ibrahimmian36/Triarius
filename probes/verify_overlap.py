#!/usr/bin/env python3
"""Session 119: overlapping (non-nested) protectors.

Protectors C_e = {x = 1 mod p_i, x = 1 mod p_j} over the pairs
e = {i,j} of a set of primes: moduli p_i p_j share factors
pairwise and never nest.

Exact quantities, by symmetric sums (never enumerating the
O(n^2) pairs of edges):
   A = sum 1/p, B = sum 1/p^2
   M = sum_e mu(C_e) = (A^2 - B)/2
   sum_{e,f} mu(C_e n C_f) = M                      (e = f)
        + sum_i (1/p_i)[(A-1/p_i)^2 - (B-1/p_i^2)]  (share one)
        + [M^2 - sum_e mu_e^2 - share]              (disjoint)
   R = that / M^2
   nu_prot = mu(x in no C_e) = mu(at most one i with x = 1 mod p_i)
           = P + sum_i (1/p_i)/(1-1/p_i) * P,  P = prod(1-1/p_i)

DENSE instance: all primes 3..N. Expect M -> inf, R bounded and
falling toward 1 (the decay is O(1/A) and A = sum 1/p grows like
log log N, so it is slow by nature), nu_prot -> 0.
SPARSE instance: least prime >= i^2. Expect M small, nu_prot
bounded below, so the summable-protectors theorem applies.
MC on Zhat checks nu_prot in both.
"""
import random, sys


def primes_upto(n):
    s=bytearray([1])*(n+1); s[0]=s[1]=0
    for i in range(2,int(n**.5)+1):
        if s[i]: s[i*i::i]=bytearray(len(s[i*i::i]))
    return [i for i in range(n+1) if s[i]]


def stats(ps):
    """all sums in closed form; no O(n^2) enumeration anywhere."""
    A=sum(1.0/p for p in ps); B=sum(1.0/p**2 for p in ps)
    C=sum(1.0/p**4 for p in ps)
    M=(A*A-B)/2.0                       # sum_e mu_e
    same=(B*B-C)/2.0                    # sum_e mu_e^2
    share=sum((1.0/p)*((A-1.0/p)**2-(B-1.0/p**2)) for p in ps)
    shareprod=sum((1.0/p**2)*((A-1.0/p)**2-(B-1.0/p**2)) for p in ps)
    disjoint=M*M-same-shareprod         # ordered disjoint pairs, mu_e mu_f
    R=(M+share+disjoint)/(M*M)
    P=1.0
    for p in ps: P*=(1-1.0/p)
    nu_prot=P*(1.0+sum((1.0/p)/(1-1.0/p) for p in ps))
    return A,M,R,nu_prot


def mc(ps, rng, n=20000):
    surv=0
    for _ in range(n):
        hits=sum(1 for p in ps if rng.randrange(p)==1%p)
        surv+= (hits<=1)
    return surv/float(n)


def main():
    rng=random.Random(31)
    out={}
    print("dense: all primes 3..N (protectors = all pairs)")
    for N in (10**3,10**4,10**5,10**6,10**7):
        ps=[p for p in primes_upto(N) if p>=3]
        A,M,R,nu=stats(ps)
        print("   N=%9d  %6d primes  A=%.3f  sum mu(C_e)=%8.3f  ratio R=%.4f  nu_prot=%.4f"
              % (N,len(ps),A,M,R,nu))
        out["dense%d"%N]=(len(ps),round(A,3),round(M,3),round(R,4),round(nu,4))
    ps=[p for p in primes_upto(1000) if p>=3]
    m=mc(ps,rng); print("   MC nu_prot at N=1000: %.4f" % m)
    out["dense_mc"]=round(m,4)

    print("sparse: p_i = least prime >= i^2, i = 2..I")
    P=primes_upto(20000)
    for I in (20,50):
        ps=[]
        for i in range(2,I+1):
            q=next((p for p in P if p>=i*i), None)
            if q and q not in ps: ps.append(q)
        A,M,R,nu=stats(ps)
        print("   I=%3d  %3d primes  A=%.3f  sum mu(C_e)=%.4f  ratio R=%.4f  nu_prot=%.4f"
              % (I,len(ps),A,M,R,nu))
        out["sparse%d"%I]=(len(ps),round(A,3),round(M,4),round(R,4),round(nu,4))
    ms=mc(ps,rng); print("   MC nu_prot at I=50: %.4f" % ms)
    out["sparse_mc"]=round(ms,4)

    ds=[out["dense%d"%N] for N in (10**3,10**4,10**5,10**6,10**7)]
    s=out["sparse50"]
    # R = 1 + 4(1+o(1))/A analytically; A = sum 1/p grows like
    # log log N, so R cannot approach 1 in any computation.  What
    # is checkable is the FIT: (R-1)A/4 falling toward 1.
    fit=[round((d[3]-1)*d[1]/4.0,4) for d in ds]
    print("   fit (R-1)A/4, should fall toward 1: %s" % fit)
    out["fit"]=fit
    ok = (all(b[3]<a[3] for a,b in zip(ds,ds[1:]))      # ratio falling
          and all(b<a for a,b in zip(fit,fit[1:]))      # fit falling
          and 1.0<fit[-1]<1.2                           # to the predicted 1
          and all(b[4]<a[4] for a,b in zip(ds,ds[1:]))  # nu_prot falling
          and abs(out["dense_mc"]-ds[0][4])<0.02        # MC matches exact
          and s[2]<1.0 and s[4]>0.5                     # sparse: summable, survivors
          and abs(out["sparse_mc"]-s[4])<0.02)
    print()
    print("[verdict]", "PASS: (R-1)A/4 falls to the predicted 1; dense protectors lose the survivors, sparse keep them with summable measure." if ok else "FAIL")
    out["ok"]=ok
    return out


if __name__=="__main__":
    r=main(); sys.exit(0 if r["ok"] else 1)
