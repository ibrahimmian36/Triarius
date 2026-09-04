#!/usr/bin/env python3
"""Session 109 (corrected). The block construction cannot test
the question: its drift grows like s*loglog, reaching only 0.05
at eight small primes, so NOTHING dies and both the protected
and unprotected versions retain survivors near 0.97. My
calibration passed that on a '> 0.45' threshold the value
cleared trivially.

Test the same claim where the drift is LARGE. The question's
misaligned branch says: unbounded prime support, no common
factor, divergent drift, residues misaligned => survivors die.

  UNPROTECTED: odd moduli, random residues. No common factor,
               unbounded prime support, drift ~ (1/2) log X.
  PROTECTED:   the same moduli doubled, residues forced even, so
               every kill lands in the evens and the odds are
               structurally unreachable.

The pair differs ONLY by the protector, so the comparison
isolates it.
"""
import sys
from math import gcd
from functools import reduce


def survivors(gen, X):
    m=bytearray(X)
    for n,a in gen:
        st=a%n
        if st==0: st=n
        while st<n: st+=n
        if st<X: m[st::n]=b"\x01"*len(m[st::n])
    out=[]; Y=16384
    while Y<=X:
        out.append(1-sum(m[Y//2:Y])/float(Y-Y//2)); Y*=2
    return out


def main():
    X=4*10**6
    import random
    rng=random.Random(17)
    odd=[n for n in range(3, X//2, 2)]
    drift=sum(1.0/n for n in odd)
    print("odd moduli 3..%d : %d classes, drift %.3f "
          "(divergent), gcd = %d"
          % (X//2, len(odd), drift, reduce(gcd, odd)))
    print()
    res={n: rng.randrange(n) for n in odd[:60000]}
    sub=odd[:60000]
    su=survivors(((n,res[n]) for n in sub), X)
    print("  UNPROTECTED (odd moduli, random residues)")
    print("     survivors: %s" % "  ".join("%.4f"%v for v in su[-6:]))
    sp=survivors(((2*n, (2*res[n]) % (2*n)) for n in sub), X)
    print("  PROTECTED   (same, doubled, kills confined to evens)")
    print("     survivors: %s" % "  ".join("%.4f"%v for v in sp[-6:]))
    print()
    calib = abs(sp[-1]-0.5) < 0.08
    print("  calibration - protected sits near 1/2 (odds are "
          "unreachable): %s (%.4f)" % (calib, sp[-1]))
    if not calib:
        print("  [ABORT] calibration failed: protected version "
              "is at %.4f, not near 1/2" % sp[-1]); sys.exit(1)
    falling = all(b <= a+1e-9 for a,b in zip(su,su[1:]))
    print("  unprotected falling: %s (%.4f -> %.4f)"
          % (falling, su[0], su[-1]))
    print()
    if su[-1] < 0.05:
        print("[verdict] REMOVING THE PROTECTOR KILLS THE "
              "SURVIVORS: %.4f against %.4f for the protected "
              "twin. The misaligned branch behaves as the "
              "question requires - divergent misaligned drift "
              "with nu > 0 needed the common factor."
              % (su[-1], sp[-1]))
    else:
        print("[verdict] unprotected system RETAINS survivors "
              "(%.4f) - check against every theorem before "
              "calling it a counterexample." % su[-1])


if __name__=="__main__":
    main()
