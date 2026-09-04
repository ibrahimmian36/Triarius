#!/usr/bin/env python3
"""Session 98: why the intermediate-modulus range resists.

Two crude bounds fail. (a) A union bound over the intermediate
classes gives about c*log X. (b) One might hope the window
[X/2, X] samples the union of those classes at its GLOBAL
density - but that union has period lcm ~ e^{X/2}, vastly longer
than the window, so there is no reason it should.

It genuinely does not, and the classical example shows it. Take
the classes to be the primes p <= X/2 with residue 0. Survivors
in [X/2, X] are exactly the PRIMES there, of density
~ 1/log X by the prime number theorem, while the global survivor
density of the same finite system is prod_{p<=X/2}(1-1/p)
~ e^{-gamma}/log(X/2) by Mertens. The window therefore
over-represents survivors by a factor tending to
e^gamma = 1.781... - the Buchstab phenomenon.

This is measured here before being claimed. Memory is capped:
last session's probe peaked at 602 MB and the user asked for
care.
"""
import math
import sys

FAIL = []
def ck(n, ok, d=""):
    print("  [%s] %s%s" % ("PASS" if ok else "FAIL", n,
                           ("  "+d) if d else ""))
    if not ok: FAIL.append(n)


def primes_upto(n):
    s = bytearray([1])*(n+1); s[0]=s[1]=0
    i = 2
    while i*i <= n:
        if s[i]: s[i*i::i] = bytearray(len(s[i*i::i]))
        i += 1
    return [i for i in range(2, n+1) if s[i]]


def main():
    print("window vs global survivor density for the sieve by")
    print("primes up to X/2, residue 0")
    print()
    EG = math.exp(0.5772156649)
    rows = []
    for X in (2*10**5, 10**6, 4*10**6, 1.6*10**7):
        X = int(X)
        half = X//2
        ps = primes_upto(half)
        # window survivors: integers in [X/2,X] with no prime
        # factor <= X/2 -- exactly the primes in (X/2, X]
        seg = bytearray([1])*(X-half+1)
        for p in ps:
            start = ((half+p-1)//p)*p
            for m in range(start, X+1, p):
                seg[m-half] = 0
        wsurv = sum(seg)/float(X-half+1)
        # global survivor density of the SAME finite system
        glob = 1.0
        for p in ps:
            glob *= (1.0 - 1.0/p)
        ratio = wsurv/glob if glob else float("nan")
        rows.append((X, wsurv, glob, ratio))
        print("  X=%-10d window %.3e   global %.3e   "
              "ratio %.4f" % (X, wsurv, glob, ratio))
        del seg, ps
    print()
    print("  e^gamma = %.4f" % EG)
    ratios = [r for _, _, _, r in rows]
    close = abs(ratios[-1] - EG) < 0.15
    converging = abs(ratios[-1]-EG) <= abs(ratios[0]-EG) + 0.02
    ck("window over-represents survivors by ~e^gamma",
       close and converging,
       "ratios %s" % ", ".join("%.4f" % r for r in ratios))
    ck("the discrepancy is a genuine constant factor, not a "
       "vanishing artifact", min(ratios) > 1.4,
       "min ratio %.4f" % min(ratios))
    print()
    if FAIL:
        print("FAILED: %s" % ", ".join(FAIL)); sys.exit(1)
    print("[verdict] window sampling fails for the intermediate "
          "range by a constant factor e^gamma, in the most "
          "natural example. The obstruction is real, not an "
          "artifact of the crude bounds.")


if __name__ == "__main__":
    main()
