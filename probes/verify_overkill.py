#!/usr/bin/env python3
"""Session 99: can intermediate-modulus classes OVER-kill?

The excess is at most the window's over-representation of the
killed set relative to its global density. Buchstab's example
runs the OTHER way - the prime sieve leaves MORE survivors in
[X/2,X] than globally, so it UNDER-kills, and under-killing
never creates an excursion.

So: can the intermediate range over-kill at all?  If not, then
with the spike bound of session 97 the whole excess is
controlled.

Sieve theory predicts it can, but barely: Buchstab's omega(u)
oscillates about e^{-gamma} with rapidly decaying amplitude, so
the window-to-global survivor ratio should cross 1 in BOTH
directions with shrinking deviation.

Method: take the classes to be PRIMES in a range, residue 0.
Coprimality makes the global survivor density the exact product
prod(1-1/p), so only the window needs sieving.

REGISTERED: OVER-KILLING OCCURS iff the ratio drops below 1 for
some range; UNDER-KILLING ONLY iff it stays above 1 everywhere -
which would be a structural finding and must be reported as
MEASURED, not proved.
"""
import math
import sys


def primes_between(lo, hi):
    s = bytearray([1])*(hi+1); s[0]=s[1]=0
    i = 2
    while i*i <= hi:
        if s[i]: s[i*i::i] = bytearray(len(s[i*i::i]))
        i += 1
    return [i for i in range(max(2,lo), hi+1) if s[i]]


def window_survivor_density(ps, X):
    """density of [X/2,X] avoiding all p in ps"""
    lo = X//2
    seg = bytearray([1])*(X-lo+1)
    for p in ps:
        start = ((lo+p-1)//p)*p
        for m in range(start, X+1, p):
            seg[m-lo] = 0
    return sum(seg)/float(X-lo+1)


def main():
    EG = math.exp(-0.5772156649)
    X = 8*10**6
    print("classes = primes in [X^a, X^b], residue 0; window "
          "[X/2, X];  X = %d" % X)
    print("global survivor density is the exact product "
          "prod(1-1/p)")
    print()
    rows = []
    for a, b in ((0.05,0.15),(0.10,0.25),(0.15,0.35),
                 (0.20,0.45),(0.25,0.50),(0.05,0.50)):
        lo = max(2, int(X**a)); hi = int(X**b)
        ps = primes_between(lo, hi)
        if not ps:
            continue
        glob = 1.0
        for p in ps:
            glob *= (1.0-1.0/p)
        win = window_survivor_density(ps, X)
        ratio = win/glob
        rows.append((a,b,len(ps),win,glob,ratio))
        print("  moduli in [X^%.2f, X^%.2f]  %5d primes   "
              "window %.5f  global %.5f  ratio %.4f%s"
              % (a,b,len(ps),win,glob,ratio,
                 "   <-- OVER-KILLS" if ratio < 1 else ""))
    print()
    ctrl_ps = primes_between(2, X//2)
    g = 1.0
    for p in ctrl_ps: g *= (1.0-1.0/p)
    w = window_survivor_density(ctrl_ps, X)
    print("  control, full prime sieve to X/2: ratio %.4f "
          "(expect ~e^gamma = %.4f)" % (w/g, 1/EG))
    print()
    over = [r for r in rows if r[5] < 1]
    if over:
        print("[verdict] OVER-KILLING OCCURS: the ratio drops "
              "below 1 for %d of %d ranges (min %.4f). "
              "Intermediate classes CAN create excursions, so "
              "the spike bound alone does not control the "
              "excess. Measured, not proved."
              % (len(over), len(rows),
                 min(r[5] for r in rows)))
    else:
        print("[verdict] NO OVER-KILLING in any tested range "
              "(min ratio %.4f). Intermediate classes only ever "
              "left MORE survivors than the global density. "
              "That is a structural finding and is reported as "
              "MEASURED, not proved - a function that oscillates "
              "may simply not have dipped in these ranges."
              % min(r[5] for r in rows))


if __name__ == "__main__":
    main()
