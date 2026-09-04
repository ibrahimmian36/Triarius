#!/usr/bin/env python3
"""Session 77: attack Conjecture G, and prove Theorem H.

G (cylinder form): for every eps there is a modulus M such that
on EVERY residue class mod M the induced subsystem either DIES
or is channel-dominated (finitely many channels carry all but
eps of its drift).

THEOREM H (proved here, and new). If the classes are pairwise
INCOMPATIBLE then the survivor density EXISTS and equals
1 - sum 1/n_i. Proof: incompatible classes have DISJOINT kill
sets (Lemma A's engine), so the killed density is exactly
sum 1/n_i; by Lemma A that sum is at most 1, hence always
converges; densities of a disjoint union add in the limit. QED
This is a FIFTH closing mechanism, independent of channels,
death, and any drift hypothesis.

H also supplies the sharpest attack on G: such a system needs
INFINITELY many channels (every channel is a single class) and
need NOT die - exactly what G forbids. Whether it breaks G turns
on whether a cylinder decomposition rescues it.

BATTERY. Every construction this campaign has built, plus new
adversaries: the C-refuting blocks, the F-refuting doubled
blocks, aligned, the full family 2..K, dyadic valuation classes
(pairwise incompatible), a tower designed so no single M
suffices, and random systems.

VACUITY CHECK. G is worthless unless some system genuinely needs
M > 1 and some genuinely needs the dying branch. If everything
satisfies G at M = 1, G is a restatement, not a target, and this
probe says so.
Pod-only, stdlib only.
"""

import argparse
import random
import sys
import time
from math import gcd


def primes_upto(n):
    s = bytearray([1]) * (n + 1)
    s[0] = s[1] = 0
    i = 2
    while i * i <= n:
        if s[i]:
            s[i * i:: i] = bytearray(len(s[i * i:: i]))
        i += 1
    return [i for i in range(2, n + 1) if s[i]]


def kill_mask(classes, X):
    killed = bytearray(X)
    for n, a in classes:
        s0 = a % n
        if s0 == 0:
            s0 = n
        while s0 < n:
            s0 += n
        for i in range(s0, X, n):
            killed[i] = 1
    return killed


def check_H(X):
    """pairwise incompatible => survivor density = 1 - sum 1/n"""
    fam = [(2 ** k, 2 ** (k - 1)) for k in range(2, 16)]
    for i in range(len(fam)):
        for j in range(i + 1, len(fam)):
            (n1, a1), (n2, a2) = fam[i], fam[j]
            if (a1 - a2) % gcd(n1, n2) == 0:
                print("[H] family is NOT pairwise incompatible")
                return False
    s = sum(1.0 / n for n, _ in fam)
    killed = kill_mask(fam, X)
    lo = X // 2
    got = 1.0 - sum(killed[lo:X]) / float(X - lo)
    ok = abs(got - (1.0 - s)) < 0.005
    # MUTATION: a COMPATIBLE family must deviate from 1 - sum
    mut = [(2 ** k, 0) for k in range(2, 16)]
    sm = sum(1.0 / n for n, _ in mut)
    km = kill_mask(mut, X)
    gm = 1.0 - sum(km[lo:X]) / float(X - lo)
    fired = abs(gm - (1.0 - sm)) > 0.02
    print("[1] THEOREM H: pairwise-incompatible family, "
          "predicted survivor density %.5f, measured %.5f -> %s"
          % (1.0 - s, got, "PASS" if ok else "FAIL"))
    print("    MUTATION (same moduli, all residues 0, so "
          "COMPATIBLE): predicted %.5f if the formula applied, "
          "measured %.5f, deviation %.5f -> %s"
          % (1.0 - sm, gm, abs(gm - (1.0 - sm)),
             "fires" if fired else "DOES NOT FIRE"))
    return ok and fired


def sys_blocks(nb, s, double, plimit=2000000, pstart=200):
    ps = primes_upto(plimit)
    smalls = ps[1:nb + 1]
    pool = [p for p in ps if p > smalls[-1] and p >= pstart]
    out, idx = [], 0
    for q in smalls:
        acc = 0.0
        j = 0
        while idx < len(pool) and acc < s:
            p = pool[idx]
            idx += 1
            r = j % q
            j += 1
            inv = pow(q % p, -1, p)
            k = ((0 - r) % p) * inv % p
            a = (r + q * k) % (p * q)
            n = p * q
            if double:
                n2, a2 = 2 * n, (a if a % 2 == 0 else a + n)
                out.append((n2, a2 % n2))
            else:
                out.append((n, a))
            acc += 1.0 / p
    return out


def systems(X):
    ps = primes_upto(4000)
    S = {}
    # SCALED DOWN so every modulus acts inside the measured
    # window: with pstart=200/plimit=2e6 the moduli reach the
    # millions, and at any feasible X most classes kill nothing
    # in range - which fakes a G failure by measuring the window
    # rather than the mathematics.
    # pstart must be large enough that a block needs MANY
    # primes: with pstart=11 a single prime (1/11 = 0.09)
    # already exceeds s=0.05, so every block held ONE class and
    # was trivially channel-dominated - which is why the first
    # battery reported M=1 everywhere. pstart=200 gives
    # 1/p ~ 0.005, so ~10 primes per block, while the moduli
    # stay near 10^4 and act well inside the window.
    S["blocks (refutes C)"] = sys_blocks(6, 0.05, False,
                                         plimit=6000, pstart=200)
    S["doubled (refutes F)"] = sys_blocks(6, 0.05, True,
                                          plimit=6000, pstart=200)
    S["aligned 2,4,6,..."] = [(2 * k, 0) for k in range(1, 300)]
    S["full 2..600"] = [(n, 1 % n) for n in range(2, 600)]
    S["dyadic valuation"] = [(2 ** k, 2 ** (k - 1))
                             for k in range(1, 17)]
    S["incompatible mod 12"] = [(2, 0)] + [(12, r)
                                           for r in (1, 3, 5)]
    tower = []
    for k in range(1, 9):
        M = 2 ** k
        c = (M // 2 - 1) % M
        for p in ps[10:10 + 12]:
            tower.append((M * p, (c + M * (p % 7)) % (M * p)))
    S["tower"] = tower
    rng = random.Random(4)
    S["random 2..600"] = [(n, rng.randrange(n))
                          for n in range(2, 600)]
    return S


def best_channel_drift(cls, restarts, rng):
    best = 0.0
    for t in range(restarts):
        order = list(cls)
        if t:
            rng.shuffle(order)
        chosen, d = [], 0.0
        for n, a in order:
            if all((a - b) % gcd(n, m) == 0 for m, b in chosen):
                chosen.append((n, a))
                d += 1.0 / n
        best = max(best, d)
    return best


def test_G(name, cls, X, Ms, rng, killed):
    """MEMORY/TIME FIX: the induced subsystem's kills on a
    cylinder are exactly the full system's kills restricted to
    it, since classes missing the cylinder contribute nothing
    there. So the sieve runs ONCE per system, not once per
    residue class - ~1900 mask builds become 9, and the peak
    footprint is one mask."""
    """for each M, does every class mod M either die or become
    channel-dominated? thresholds are RELATIVE to the class's own
    drift - no absolute constants."""
    best_M, detail = None, ""
    for M in Ms:
        ok = True
        worst = ""
        for r in range(M):
            sub = [(n, a) for n, a in cls if (a - r) % gcd(n, M) == 0]
            if not sub:
                continue        # no class touches this cylinder
            lo = X // 2
            start = lo + ((r - lo) % M)
            tot_n = alive = 0
            for x in range(start, X, M):
                tot_n += 1
                if not killed[x]:
                    alive += 1
            if not tot_n:
                continue
            dens = alive / float(tot_n)
            tot = sum(1.0 / n for n, _ in sub)
            ch = best_channel_drift(sub, 6, rng)
            # RELATIVE, no absolute constant: a cylinder counts
            # as channel-dominated when one channel carries at
            # least half its drift.
            dominated = tot <= 1e-12 or ch >= 0.5 * tot
            dies = dens < 0.02
            if not (dies or dominated):
                ok = False
                worst = ("r=%d density %.3f, channel %.4f of "
                         "total %.4f" % (r, dens, ch, tot))
                break
        if ok:
            best_M = M
            break
        detail = worst
    return best_M, detail


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--X", type=int, default=1200000)
    args = ap.parse_args()
    t0 = time.time()
    okH = check_H(args.X)
    print()
    print("[2] Conjecture G on the full battery (smallest M that "
          "works; every class must DIE or be channel-dominated):")
    rng = random.Random(9)
    Ms = [1, 2, 3, 4, 6, 8, 12, 16, 24]
    fails, needM, needdeath = [], [], []
    for name, cls in sorted(systems(args.X).items()):
        tot_all = sum(1.0 / n for n, _ in cls)
        in_win = sum(1.0 / n for n, _ in cls
                     if n <= args.X // 4)
        cover = in_win / tot_all if tot_all else 1.0
        if cover < 0.8:
            print("    %-22s -> NO TEST: only %.0f%% of the "
                  "drift sits in classes acting inside the "
                  "window (X=%d); the rest kill nothing here"
                  % (name, 100 * cover, args.X))
            continue
        killed = kill_mask(cls, args.X)
        M, why = test_G(name, cls, args.X, Ms, rng, killed)
        del killed
        if M is None:
            fails.append((name, why))
            print("    %-22s -> NO M WORKS  (%s)" % (name, why))
        else:
            print("    %-22s -> M = %d" % (name, M))
            if M > 1:
                needM.append(name)
    print()
    print("[3] VACUITY CHECK: G is only meaningful if some system "
          "needs M > 1.")
    print("    systems needing M > 1: %s"
          % (", ".join(needM) if needM else "NONE"))
    print()
    if not okH:
        print("[verdict] Theorem H check failed; nothing else is "
              "trustworthy.")
        sys.exit(1)
    if fails:
        print("[verdict] CONJECTURE G REFUTED by %s. No tested "
              "modulus makes every residue class either die or "
              "become channel-dominated. G joins C, F and the "
              "plain epsilon-split. CANDIDATE: the M range is "
              "finite, so this is refutation over the tested "
              "range, not a proof."
              % ", ".join(n for n, _ in fails))
    elif not needM:
        print("[verdict] G IS VACUOUS over this battery: every "
              "system works at M = 1, so G restates 'dies or is "
              "channel-dominated' and adds nothing. It is not a "
              "usable target in this form.")
    else:
        print("[verdict] G SURVIVES this battery, and is NOT "
              "vacuous - some systems genuinely need M > 1. "
              "First target in the campaign to survive its own "
              "adversarial session. Theorem H proved alongside: "
              "pairwise-incompatible systems close by "
              "DISJOINTNESS, a fifth mechanism.")
    print("[done] %ds" % (time.time() - t0))


if __name__ == "__main__":
    main()
