#!/usr/bin/env python3
"""Theorem 4.4's 'up to a finite set' holds because c is an
INTEGER: the exceptional integers are the representatives of c
modulo each n, and for every n > c that representative is c
itself, so at most c+1 distinct values arise however many
classes there are.

MUTATION: with residues that genuinely vary with n the
exceptional set is unbounded. (The first mutation attempted here
used a_n = (7n+3) mod n, which is just 3 mod n - itself aligned
to c = 3 - so it tested nothing and did not fire.)
"""
Ns = [50, 200, 1000, 5000, 20000]


def exceptional(residue, N):
    """distinct integers below their class's truncation start"""
    return {residue(n) % n for n in range(2, N + 2)}


c = 37
print("aligned to the integer c = %d  (bound c+1 = %d):"
      % (c, c + 1))
sizes = []
for N in Ns:
    e = exceptional(lambda n: c, N)
    sizes.append(len(e))
    print("   %6d classes -> %3d distinct exceptional integers"
          % (N, len(e)))
flat = len(set(sizes)) == 1 and sizes[0] <= c + 1

print()
print("MUTATION - residues a_n = n-1, genuinely varying:")
msizes = []
for N in Ns:
    e = exceptional(lambda n: n - 1, N)
    msizes.append(len(e))
    print("   %6d classes -> %6d distinct exceptional integers"
          % (N, len(e)))
grows = all(y > x for x, y in zip(msizes, msizes[1:]))

print()
print("aligned set is bounded and constant: %s" % flat)
print("mutation set grows without bound:    %s" % grows)
print("-> %s" % ("integer alignment is LOAD-BEARING for the "
                 "finiteness claim" if (flat and grows)
                 else "CHECK FAILED"))
