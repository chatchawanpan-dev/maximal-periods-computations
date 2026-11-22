"""Checks for recurrences of order k >= 3.

Part 1 (brute force, all fields q <= QMAX): for every monic f of degree k with
f(0) != 0, follow the companion map on F_q^k \\ {0} and record the set of cycle
lengths.  Check  #P(f) <= 2^k - 1, the structure of P(f) for split f with
distinct roots (P = set of subset lcms of root orders), the characterisation of
equality, the counts of Theorems A and B, and the invariance in ord(det).

Part 2 (brute force, one determinant): the same for a single delta in a larger
prime field, e.g. q = 31, k = 3, delta primitive.

Usage: python3 check_korder.py
"""
import itertools
from collections import Counter, defaultdict
from math import gcd
from gf import field, prime_powers, divisors, lcm
from formulas import T_max, T_max_L


def companion_cycles(F, coeffs):
    """coeffs = (c_{k-1}, ..., c_0); recurrence u_{j+k} = sum c_i u_{j+i}.
    State (u_{j+k-1}, ..., u_j) encoded in base q, most significant first."""
    q = F.q
    k = len(coeffs)
    add, mul = F.add, F.mul
    Q = q ** k
    top = q ** (k - 1)
    seen = bytearray(Q)
    seen[0] = 1
    cyc = Counter()
    for s in range(1, Q):
        if seen[s]:
            continue
        cur = s
        length = 0
        while not seen[cur]:
            seen[cur] = 1
            # digits
            x = cur
            digs = []
            for _ in range(k):
                digs.append(x % q)
                x //= q
            # digs[0] = u_j, ..., digs[k-1] = u_{j+k-1}
            new = 0
            for i in range(k):
                new = add[new][mul[coeffs[k - 1 - i]][digs[i]]]
            cur = new * top + cur // q
            length += 1
        cyc[length] += 1
    return cyc


def roots_of(F, coeffs):
    """Roots (with multiplicity) of f = X^k - sum c_i X^i if f splits, else None."""
    q = F.q
    k = len(coeffs)
    # polynomial as list low->high: f = X^k - c_{k-1} X^{k-1} - ... - c_0
    poly = [F.neg[c] for c in reversed(coeffs)] + [1]
    roots = []
    for _ in range(k):
        found = None
        for x in range(1, q):
            # Horner
            v = 0
            for a in reversed(poly):
                v = F.add[F.mul[v][x]][a]
            if v == 0:
                found = x
                break
        if found is None:
            return None
        roots.append(found)
        # synthetic division by (X - found)
        n = len(poly) - 1
        newp = [0] * n
        acc = 0
        for i in range(n, 0, -1):
            acc = F.add[F.mul[acc][found]][poly[i]]
            newp[i - 1] = acc
        poly = newp
    return roots


def subset_lcms(orders):
    k = len(orders)
    out = set()
    for mask in range(1, 1 << k):
        L = 1
        for i in range(k):
            if mask >> i & 1:
                L = lcm(L, orders[i])
        out.add(L)
    return out


def private_primes_ok(orders):
    k = len(orders)
    for i in range(k):
        rest = 1
        for j in range(k):
            if j != i:
                rest = lcm(rest, orders[j])
        if rest % orders[i] == 0:
            return False
    return True


def run_field(q, k, deltas=None, verbose=True):
    F = field(q)
    N = q - 1
    maxP = 2 ** k - 1
    struct_by_order = defaultdict(list)
    deltas = range(1, q) if deltas is None else deltas
    for delta in deltas:
        r = F.order[delta]
        # det of companion = (-1)^(k+1) c_0  ->  c_0 = (-1)^(k+1) delta
        c0 = delta if (k + 1) % 2 == 0 else F.neg[delta]
        count_max = 0
        byL = Counter()
        structs = []
        hist = Counter()
        for rest in itertools.product(range(q), repeat=k - 1):
            coeffs = tuple(rest) + (c0,)
            cyc = companion_cycles(F, coeffs)
            P = set(cyc)
            assert sum(l * c for l, c in cyc.items()) == q ** k - 1
            hist[len(P)] += 1
            structs.append(tuple(sorted(cyc.items())))
            assert len(P) <= maxP
            rts = roots_of(F, coeffs)
            if rts is not None and len(set(rts)) == k:
                orders = [F.order[x] for x in rts]
                assert P == subset_lcms(orders), (q, coeffs, P, orders)
                assert (len(P) == maxP) == private_primes_ok(orders)
                # product of roots is delta
                prod = 1
                for x in rts:
                    prod = F.mul[prod][x]
                assert prod == delta
            if len(P) == maxP:
                assert rts is not None and len(set(rts)) == k
                count_max += 1
                byL[max(P)] += 1
        assert count_max == T_max(q, r, k), (q, k, delta, count_max, T_max(q, r, k))
        for L in divisors(N):
            assert byL[L] == T_max_L(q, r, k, L), (q, k, delta, L, byL[L], T_max_L(q, r, k, L))
        struct_by_order[r].append(tuple(sorted(structs)))
        if verbose:
            print(f"  q={q} k={k} delta={delta} (ord {r}): max-count {count_max}  #P histogram {dict(sorted(hist.items()))}")
    for r, lst in struct_by_order.items():
        assert all(x == lst[0] for x in lst), (q, k, r)


if __name__ == "__main__":
    import sys
    print("Part 1: k = 3, all prime powers q <= 13, all determinants")
    for q in prime_powers(13):
        run_field(q, 3, verbose=False)
        print(f"q={q}: ok")
    print("Part 1b: k = 4, q in {2,3,4,5}")
    for q in (2, 3, 4, 5):
        run_field(q, 4, verbose=False)
        print(f"q={q}: ok")
    print("Part 2: k = 3, q = 31, primitive delta = 3 and delta of order 15, 10, 6")
    F = field(31)
    ds = [3] + [next(x for x in range(1, 31) if F.order[x] == o) for o in (15, 10, 6)]
    run_field(31, 3, deltas=ds)
    print("all checks passed")
