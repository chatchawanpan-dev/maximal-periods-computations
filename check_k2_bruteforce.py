"""Brute-force check of every second-order statement in the revised manuscript.

For each prime power q <= QMAX, each delta in F_q^*, and each t in F_q, follow the
permutation (x, y) -> (t x - delta y, x) of F_q^2 \\ {0} and record the multiset of
cycle lengths.  Compare with the closed formulas in formulas.py.

Usage: python3 check_k2_bruteforce.py [QMAX]
"""
import sys
from collections import Counter, defaultdict
from gf import field, prime_powers, divisors, lcm
from formulas import (T_max, T_max_L, census_k2, pair_count, Pi_N, D_N)
from math import gcd, sqrt


def cycle_structure(F, t, delta):
    q = F.q
    add, mul, neg = F.add, F.mul, F.neg
    nd = neg[delta]
    seen = bytearray(q * q)
    seen[0] = 1
    cyc = Counter()
    for s in range(1, q * q):
        if seen[s]:
            continue
        length = 0
        cur = s
        while not seen[cur]:
            seen[cur] = 1
            x, y = divmod(cur, q)
            cur = add[mul[t][x]][mul[nd][y]] * q + x
            length += 1
        cyc[length] += 1
    return tuple(sorted(cyc.items()))


def main(QMAX):
    total_mats = 0
    total_pairs = 0
    fields = 0
    for q in prime_powers(QMAX):
        F = field(q)
        N = q - 1
        fields += 1
        structures_by_order = defaultdict(list)
        trace_three = Counter()  # t -> number of delta with 3 periods (fixed trace)
        for delta in range(1, q):
            r = F.order[delta]
            rho = sum(1 for x in range(1, q) if F.mul[x][x] == delta)
            counts = Counter()
            byL = Counter()
            bypair = Counter()
            structs = []
            for t in range(q):
                cs = cycle_structure(F, t, delta)
                structs.append(cs)
                lengths = [l for l, c in cs]
                P = len(lengths)
                counts[P] += 1
                total_mats += 1
                if P == 3:
                    L = max(lengths)
                    byL[L] += 1
                    m, n = sorted(lengths)[:2]
                    bypair[(m, n)] += 1
                    trace_three[t] += 1
                    # orbit counts for three periods: N/m, N/n, N^2/L
                    d = dict(cs)
                    assert d[m] == N // m and d[n] == N // n and d[L] == N * N // L, (q, t, delta, cs)
                assert P <= 3
            structures_by_order[r].append(tuple(sorted(structs)))
            total_pairs += 1
            # census
            T1, T2, T3 = census_k2(q, r, rho)
            assert (counts[1], counts[2], counts[3]) == (T1, T2, T3), (q, delta, counts, (T1, T2, T3))
            assert T3 == T_max(q, r, 2)
            # largest period
            for L in divisors(N):
                assert byL[L] == T_max_L(q, r, 2, L), (q, delta, L, byL[L], T_max_L(q, r, 2, L))
            # prescribed pair
            for m in divisors(N):
                for n in divisors(N):
                    if m < n and m % n and n % m:
                        pc = pair_count(N, r, m, n)
                        assert bypair[(m, n)] == pc, (q, delta, m, n, bypair[(m, n)], pc)
            # primitive determinant, odd q
            if q % 2 and r == N:
                from formulas import factor_int
                phi = N
                for l in factor_int(N):
                    phi = phi // l * (l - 1)
                assert (counts[1], counts[2], counts[3]) == ((q + 1) // 2, phi, N // 2 - phi)
        # invariance: the multiset of cycle structures depends only on ord(delta)
        for r, lst in structures_by_order.items():
            assert all(x == lst[0] for x in lst), (q, r)
        # fixed trace: compare with main term and explicit bound
        assert trace_three[0] == 0
        tot = sum(trace_three.values())
        # average identity: sum over delta of T3 equals N^2 Pi_N / 2
        assert tot == Pi_N(N) * N * N / 2, (q, tot, Pi_N(N) * N * N / 2)
        main = (q - 2) * Pi_N(N) / 2
        bound = 1.5 * D_N(N) * sqrt(q)
        worst = max(abs(trace_three[t] - main) for t in range(1, q))
        assert worst <= bound, (q, worst, bound)
        print(f"q={q:3d}: ok   fixed-trace max |T3tr - main| = {float(worst):7.2f}  (main {float(main):7.2f}, bound {bound:8.1f})")
    print(f"checked {fields} fields, {total_pairs} (q, delta) pairs, {total_mats} matrices")


if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 32)
