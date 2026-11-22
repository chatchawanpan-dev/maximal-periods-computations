"""Check the prescribed-orders formula (Proposition on prescribed root orders):
the number of ordered k-tuples in C_N with product of order r (here: sum s = N/r in
Z/N) and orders exactly (m_1, ..., m_k) equals

    sum_{d_i | m_i} prod mu(m_i/d_i) * [r | lcm(d)] * prod(d_i) / lcm(d).

Also checks the special case omega(r) = k of Theorem A:
    T_k = ( (r / rad r) * (N/r)' )^(k-1),
where (N/r)' is the largest divisor of N/r coprime to r.
"""
import itertools
import numpy as np
from math import gcd, prod
from collections import Counter
from gf import divisors, lcm, factor_int
from formulas import mobius, T_max


def formula(r, ms):
    total = 0
    for ds in itertools.product(*[divisors(m) for m in ms]):
        L = 1
        for d in ds:
            L = lcm(L, d)
        if L % r:
            continue
        total += prod(mobius(m // d) for m, d in zip(ms, ds)) * prod(ds) // L
    return total


def enumerate_orders(N, k, s):
    ordtab = np.array([N // gcd(i, N) for i in range(N)], dtype=np.int64)
    idx = np.arange(N ** (k - 1), dtype=np.int64)
    cols = []
    x = idx.copy()
    for _ in range(k - 1):
        cols.append(x % N)
        x //= N
    cols.append((s - sum(cols)) % N)
    ords = np.stack([ordtab[c] for c in cols], axis=1)
    return Counter(map(tuple, ords.tolist()))


def coprime_part(n, r):
    for l in factor_int(r):
        while n % l == 0:
            n //= l
    return n


if __name__ == "__main__":
    for N, k in ((60, 2), (120, 2), (60, 3), (120, 3), (90, 3)):
        for r in divisors(N):
            cnt = enumerate_orders(N, k, N // r)
            for ms in itertools.product(divisors(N), repeat=k):
                assert cnt.get(ms, 0) == formula(r, ms), (N, k, r, ms, cnt.get(ms, 0), formula(r, ms))
        print(f"N={N} k={k}: prescribed-orders formula agrees for all r and all order tuples")
    # special case omega(r) = k
    for N in (30, 60, 120, 180, 210, 420, 2310, 4620):
        for k in (2, 3, 4, 5):
            for r in divisors(N):
                if len(factor_int(r)) == k:
                    rad = prod(factor_int(r))
                    val = ((r // rad) * coprime_part(N // r, r)) ** (k - 1)
                    assert T_max(N + 1, r, k) == val, (N, k, r)
    print("special case omega(r) = k agrees with Theorem A")
