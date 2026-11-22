"""Direct enumeration of root tuples in a cyclic group C_N (written additively as
Z/N) with prescribed sum s, for checking Theorems A and B in ranges where brute
force on states is too slow.  A tuple (i_1, ..., i_k) with i_1 + ... + i_k = s
corresponds to roots g^{i_j} of a split polynomial with determinant g^s; distinct
unordered tuples give distinct polynomials.  The check computes the orders
N / gcd(i_j, N) and tests directly whether all 2^k - 1 subset lcms are distinct.

Usage: python3 check_rootenum.py
"""
import numpy as np
from math import gcd
from itertools import combinations
from collections import Counter
from gf import divisors
from formulas import T_max, T_max_L


def subset_lcm_injective(ords):
    """ords: array of shape (M, k). Returns boolean mask: all subset lcms distinct."""
    M, k = ords.shape
    lcms = []
    for size in range(1, k + 1):
        for S in combinations(range(k), size):
            L = ords[:, S[0]].copy()
            for j in S[1:]:
                L = np.lcm(L, ords[:, j])
            lcms.append(L)
    A = np.stack(lcms, axis=1)  # (M, 2^k-1)
    A.sort(axis=1)
    return np.all(A[:, 1:] != A[:, :-1], axis=1), A[:, -1]


def count(N, k, s, chunk=2_000_000):
    ordtab = np.array([N // gcd(i, N) for i in range(N)], dtype=np.int64)
    total = 0
    byL = Counter()
    # enumerate (i_1, ..., i_{k-1}) as integers base N
    T = N ** (k - 1)
    for start in range(0, T, chunk):
        idx = np.arange(start, min(T, start + chunk), dtype=np.int64)
        cols = []
        x = idx.copy()
        for _ in range(k - 1):
            cols.append(x % N)
            x //= N
        last = (s - sum(cols)) % N
        cols.append(last)
        I = np.stack(cols, axis=1)
        ords = ordtab[I]
        ok, Lmax = subset_lcm_injective(ords)
        total += int(ok.sum())
        vals, cnts = np.unique(Lmax[ok], return_counts=True)
        for v, c in zip(vals, cnts):
            byL[int(v)] += int(c)
    from math import factorial
    f = factorial(k)
    assert total % f == 0
    return total // f, {L: c // f for L, c in byL.items()}


def check(N, k, rs=None):
    q = N + 1  # formulas only depend on N = q - 1
    for r in (divisors(N) if rs is None else rs):
        s = N // r  # an element of order r
        c, byL = count(N, k, s)
        assert c == T_max(q, r, k), (N, k, r, c, T_max(q, r, k))
        for L in divisors(N):
            assert byL.get(L, 0) == T_max_L(q, r, k, L), (N, k, r, L, byL.get(L, 0))
        if c:
            print(f"  N={N} k={k} r={r}: {c} polynomials; by largest period {dict(sorted(byL.items()))}")
    print(f"N={N}, k={k}: r in {'all divisors' if rs is None else rs} agree", flush=True)


if __name__ == "__main__":
    for N in (30, 60, 120, 180, 210):
        check(N, 3)
    check(210, 4, rs=[210, 105, 42])
    check(420, 4, rs=[420, 210])
    print("all checks passed")
