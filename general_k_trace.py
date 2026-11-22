"""Fixed-trace counts for order k (sum of roots prescribed) and the density
constants c_k.  T_k^tr(t) = #{monic f of degree k, roots in F_q (distinct), sum t,
with 2^k-1 periods} computed by enumerating root tuples in a prime field.
Compares with the main term ((N^k - (-1)^k)/q) P_k(N) / k! and the bound
q^{(k-1)/2} sum_{u=1}^k binom(k,u) prod (1 + 2 u e).
"""
import numpy as np
from math import gcd, comb, factorial, sqrt, prod
from fractions import Fraction
from gf import factor_int, lcm


def p_local(l, e, k):
    return Fraction((l - 1) * (l ** (k * e) - 1), l ** (k * e) * (l ** k - 1))


def P_k(N, k):
    f = factor_int(N)
    s = Fraction(0)
    for j in range(k + 1):
        pr = Fraction(1)
        for l, e in f.items():
            pr *= 1 - j * p_local(l, e, k)
        s += (-1) ** j * comb(k, j) * pr
    return s


def bound(N, q, k):
    f = factor_int(N)
    return q ** ((k - 1) / 2) * sum(comb(k, u) * prod(1 + 2 * u * e for e in f.values()) for u in range(1, k + 1))


def orders_mod_p(p):
    N = p - 1
    g = next(g for g in range(2, p) if all(pow(g, N // l, p) != 1 for l in factor_int(N)))
    o = np.zeros(p, dtype=np.int64)
    x = 1
    for i in range(N):
        o[x] = N // gcd(i, N)
        x = x * g % p
    return o


def private_ok(O):
    # O: (M, k) array of orders; True where every column has a private prime
    M, k = O.shape
    ok = np.ones(M, dtype=bool)
    for i in range(k):
        rest = np.ones(M, dtype=np.int64)
        for j in range(k):
            if j != i:
                rest = np.lcm(rest, O[:, j])
        ok &= (rest % O[:, i]) != 0
    return ok


def trace_counts_k3(p):
    o = orders_mod_p(p)
    N = p - 1
    res = {}
    a = np.arange(1, p)
    A, B = np.meshgrid(a, a, indexing='ij')
    A = A.ravel(); B = B.ravel()
    for t in range(p):
        C = (t - A - B) % p
        m = C != 0
        O = np.stack([o[A[m]], o[B[m]], o[C[m]]], axis=1)
        cnt = int(private_ok(O).sum())
        assert cnt % 6 == 0
        res[t] = cnt // 6
    return res


if __name__ == "__main__":
    # density constants
    from itertools import count
    def primes(n):
        s = bytearray([1]) * (n + 1); s[0] = s[1] = 0
        for i in range(2, int(n ** .5) + 1):
            if s[i]: s[i*i::i] = bytearray(len(s[i*i::i]))
        return [i for i in range(n + 1) if s[i]]
    PR = primes(10 ** 6)
    for k in (3, 4, 5):
        c = 0.0
        for j in range(k + 1):
            pr = 1.0
            for l in PR:
                pr *= 1 - j * (l - 1) / (l ** k - 1)
            c += (-1) ** j * comb(k, j) * pr
        print(f"c_{k} ~ {c:.6f};  c_k/k! ~ {c / factorial(k):.6f}")
    # fixed trace k = 3
    for p in (31, 61, 151, 181, 211):
        N = p - 1
        res = trace_counts_k3(p)
        main = Fraction(N ** 3 + 1, p) * P_k(N, 3) / 6
        bd = bound(N, p, 3) / 6
        dev = max(abs(res[t] - float(main)) for t in range(1, p))
        total = sum(res.values())
        exact_total = Fraction(N ** 3) * P_k(N, 3) / 6
        assert total == exact_total, (p, total, exact_total)
        assert dev <= bd
        print(f"p={p}: main {float(main):9.2f}  range [{min(res[t] for t in range(1,p))}, {max(res[t] for t in range(1,p))}]  max dev {dev:7.2f}  bound {bd:10.1f}  T(0)={res[0]}  total ok")
