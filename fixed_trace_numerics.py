"""Fixed-trace counts n3tr(q, t) = #{delta : X^2 - t X + delta has three periods}
for prime q, computed directly from orders of the roots lambda, t - lambda.
(The equivalence between three periods and incomparable root orders is checked by
brute force in check_k2_bruteforce.py.)  Reports the largest deviation from the
mean N*Pi_N/2 and the theoretical bound 1.5 D(N) sqrt(q), and the traces with no
three-period determinant.
"""
import numpy as np
from math import gcd, sqrt
from formulas import Pi_N, D_N
from gf import factor_int


def primes_upto(n):
    s = bytearray([1]) * (n + 1)
    s[0] = s[1] = 0
    for i in range(2, int(n ** 0.5) + 1):
        if s[i]:
            s[i * i::i] = bytearray(len(s[i * i::i]))
    return [i for i in range(n + 1) if s[i]]


def orders_mod_p(p):
    N = p - 1
    g = next(g for g in range(2, p) if all(pow(g, N // l, p) != 1 for l in factor_int(N)))
    ordv = np.zeros(p, dtype=np.int64)
    x = 1
    for i in range(N):
        ordv[x] = N // gcd(i, N)
        x = x * g % p
    return ordv


def counts(p):
    ordv = orders_mod_p(p)
    lam = np.arange(1, p)
    out = np.zeros(p, dtype=np.int64)
    for t in range(1, p):
        mu = (t - lam) % p
        ok = mu != 0
        a = ordv[lam[ok]]
        b = ordv[mu[ok]]
        inc = (b % a != 0) & (a % b != 0)
        out[t] = int(inc.sum()) // 2
    return out


if __name__ == "__main__":
    import sys
    X = int(sys.argv[1]) if len(sys.argv) > 1 else 1000
    worst_ratio = 0.0
    zero_cases = []
    rows = []
    for p in primes_upto(X):
        if p < 5:
            continue
        N = p - 1
        c = counts(p)
        mean = float(Pi_N(N)) * N / 2
        dev = max(abs(c[t] - mean) for t in range(1, p))
        bound = 1.5 * D_N(N) * sqrt(p)
        assert dev < bound
        worst_ratio = max(worst_ratio, dev / sqrt(p))
        assert abs(c[1:].mean() - mean) < 1e-9  # exact average over t != 0
        if len(factor_int(N)) >= 2:
            z = [t for t in range(1, p) if c[t] == 0]
            if z:
                zero_cases.append((p, len(z)))
        rows.append((p, mean, c[1:].min(), c[1:].max(), dev))
    print(f"primes 5 <= p <= {X}: max deviation / sqrt(p) = {worst_ratio:.3f}")
    big = [r for r in rows if r[0] > X // 2]
    print("sample (p, mean, min_t, max_t, max|dev|):")
    for r in big[-8:]:
        print(f"  {r[0]:5d} {r[1]:9.2f} {r[2]:6d} {r[3]:6d} {r[4]:8.2f}")
    print(f"primes with omega(p-1) >= 2 having some t != 0 with no three-period delta: {len(zero_cases)}")
    print("largest such primes:", zero_cases[-10:])
