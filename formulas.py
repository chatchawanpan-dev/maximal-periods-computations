"""Closed formulas from the revised manuscript, in exact rational arithmetic."""
from fractions import Fraction
from math import comb, gcd, factorial
from gf import factor_int, divisors, lcm


def mobius(n):
    f = factor_int(n)
    if any(e > 1 for e in f.values()):
        return 0
    return -1 if len(f) % 2 else 1


def vl(n, l):
    e = 0
    while n % l == 0:
        n //= l
        e += 1
    return e


def primes(n):
    return sorted(factor_int(n).keys())


def omega_weights(N, r):
    """w_l = l^{-1-v_l(N/r)} for primes l | r."""
    return {l: Fraction(1, l ** (1 + vl(N // r, l))) for l in primes(r)}


def T_max(q, r, k):
    """Number of monic degree-k f over F_q with det = delta (ord delta = r)
    having 2^k - 1 distinct periods (Theorem A)."""
    N = q - 1
    w = omega_weights(N, r)
    s = Fraction(0)
    for j in range(k + 1):
        prod = Fraction(1)
        for l, wl in w.items():
            prod *= 1 - j * wl ** (k - 1)
        s += (-1) ** j * comb(k, j) * prod
    val = Fraction(N ** (k - 1), factorial(k)) * s
    assert val.denominator == 1, (q, r, k, val)
    return int(val)


def E_set(L, r):
    return [l for l in primes(r) if vl(L, l) == vl(r, l)]


def M_k(L, r, k):
    val = Fraction(L ** (k - 1))
    for l in primes(L // r):
        val *= 1 - Fraction(1, l ** (k - 1))
    return val


def T_max_L(q, r, k, L):
    """Theorem B: count with 2^k-1 periods and largest period L."""
    N = q - 1
    if r == 0 or L % r or N % L:
        return 0
    E = E_set(L, r)
    s = Fraction(0)
    for j in range(k + 1):
        prod = Fraction(1)
        for l in E:
            prod *= 1 - Fraction(j, l ** (k - 1))
        s += (-1) ** j * comb(k, j) * prod
    val = M_k(L, r, k) / factorial(k) * s
    assert val.denominator == 1, (q, r, k, L, val)
    return int(val)


def census_k2(q, r, rho):
    """(T1, T2, T3) for k = 2 (Corollary)."""
    N = q - 1
    w = omega_weights(N, r)
    A = Fraction(1)
    B = Fraction(1)
    for wl in w.values():
        A *= 1 - wl
        B *= 1 - 2 * wl
    T1 = Fraction(q + 1, 1) + N * B - 2 * rho
    T1 /= 2
    T2 = rho + N * (A - B)
    T3 = Fraction(N, 2) * (1 - 2 * A + B)
    for x in (T1, T2, T3):
        assert x.denominator == 1
    return int(T1), int(T2), int(T3)


def pair_count(N, r, m, n):
    """Number of traces with period set {m, n, lcm(m,n)} (m, n incomparable)."""
    s = 0
    for d in divisors(m):
        for e in divisors(n):
            if lcm(d, e) % r == 0:
                s += mobius(m // d) * mobius(n // e) * gcd(d, e)
    return s


def Pi_N(N):
    """Probability that two independent uniform elements of C_N have orders
    incomparable under divisibility."""
    alpha = Fraction(1)
    beta = Fraction(1)
    for l, e in factor_int(N).items():
        b = Fraction(l - 1, l + 1) + Fraction(2, (l + 1) * l ** (2 * e))
        beta *= b
        alpha *= (1 + b) / 2
    return 1 - 2 * alpha + beta


def D_N(N):
    val = 1
    for l, e in factor_int(N).items():
        val *= 1 + 4 * e
    return val
