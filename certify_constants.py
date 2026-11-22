"""Certified enclosures of c_3 and c_4 (Theorem 1.4).

F((s_l)_{l <= y}) is computed with the recursion over primes
    P'_m = P_m (1 - (k - m) s_l) + P_{m-1} (k - m + 1) s_l,   0 <= m <= k,
starting from P_0 = 1, where P_m is the probability that exactly m given indices
are covered; then F = P_k.  The recursion is run in interval arithmetic
(mpmath.iv, 60 bits), so the output interval contains the exact value of
F((s*_l)_{l <= y}).  The primes l > y raise the covering probability by at most
k * sum_{l > y} s*_l <= k * sum_{odd n > y} n^{1-k} <= k / (2 (k - 2) (y - 1)^{k-2})
(Lemma 6.1), which gives an enclosure of c_k.
"""
from mpmath import iv, mpf, mp
import sys, math

iv.prec = 60


def primes(n):
    s = bytearray([1]) * (n + 1); s[0] = s[1] = 0
    for i in range(2, int(n ** .5) + 1):
        if s[i]:
            s[i * i::i] = bytearray(len(s[i * i::i]))
    return [i for i in range(n + 1) if s[i]]


def enclose(k, y, PR):
    P = [iv.mpf(1)] + [iv.mpf(0)] * k
    for l in PR:
        s = iv.mpf(l - 1) / (iv.mpf(l) ** k - 1)
        new = [P[m] * (1 - (k - m) * s) for m in range(k + 1)]
        for m in range(1, k + 1):
            new[m] += P[m - 1] * (k - m + 1) * s
        P = new
    lo = P[k].a
    tail = iv.mpf(k) / (2 * (k - 2) * iv.mpf(y - 1) ** (k - 2))
    hi = (P[k] + tail).b
    return lo, hi, P[k]


if __name__ == "__main__":
    y = int(sys.argv[1]) if len(sys.argv) > 1 else 2_000_000
    PR = primes(y)
    print(f"y = {y}, pi(y) = {len(PR)}")
    for k in (3, 4):
        lo, hi, Pk = enclose(k, y, PR)
        fk = math.factorial(k)
        print(f"k={k}: F_y in [{mp.nstr(Pk.a, 15)}, {mp.nstr(Pk.b, 15)}] (width {mp.nstr(Pk.b - Pk.a, 3)})")
        print(f"      c_k in [{mp.nstr(lo, 15)}, {mp.nstr(hi, 15)}];  c_k/k! < {mp.nstr(hi / fk, 10)}")
    print("checks: c_3/3! < 0.00123 and c_4/4! < 6e-8 are asserted below")
    lo3, hi3, _ = enclose(3, y, PR)
    lo4, hi4, _ = enclose(4, y, PR)
    assert hi3 / 6 < mpf("0.00123") and hi4 / 24 < mpf("6e-8")
    assert mpf("0.007365") <= lo3 and hi3 < mpf("0.007366")
    assert mpf("1.3478e-6") <= lo4 and hi4 < mpf("1.3479e-6")
    print("stated digits c_3 = 0.007365..., c_4 = 1.3478...e-6 certified")
