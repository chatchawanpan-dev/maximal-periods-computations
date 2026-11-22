"""(a) First prime q for which the bound (3/2) D(N) sqrt(q) of Theorem 1.3 (k = 2)
is smaller than the mean N theta_2(N)/2.  (b) Largest proportions nu_2(q), nu_3(q)
over primes q below a limit."""
import numpy as np, sys
from math import sqrt

LIM = int(sys.argv[1]) if len(sys.argv) > 1 and sys.argv[0].endswith("threshold_and_maxdensity.py") else 12_000_000
spf = np.zeros(LIM + 1, dtype=np.int32)
for i in range(2, int(LIM ** 0.5) + 1):
    if spf[i] == 0:
        spf[i * i::i][spf[i * i::i] == 0] = i
isprime = (spf == 0); isprime[:2] = False

def fac(n):
    f = {}
    while n > 1:
        p = int(spf[n]) or n
        e = 0
        while n % p == 0:
            n //= p; e += 1
        f[p] = e
    return f

def theta(N, k):
    f = fac(N)
    s = 0.0
    from math import comb
    for j in range(k + 1):
        pr = 1.0
        for l, e in f.items():
            kap = (l - 1) * (1 - l ** (-k * e)) / (l ** k - 1)
            pr *= 1 - j * kap
        s += (-1) ** j * comb(k, j) * pr
    return s

if __name__ == "__main__":
    first = None
    best2 = (0, 0); best3 = (0, 0)
    for q in np.nonzero(isprime)[0]:
        q = int(q)
        if q < 5: continue
        N = q - 1
        f = fac(N)
        D = 1
        for e in f.values(): D *= 1 + 4 * e
        t2 = theta(N, 2)
        if first is None and 1.5 * D * sqrt(q) < N * t2 / 2:
            first = (q, f)
        nu2 = (N / q) * t2 / 2
        if nu2 > best2[0]: best2 = (nu2, q)
        if len(f) >= 3:
            nu3 = (N / q) ** 2 * theta(N, 3) / 6
            if nu3 > best3[0]: best3 = (nu3, q)
    print("first prime with nontrivial k=2 bound:", first)
    print("max nu_2 over primes <", LIM, ":", best2, fac(best2[1] - 1))
    print("max nu_3 over primes <", LIM, ":", best3, fac(best3[1] - 1))
