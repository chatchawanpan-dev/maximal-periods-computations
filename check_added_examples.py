"""Checks for the statements added in v5.1 (proof review of 2026-11-22).
1. Over F_211, delta = 7: the cubics X^3-3X^2+144X-7 and X^3-3X^2+65X-7 have
   trace 3, determinant 7 and seven periods (root orders 6,5,7 and 14,21,35);
   T_3(7) from Theorem 1.2, and the number of distinct traces among those cubics.
2. Over F_5, X^2-X-1 = (X-3)^2 and ord(3) = 4 = q-1; the Fibonacci recurrence
   mod 5 has the periods 4 and 20.
3. Roots of orders 6, 10, 15 (e.g. in F_31) give only the four periods 6,10,15,30.
Prime fields only; cycle lengths are computed state by state.
Usage: python3 check_added_examples.py
"""
from itertools import combinations, product
from fractions import Fraction
from math import gcd, comb

def lcm(a, b): return a * b // gcd(a, b)

def order(x, p):
    n, y = 1, x % p
    while y != 1:
        y = y * x % p; n += 1
    return n

def primes_of(n):
    out, d = [], 2
    while d * d <= n:
        if n % d == 0:
            out.append(d)
            while n % d == 0: n //= d
        d += 1
    if n > 1: out.append(n)
    return out

def periods(p, coeffs):
    """coeffs = (c_{k-1},...,c_0) of u_{j+k} = c_{k-1}u_{j+k-1}+...+c_0 u_j;
    returns {cycle length: number of cycles} on nonzero states."""
    k = len(coeffs)
    seen, out = set(), {}
    for s in product(range(p), repeat=k):
        if not any(s) or s in seen: continue
        cur, n = s, 0
        while True:
            seen.add(cur)
            nxt = sum(c * x for c, x in zip(coeffs, cur)) % p
            cur = (nxt,) + cur[:-1]; n += 1
            if cur == s: break
        out[n] = out.get(n, 0) + 1
    return out

def private_ok(ms):
    k = len(ms)
    return len({__import__('functools').reduce(lcm, [ms[i] for i in S], 1)
                for r in range(1, k + 1) for S in combinations(range(k), r)}) == 2 ** k - 1

# 1. F_211
p, delta = 211, 7
assert order(delta, p) == 210
for roots in [(15, 55, 144), (63, 178, 184)]:
    t = sum(roots) % p; d = roots[0] * roots[1] * roots[2] % p
    s = (roots[0]*roots[1] + roots[0]*roots[2] + roots[1]*roots[2]) % p
    ms = [order(x, p) for x in roots]
    print("roots", roots, "trace", t, "det", d, "s", s, "orders", ms, "private", private_ok(ms))
    assert t == 3 and d == 7 and private_ok(ms)
assert {(15*55 + 15*144 + 55*144) % p, (63*178 + 63*184 + 178*184) % p} == {144, 65}
# all cubics with determinant 7 and seven periods
N = p - 1
orders = {x: order(x, p) for x in range(1, p)}
cubics = []
for a in range(1, p):
    for b in range(a + 1, p):
        c = delta * pow(a * b, -1, p) % p
        if c <= b: continue
        if private_ok([orders[a], orders[b], orders[c]]):
            cubics.append((a, b, c))
traces = {}
for a, b, c in cubics:
    traces[(a + b + c) % p] = traces.get((a + b + c) % p, 0) + 1
k = 3; r = 210
w = {l: Fraction(1, l ** (1 + 0)) for l in primes_of(r)}  # N/r = 1
T = Fraction(N ** (k - 1), 6) * sum((-1) ** j * comb(k, j) *
     __import__('math').prod([1 - j * w[l] ** (k - 1) for l in w]) for j in range(k + 1))
mult = {}
for v in traces.values(): mult[v] = mult.get(v, 0) + 1
print("T_3(7) formula", T, "enumerated", len(cubics), "distinct traces", len(traces), "multiplicities", mult)
assert T == len(cubics)
# 2. F_5
p = 5
assert (pow(3, 2, p), (-6) % p, 9 % p) == (4, 4, 4)  # (X-3)^2 = X^2 - 6X + 9 = X^2 - X - 1 mod 5
assert order(3, 5) == 4
print("Fibonacci mod 5 cycles:", periods(5, (1, 1)))
assert set(periods(5, (1, 1))) == {4, 20}
# 3. orders 6, 10, 15 in F_31
p = 31
roots = [next(x for x in range(1, p) if order(x, p) == m) for m in (6, 10, 15)]
a, b, c = roots
coeffs = ((a + b + c) % p, -(a*b + a*c + b*c) % p, (a * b * c) % p)
per = periods(p, coeffs)
print("orders 6,10,15 roots", roots, "periods", sorted(per))
assert sorted(per) == [6, 10, 15, 30]
print("all checks passed")
