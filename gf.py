"""Small finite fields GF(p^n) with q <= ~256, elements encoded as integers 0..q-1.

Element x encodes the polynomial sum c_i X^i with x = sum c_i p^i.  Arithmetic is
done through addition and multiplication tables built from a fixed irreducible
polynomial found by brute force.  Only used for verification of the manuscript.
"""
from functools import lru_cache
from math import gcd


def factor_int(n):
    f = {}
    d = 2
    while d * d <= n:
        while n % d == 0:
            f[d] = f.get(d, 0) + 1
            n //= d
        d += 1
    if n > 1:
        f[n] = f.get(n, 0) + 1
    return f


def prime_power(q):
    f = factor_int(q)
    if len(f) != 1:
        return None
    (p, n), = f.items()
    return p, n


def divisors(n):
    ds = [1]
    for p, e in factor_int(n).items():
        ds = [d * p ** i for d in ds for i in range(e + 1)]
    return sorted(ds)


def lcm(a, b):
    return a // gcd(a, b) * b


class GF:
    def __init__(self, q):
        pp = prime_power(q)
        assert pp, q
        self.q = q
        self.p, self.n = pp
        p, n = self.p, self.n
        if n == 1:
            self.mod = None
        else:
            self.mod = self._find_irreducible()
        # tables
        self.add = [[self._add(a, b) for b in range(q)] for a in range(q)]
        self.mul = [[0] * q for _ in range(q)]
        for a in range(q):
            for b in range(a, q):
                c = self._mul(a, b)
                self.mul[a][b] = c
                self.mul[b][a] = c
        self.neg = [self._neg(a) for a in range(q)]
        self.inv = [0] * q
        for a in range(1, q):
            for b in range(1, q):
                if self.mul[a][b] == 1:
                    self.inv[a] = b
                    break
        # multiplicative orders
        self.order = [0] * q
        for a in range(1, q):
            x, k = a, 1
            while x != 1:
                x = self.mul[x][a]
                k += 1
            self.order[a] = k
        self.gen = next(a for a in range(1, q) if self.order[a] == q - 1)
        # discrete log w.r.t. gen
        self.log = [None] * q
        x = 1
        for i in range(q - 1):
            self.log[x] = i
            x = self.mul[x][self.gen]
        self.exp = [0] * (q - 1)
        for a in range(1, q):
            self.exp[self.log[a]] = a

    # polynomial helpers on coefficient lists (low degree first)
    def _digits(self, a):
        p, n = self.p, self.n
        return [(a // p ** i) % p for i in range(n)]

    def _undigits(self, ds):
        p = self.p
        return sum(c * p ** i for i, c in enumerate(ds))

    def _add(self, a, b):
        p = self.p
        return self._undigits([(x + y) % p for x, y in zip(self._digits(a), self._digits(b))])

    def _neg(self, a):
        p = self.p
        return self._undigits([(-x) % p for x in self._digits(a)])

    def _polymulmod(self, A, B, M):
        p = self.p
        prod = [0] * (len(A) + len(B) - 1)
        for i, x in enumerate(A):
            for j, y in enumerate(B):
                prod[i + j] = (prod[i + j] + x * y) % p
        # reduce modulo monic M of degree n
        n = len(M) - 1
        for d in range(len(prod) - 1, n - 1, -1):
            c = prod[d]
            if c:
                for i in range(n + 1):
                    prod[d - n + i] = (prod[d - n + i] - c * M[i]) % p
        res = prod[:n] + [0] * (n - len(prod[:n]))
        return res

    def _mul(self, a, b):
        p, n = self.p, self.n
        if n == 1:
            return (a * b) % p
        return self._undigits(self._polymulmod(self._digits(a), self._digits(b), self.mod))

    def _find_irreducible(self):
        p, n = self.p, self.n
        import itertools
        for coeffs in itertools.product(range(p), repeat=n):
            M = list(coeffs) + [1]
            if M[0] == 0:
                continue
            # irreducible iff no roots for n<=3; general test: X^(p^n) = X mod M and
            # gcd conditions; brute force: check no factor of degree <= n/2
            if self._is_irreducible(M):
                return M
        raise ValueError

    def _is_irreducible(self, M):
        import itertools
        p = self.p
        n = len(M) - 1
        for d in range(1, n // 2 + 1):
            for coeffs in itertools.product(range(p), repeat=d):
                D = list(coeffs) + [1]
                if self._polydivides(D, M):
                    return False
        return True

    def _polydivides(self, D, M):
        p = self.p
        R = M[:]
        dD = len(D) - 1
        for d in range(len(R) - 1, dD - 1, -1):
            c = R[d]
            if c:
                for i in range(dD + 1):
                    R[d - dD + i] = (R[d - dD + i] - c * D[i]) % p
        return all(x == 0 for x in R[:dD])

    def sub(self, a, b):
        return self.add[a][self.neg[b]]

    def pow(self, a, e):
        if a == 0:
            return 0 if e > 0 else 1
        return self.exp[(self.log[a] * e) % (self.q - 1)]


@lru_cache(maxsize=None)
def field(q):
    return GF(q)


def prime_powers(limit):
    return [q for q in range(2, limit + 1) if prime_power(q)]
