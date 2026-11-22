"""Full check of Lemma 2.1: for every monic f of degree k with f(0) != 0 over small
fields, compare the cycle structure of C_f (followed state by state) with
sum over monic divisors d != 1 of f of Phi_q(d)/ord(d) cycles of length ord(d),
where divisors, ord(d) and Phi_q(d) are all computed by brute force."""
import itertools
from collections import Counter
from gf import field
from check_korder import companion_cycles

# polynomials: tuples of coefficients low -> high, monic


def pdivmod(F, A, B):
    A = list(A)
    db = len(B) - 1
    inv = F.inv[B[-1]]
    Q = [0] * max(len(A) - db, 1)
    for i in range(len(A) - 1, db - 1, -1):
        c = F.mul[A[i]][inv]
        if c:
            Q[i - db] = c
            for j in range(db + 1):
                A[i - db + j] = F.sub(A[i - db + j], F.mul[c][B[j]])
    R = A[:db] if db > 0 else [0]
    while len(R) > 1 and R[-1] == 0:
        R.pop()
    return Q, R


def is_zero(R):
    return all(x == 0 for x in R)


def pmul(F, A, B):
    out = [0] * (len(A) + len(B) - 1)
    for i, a in enumerate(A):
        if a:
            for j, b in enumerate(B):
                out[i + j] = F.add[out[i + j]][F.mul[a][b]]
    return out


def pgcd(F, A, B):
    A = list(A); B = list(B)
    while len(B) > 1 or B[0] != 0:
        _, R = pdivmod(F, A, B)
        A, B = B, R
    return A


def ord_poly(F, d):
    """least j >= 1 with X^j = 1 mod d (d monic, d(0) != 0, deg d >= 1)."""
    x = [0, 1]
    cur = pdivmod(F, x, d)[1] if len(d) > 2 else pdivmod(F, x, d)[1]
    j = 1
    while not (len(cur) == 1 and cur[0] == 1):
        cur = pdivmod(F, pmul(F, cur, [0, 1]), d)[1]
        j += 1
    return j


def phi_poly(F, d):
    n = len(d) - 1
    cnt = 0
    for h in itertools.product(range(F.q), repeat=n):
        h = list(h)
        if is_zero(h):
            continue
        while len(h) > 1 and h[-1] == 0:
            h.pop()
        g = pgcd(F, d, h)
        if len(g) == 1:
            cnt += 1
    return cnt


def predicted(F, f):
    k = len(f) - 1
    cyc = Counter()
    for deg in range(1, k + 1):
        for low in itertools.product(range(F.q), repeat=deg):
            d = list(low) + [1]
            if d[0] == 0:
                continue
            _, R = pdivmod(F, f, d)
            if is_zero(R):
                o = ord_poly(F, d)
                ph = phi_poly(F, d)
                assert ph % o == 0
                cyc[o] += ph // o
    return cyc


if __name__ == "__main__":
    for k, qs in ((2, (2, 3, 4, 5, 7, 8, 9)), (3, (2, 3, 4, 5, 7)), (4, (2, 3))):
        for q in qs:
            F = field(q)
            n = 0
            for coeffs in itertools.product(range(q), repeat=k):
                if coeffs[-1] == 0:
                    continue
                # f = X^k - c_{k-1}X^{k-1} - ... - c_0, low -> high
                f = [F.neg[c] for c in reversed(coeffs)] + [1]
                obs = companion_cycles(F, coeffs)
                pre = predicted(F, f)
                assert obs == pre, (q, k, coeffs, obs, pre)
                n += 1
            print(f"k={k} q={q}: Lemma 2.1 agrees for all {n} polynomials")
