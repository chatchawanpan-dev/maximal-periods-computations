"""State-level checks with the C cycle program (prime fields only):
 - k = 3, q = 31, every determinant;  k = 3, q = 61, determinants 2 (order 60) and 4 (order 30);
 - the quartic over F_211 with roots -1, 196, 188, 148 (determinant 2, primitive).
For every polynomial: #P <= 2^k - 1; #P = 2^k - 1 exactly for split f with private
primes, with N^{|S|}/L_S cycles of length L_S; counts per determinant match
Theorems 1.2 and 3.2."""
import itertools, subprocess
from collections import Counter
from gf import field, divisors, lcm
from formulas import T_max, T_max_L
from check_korder import roots_of, private_primes_ok

def cycles(q, coeffs):
    out = subprocess.run(["./cycles", str(q), str(len(coeffs))] + [str(c) for c in coeffs],
                         capture_output=True, text=True, check=True).stdout.split()
    return Counter({int(out[i]): int(out[i + 1]) for i in range(0, len(out), 2)})

def check(q, k, deltas):
    F = field(q); N = q - 1
    for delta in deltas:
        r = F.order[delta]
        c0 = delta if (k + 1) % 2 == 0 else F.neg[delta]
        cnt = 0; byL = Counter()
        for rest in itertools.product(range(q), repeat=k - 1):
            coeffs = tuple(rest) + (c0,)
            cyc = cycles(q, coeffs)
            assert sum(l * c for l, c in cyc.items()) == q ** k - 1
            P = set(cyc)
            assert len(P) <= 2 ** k - 1
            if len(P) == 2 ** k - 1:
                rts = roots_of(F, coeffs)
                o = [F.order[x] for x in rts]
                assert private_primes_ok(o)
                for mask in range(1, 1 << k):
                    L = 1
                    for i in range(k):
                        if mask >> i & 1: L = lcm(L, o[i])
                    assert cyc[L] == N ** bin(mask).count("1") // L
                cnt += 1; byL[max(P)] += 1
        assert cnt == T_max(q, r, k)
        for L in divisors(N):
            assert byL[L] == T_max_L(q, r, k, L)
        print(f"q={q} k={k} delta={delta} (ord {r}): {cnt} with 2^k-1 periods; by largest period {dict(byL)}", flush=True)

if __name__ == "__main__":
    check(31, 3, range(1, 31))
    check(61, 3, [2, 4])
    # the quartic over F_211
    F = field(211)
    roots = [210, 196, 188, 148]
    f = [1]
    for x in roots:
        f = [(a - x * b) % 211 for a, b in zip(f + [0], [0] + f)]
    coeffs = tuple((-c) % 211 for c in f[1:])
    cyc = cycles(211, coeffs)
    o = [F.order[x] for x in roots]
    assert sorted(o) == [2, 3, 5, 7]
    assert len(cyc) == 15
    for mask in range(1, 16):
        L = 1
        for i in range(4):
            if mask >> i & 1: L = lcm(L, o[i])
        assert cyc[L] == 210 ** bin(mask).count("1") // L
    print("F_211 quartic: 15 periods, cycle counts N^|S|/L_S confirmed:", sorted(cyc.items()))
