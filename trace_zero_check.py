"""Check the t = 0 remark: for k >= 2,
   |k! T_k^tr(0) - M_0 theta_k(N)| <= B_k(N) q^{k/2},
   M_0 = ((q-1)^k - (-1)^k)/q + (-1)^k,
and T_2^tr(0) = 0.  Prime fields only; k = 2, 3 by full enumeration, k = 4 for q = 211."""
import numpy as np
from math import comb, prod, factorial
from fractions import Fraction
from gf import factor_int
from general_k_trace import P_k, orders_mod_p, private_ok


def B(N, k):
    f = factor_int(N)
    return sum(comb(k, u) * prod(1 + 2 * u * e for e in f.values()) for u in range(1, k + 1))


def T_at_zero(p, k):
    o = orders_mod_p(p)
    a = np.arange(1, p, dtype=np.int64)
    cnt = 0
    if k == 2:
        A = a; Bv = (-A) % p
        O = np.stack([o[A], o[Bv]], axis=1)
        cnt = int(private_ok(O).sum())
    elif k == 3:
        A, Bv = np.meshgrid(a, a, indexing='ij'); A = A.ravel(); Bv = Bv.ravel()
        C = (-A - Bv) % p; m = C != 0
        cnt = int(private_ok(np.stack([o[A[m]], o[Bv[m]], o[C[m]]], axis=1)).sum())
    elif k == 4:
        A, Bv = np.meshgrid(a, a, indexing='ij'); A = A.ravel(); Bv = Bv.ravel()
        for c in range(1, p):
            D = (-A - Bv - c) % p; m = D != 0
            Cc = np.full(int(m.sum()), c, dtype=np.int64)
            cnt += int(private_ok(np.stack([o[A[m]], o[Bv[m]], o[Cc], o[D[m]]], axis=1)).sum())
    assert cnt % factorial(k) == 0
    return cnt // factorial(k)


if __name__ == "__main__":
    for k, ps in ((2, (31, 61, 211, 421)), (3, (31, 61, 151, 181, 211, 421)), (4, (211,))):
        for p in ps:
            N = p - 1
            T0 = T_at_zero(p, k)
            M0 = Fraction((p - 1) ** k - (-1) ** k, p) + (-1) ** k
            main = M0 * P_k(N, k)
            dev = abs(factorial(k) * T0 - main)
            bd = B(N, k) * p ** (k / 2)
            assert dev <= bd, (k, p)
            if k == 2:
                assert T0 == 0
            print(f"k={k} q={p}: T(0)={T0}  k!T(0)={factorial(k)*T0}  main={float(main):.2f}  "
                  f"dev={float(dev):.2f}  dev/q^(k/2)={float(dev)/p**(k/2):.3f}  bound/q^(k/2)={B(N,k)}")
    print("t = 0 checks passed")
