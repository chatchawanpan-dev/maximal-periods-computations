"""c_k = sup_N P_k(N) = probability that k indices are all covered when each prime l
independently covers one uniformly chosen index with probability k p_l,
p_l = (l-1)/(l^k-1).  Computed by a DP over primes (no cancellation), with the
tail beyond X bounded by k * sum_{odd n > X} n^(1-k) <= k / (2 (k-2) (X-1)^(k-2)) for k >= 3 which can
only increase the covering probability.
"""
from fractions import Fraction
from mpmath import mp, mpf
mp.dps = 40

def primes(n):
    s = bytearray([1]) * (n + 1); s[0] = s[1] = 0
    for i in range(2, int(n ** .5) + 1):
        if s[i]: s[i*i::i] = bytearray(len(s[i*i::i]))
    return [i for i in range(n + 1) if s[i]]

X = 2_000_000
PR = primes(X)
for k in (2, 3, 4, 5):
    dp = [mpf(0)] * (k + 1); dp[0] = mpf(1)
    for l in PR:
        p = mpf(l - 1) / (mpf(l) ** k - 1)
        new = [mpf(0)] * (k + 1)
        for n in range(k + 1):
            if dp[n] == 0: continue
            new[n] += dp[n] * (1 - (k - n) * p)      # no new index (none chosen, or an old one)
            if n < k:
                new[n + 1] += dp[n] * (k - n) * p
        dp = new
    tail = (mpf(k) / (2 * (k - 2) * mpf(X - 1) ** (k - 2))) if k > 2 else mpf("inf")  # >= k * sum over odd n > X of n^(1-k) >= sum_{l > X} k p_l
    print(f"k={k}: c_k(truncated) = {mp.nstr(dp[k], 12)}, tail mass <= {mp.nstr(tail, 3)}, c_k/k! = {mp.nstr(dp[k]/mp.factorial(k), 8)}")
