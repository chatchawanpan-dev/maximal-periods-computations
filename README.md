# Programs for "Linear recurrences over finite fields with the maximal number of periods"

C. Panraksa and N. Somasunderam, *Linear recurrences over finite fields with the maximal number of periods*, manuscript, 2026.

These programs check the examples, counts and constants of the paper. They need Python 3 (standard library; `mpmath` for
`certify_constants.py` and `density_constants.py`, `numpy` for
`threshold_and_maxdensity.py`) and one C program. Run each script from this
directory, e.g. `python3 lemma21_check.py`. Compile the cycle program first:
`gcc -O2 -o cycles cycles.c`. The `logs/` folder holds the output of our runs.

| Statement in the paper | Program |
|---|---|
| Lemma 2.1 (cycle structure from the divisors of f), all small fields | `lemma21_check.py` |
| Theorem 1.1, Theorem 1.2, Proposition 2.3 (invariance), k = 3, 4, state by state | `check_korder.py`, `check_korder_c.py` (uses `cycles`) |
| Theorem 1.2 and Proposition 3.5 by enumerating root tuples; Examples 3.4 and 3.6 | `check_rootenum.py` |
| Remark 3.8 (prescribed root orders) and the case omega(r) = k of Theorem 1.2 | `check_prescribed_orders.py` |
| Section 4: Corollary 4.1, Examples 4.2 and 3.6, Remark 4.3 (q = 17), Corollary 5.3, every q <= 101 | `check_k2_bruteforce.py 101` |
| Theorem 1.3 and Corollary 5.3: fixed-trace counts, the value 0.8317... at q = 1723 | `fixed_trace_numerics.py 2000`, `general_k_trace.py` |
| Proposition 5.4 (trace 0) | `trace_zero_check.py` |
| Theorem 1.4: certified enclosures of c_3 and c_4 | `certify_constants.py` (about 75 s) |
| Theorem 1.4: c_k to high precision | `density_constants.py` |
| Section 6: largest pi_2(q), pi_3(q) for prime q < 1.2e7 | `threshold_and_maxdensity.py` |
| Introduction, Example 3.4 (F_211 cubics with trace 3), Remark 4.3 (q = 5) | `check_added_examples.py` |

`gf.py` (small finite fields) and `formulas.py` (the closed formulas of the
paper, in exact arithmetic) are shared modules.

Some comments use earlier names: Theorem A is Theorem 1.2, Theorem B and
"Theorem 3.2" are Proposition 3.5, P_k(N) and Pi_N are theta_k(N), F is
Theta_k, nu_k(q) is pi_k(q), the cutoff y is z, "Theorems A and B" in
`check_korder.py` refer to Theorems 1.1 and 1.2, and "the t = 0 remark" is
Proposition 5.4.

## License

MIT; see `LICENSE`.
