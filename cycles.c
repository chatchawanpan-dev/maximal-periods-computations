/* Cycle decomposition of the companion map of a recurrence over a PRIME field F_p.
 * u_{j+k} = c_{k-1} u_{j+k-1} + ... + c_0 u_j, state (u_{j+k-1},...,u_j) encoded
 * base p (most significant first).  Prints "length count" pairs for the cycles on
 * F_p^k \ {0}.
 * Usage: cycles p k c_{k-1} ... c_0
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>

int main(int argc, char **argv) {
    if (argc < 4) { fprintf(stderr, "usage\n"); return 1; }
    long p = atol(argv[1]);
    int k = atoi(argv[2]);
    if (argc != 3 + k) { fprintf(stderr, "need k coefficients\n"); return 1; }
    long c[16];
    for (int i = 0; i < k; i++) c[i] = atol(argv[3 + i]);  /* c[0] = c_{k-1}, ..., c[k-1] = c_0 */
    uint64_t Q = 1; for (int i = 0; i < k; i++) Q *= (uint64_t)p;
    uint64_t top = Q / (uint64_t)p;
    uint8_t *seen = calloc((Q + 7) / 8, 1);
    if (!seen) { fprintf(stderr, "oom\n"); return 1; }
    /* histogram of cycle lengths: store in a simple open-addressing table */
    #define HS 4096
    uint64_t hl[HS]; uint64_t hc[HS]; memset(hl, 0, sizeof hl); memset(hc, 0, sizeof hc);
    seen[0] |= 1;
    long dig[16];
    for (uint64_t s = 1; s < Q; s++) {
        if (seen[s >> 3] & (1u << (s & 7))) continue;
        uint64_t cur = s, len = 0;
        while (!(seen[cur >> 3] & (1u << (cur & 7)))) {
            seen[cur >> 3] |= (1u << (cur & 7));
            uint64_t x = cur;
            for (int i = 0; i < k; i++) { dig[i] = (long)(x % p); x /= p; }   /* dig[0]=u_j ... dig[k-1]=u_{j+k-1} */
            long nw = 0;
            for (int i = 0; i < k; i++) nw = (nw + c[k - 1 - i] * dig[i]) % p;
            cur = (uint64_t)nw * top + cur / p;
            len++;
        }
        uint64_t h = (len * 2654435761u) % HS;
        while (hl[h] && hl[h] != len) h = (h + 1) % HS;
        hl[h] = len; hc[h]++;
    }
    for (int h = 0; h < HS; h++) if (hl[h]) printf("%llu %llu\n", (unsigned long long)hl[h], (unsigned long long)hc[h]);
    free(seen);
    return 0;
}
