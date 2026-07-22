/* clawsweep.c — test every connected claw-free graph (graph6 on stdin) for
 * Schur positivity of its chromatic symmetric function. Exact integer
 * arithmetic throughout (int64/__int128 guards).
 *
 * Pipeline per graph: parse; claw-free filter; a_lam = #ordered stable
 * partitions with part sizes lam (mask DP); back-substitute K^T c = a in
 * lex-descending partition order (refines dominance; K unitriangular);
 * report if any Schur coefficient c_lam < 0.
 *
 * Usage: clawsweep <n> [kostka_data.txt path]     (graph6 lines on stdin)
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

#define MAXN 12
#define MAXP 128
static int N;                       /* graph order */
static uint16_t adj[MAXN];
static uint32_t FULL;

static int NP;                      /* number of partitions of N */
static int plen[MAXP], parts[MAXP][MAXN];
static long long K[MAXP][MAXP];

static void load_kostka(const char *path, int n) {
    FILE *f = fopen(path, "r");
    if (!f) { perror("kostka"); exit(1); }
    char tag[4]; int nn, np;
    while (fscanf(f, "%3s", tag) == 1) {
        if (!strcmp(tag, "N")) {
            if (fscanf(f, "%d %d", &nn, &np) != 2) exit(1);
            if (nn != n) { /* skip this block */
                char line[4096];
                fgets(line, sizeof line, f);
                for (int i = 0; i < 2 * np; i++) fgets(line, sizeof line, f);
            } else {
                NP = np;
                for (int i = 0; i < np; i++) {
                    if (fscanf(f, "%3s", tag) != 1 || strcmp(tag, "P")) exit(1);
                    int c = 0, x; long pos;
                    /* read ints until next tag */
                    while (fscanf(f, "%d", &x) == 1) parts[i][c++] = x;
                    plen[i] = c;
                    /* the failed fscanf consumed nothing; continue */
                    (void)pos;
                }
                for (int i = 0; i < np; i++) {
                    if (fscanf(f, "%3s", tag) != 1 || strcmp(tag, "K")) exit(1);
                    for (int j = 0; j < np; j++)
                        if (fscanf(f, "%lld", &K[i][j]) != 1) exit(1);
                }
                fclose(f);
                return;
            }
        }
    }
    fprintf(stderr, "no kostka block for n=%d\n", n);
    exit(1);
}

static int parse_g6(const char *s) {
    int n = s[0] - 63;
    if (n != N) return -1;
    memset(adj, 0, sizeof adj);
    int bitpos = 0;
    const char *p = s + 1;
    for (int j = 1; j < n; j++)
        for (int i = 0; i < j; i++, bitpos++)
            if (((p[bitpos / 6] - 63) >> (5 - bitpos % 6)) & 1) {
                adj[i] |= (uint16_t)(1u << j);
                adj[j] |= (uint16_t)(1u << i);
            }
    return n;
}

static inline int pc(uint32_t x) { return __builtin_popcount(x); }

static int clawfree(void) {
    for (int v = 0; v < N; v++) {
        int nb[MAXN], d = 0;
        uint32_t m = adj[v];
        while (m) { nb[d++] = __builtin_ctz(m); m &= m - 1; }
        for (int a = 0; a < d; a++)
            for (int b = a + 1; b < d; b++) {
                if (adj[nb[a]] & (1u << nb[b])) continue;
                for (int c = b + 1; c < d; c++)
                    if (!(adj[nb[a]] & (1u << nb[c])) &&
                        !(adj[nb[b]] & (1u << nb[c]))) return 0;
            }
    }
    return 1;
}

/* independent sets grouped by size */
static uint16_t *ind_sets[MAXN + 1];
static int ind_cnt[MAXN + 1], ind_cap[MAXN + 1];

static void build_ind(void) {
    for (int s = 0; s <= N; s++) ind_cnt[s] = 0;
    for (uint32_t m = 1; m <= FULL; m++) {
        int ok = 1;
        uint32_t mm = m;
        while (mm) {
            int v = __builtin_ctz(mm); mm &= mm - 1;
            if (adj[v] & m) { ok = 0; break; }
        }
        if (!ok) continue;
        int s = pc(m);
        if (ind_cnt[s] == ind_cap[s]) {
            ind_cap[s] = ind_cap[s] ? 2 * ind_cap[s] : 256;
            ind_sets[s] = realloc(ind_sets[s], ind_cap[s] * sizeof(uint16_t));
        }
        ind_sets[s][ind_cnt[s]++] = (uint16_t)m;
    }
}

static long long f_dp[1 << MAXN], g_dp[1 << MAXN];

static long long ordered_count(int pi) {
    /* count ordered partitions with sizes parts[pi][0..plen-1] */
    memset(f_dp, 0, sizeof(long long) << N);
    f_dp[0] = 1;
    int used = 0;
    for (int step = 0; step < plen[pi]; step++) {
        int s = parts[pi][step];
        memset(g_dp, 0, sizeof(long long) << N);
        int found = 0;
        for (uint32_t m = 0; m <= FULL; m++) {
            if (!f_dp[m] || pc(m) != used) continue;
            for (int t = 0; t < ind_cnt[s]; t++) {
                uint32_t I = ind_sets[s][t];
                if (I & m) continue;
                g_dp[m | I] += f_dp[m];
                found = 1;
            }
        }
        if (!found) return 0;
        memcpy(f_dp, g_dp, sizeof(long long) << N);
        used += s;
    }
    return f_dp[FULL];
}

int main(int argc, char **argv) {
    if (argc < 2) { fprintf(stderr, "usage: clawsweep n [datafile]\n"); return 1; }
    N = atoi(argv[1]);
    FULL = (1u << N) - 1;
    load_kostka(argc > 2 ? argv[2] : "kostka_data.txt", N);

    char line[64];
    long long total = 0, cf = 0, witnesses = 0;
    long long a[MAXP], c[MAXP];
    while (fgets(line, sizeof line, stdin)) {
        size_t L = strlen(line);
        while (L && (line[L-1] == '\n' || line[L-1] == '\r')) line[--L] = 0;
        if (!L) continue;
        if (parse_g6(line) < 0) continue;
        total++;
        if (!clawfree()) continue;
        cf++;
        build_ind();
        int neg = 0;
        for (int i = 0; i < NP; i++) a[i] = ordered_count(i);
        /* back-substitute K^T c = a in lex-descending order */
        for (int lam = 0; lam < NP; lam++) {
            __int128 v = a[lam];
            for (int mu = 0; mu < lam; mu++)
                if (K[mu][lam] && c[mu]) v -= (__int128)K[mu][lam] * c[mu];
            c[lam] = (long long)v;
            if (c[lam] < 0) neg = 1;
        }
        if (neg) {
            witnesses++;
            printf("WITNESS %s negs:", line);
            for (int i = 0; i < NP; i++)
                if (c[i] < 0) {
                    printf(" [");
                    for (int j = 0; j < plen[i]; j++) printf("%d%s", parts[i][j], j+1<plen[i]?",":"");
                    printf("]=%lld", c[i]);
                }
            printf("\n");
            fflush(stdout);
        }
    }
    fprintf(stderr, "done n=%d: total=%lld clawfree=%lld witnesses=%lld\n",
            N, total, cf, witnesses);
    return 0;
}
