#!/usr/bin/env python3
"""Fully independent verification of the claimed claw-free Schur-positivity
counterexample: G = line graph of (4-cycle + triangles at two opposite
vertices + pendants at the other two); claim [s_3333] X_G = -64.

Method (no external math libraries):
  1. a_lam = [m_lam] X_G by DP over vertex subsets: number of ordered
     partitions of V into independent sets with sizes lam_1, lam_2, ...
  2. Kostka matrix K_{mu,lam} = #SSYT(shape mu, content lam) via the
     horizontal-strip recursion.
  3. Solve K^T c = a exactly (Fractions) to get Schur coefficients c_mu.
Self-tests: X_{K_n} = n! * s_{1^n}; claw K_{1,3} must be non-Schur-positive;
m-coefficients cross-checked against brute-force colorings on 6 vertices.
"""
from fractions import Fraction
from functools import lru_cache
from itertools import combinations, permutations, product


# ---------- partitions ----------
def partitions(n, maxpart=None):
    if maxpart is None:
        maxpart = n
    if n == 0:
        yield ()
        return
    for p in range(min(n, maxpart), 0, -1):
        for rest in partitions(n - p, p):
            yield (p,) + rest


# ---------- chromatic symmetric function in the m-basis ----------
def csf_m_coeffs(n, edges):
    adjmask = [0] * n
    for a, b in edges:
        adjmask[a] |= 1 << b
        adjmask[b] |= 1 << a
    full = (1 << n) - 1
    # independent sets grouped by size
    ind_by_size = {}
    for mask in range(1, full + 1):
        ok = True
        mm = mask
        while mm:
            v = (mm & -mm).bit_length() - 1
            mm &= mm - 1
            if adjmask[v] & mask:
                ok = False
                break
        if ok:
            ind_by_size.setdefault(mask.bit_count(), []).append(mask)

    def ordered_count(lam):
        # DP: f[used_mask] over prefix of lam
        f = {0: 1}
        for size in lam:
            g = {}
            for used, cnt in f.items():
                for I in ind_by_size.get(size, []):
                    if I & used:
                        continue
                    g[used | I] = g.get(used | I, 0) + cnt
            f = g
            if not f:
                return 0
        return f.get(full, 0)

    return {lam: ordered_count(lam) for lam in partitions(n)}


# ---------- Kostka numbers ----------
@lru_cache(maxsize=None)
def kostka(shape, content):
    # SSYT of given shape and content, peeling the largest letter as a
    # horizontal strip from the end of the content.
    if sum(shape) != sum(content):
        return 0
    if not content:
        return 1 if not shape else 0
    k = content[-1]
    total = 0
    for smaller in horizontal_strip_predecessors(shape, k):
        total += kostka(smaller, content[:-1])
    return total


def horizontal_strip_predecessors(shape, k):
    """All partitions nu with shape/nu a horizontal strip of size k."""
    rows = len(shape)
    def rec(i, remaining, acc):
        if i == rows:
            if remaining == 0:
                nu = tuple(x for x in acc if x > 0)
                if all(nu[j] >= nu[j + 1] for j in range(len(nu) - 1)):
                    yield nu
            return
        lower = shape[i + 1] if i + 1 < rows else 0
        # nu_i between lower..shape_i, and horizontal strip: nu_i >= shape_{i+1}
        for nu_i in range(shape[i], lower - 1, -1):
            take = shape[i] - nu_i
            if take > remaining:
                continue
            # strip condition: nu_i >= shape[i+1] (cells removed in row i must
            # not sit above cells of row i+1) -- guaranteed since nu_i >= lower
            yield from rec(i + 1, remaining - take, acc + [nu_i])
    yield from rec(0, k, [])


def schur_coeffs(n, m_coeffs):
    plist = sorted(partitions(n))
    idx = {lam: i for i, lam in enumerate(plist)}
    N = len(plist)
    # K[mu][lam]
    K = [[Fraction(kostka(mu, lam)) for lam in plist] for mu in plist]
    a = [Fraction(m_coeffs.get(lam, 0)) for lam in plist]
    # solve K^T c = a
    M = [[K[mu][lam] for mu in range(N)] for lam in range(N)]  # rows lam, cols mu
    # gaussian elimination
    c = a[:]
    Mat = [row[:] + [c[i]] for i, row in enumerate(M)]
    for col in range(N):
        piv = next(r for r in range(col, N) if Mat[r][col] != 0)
        Mat[col], Mat[piv] = Mat[piv], Mat[col]
        pv = Mat[col][col]
        Mat[col] = [x / pv for x in Mat[col]]
        for r in range(N):
            if r != col and Mat[r][col] != 0:
                f = Mat[r][col]
                Mat[r] = [x - f * y for x, y in zip(Mat[r], Mat[col])]
    return {plist[mu]: Mat[mu][N] for mu in range(N)}


# ---------- brute force m-coeff cross-check ----------
def brute_m_coeffs(n, edges, lam):
    # count proper colorings with color class sizes exactly lam (ordered)
    k = len(lam)
    cnt = 0
    for coloring in product(range(k), repeat=n):
        sizes = [0] * k
        for c in coloring:
            sizes[c] += 1
        if tuple(sizes) != lam:
            continue
        if all(coloring[a] != coloring[b] for a, b in edges):
            cnt += 1
    return cnt


def run_selftests():
    # K_n -> n! s_{1^n}
    for n in (3, 4, 5):
        edges = list(combinations(range(n), 2))
        c = schur_coeffs(n, csf_m_coeffs(n, edges))
        expect = {tuple([1] * n): Fraction(__import__('math').factorial(n))}
        got = {k: v for k, v in c.items() if v != 0}
        assert got == expect, (n, got)
    print("selftest K_n OK (n=3,4,5)")
    # claw K_{1,3}: known to be non-Schur-positive
    claw = [(0, 1), (0, 2), (0, 3)]
    c = schur_coeffs(4, csf_m_coeffs(4, claw))
    neg = {k: v for k, v in c.items() if v < 0}
    assert neg, "claw should have a negative Schur coefficient"
    print("selftest claw non-s-positive OK:", dict(neg))
    # brute-force m-coeff cross-check on a 6-vertex graph
    edges6 = [(0,1),(1,2),(2,3),(3,4),(4,5),(5,0),(0,3)]
    mc = csf_m_coeffs(6, edges6)
    for lam in [(2,2,2), (3,2,1), (2,2,1,1), (4,2)]:
        assert mc[lam] == brute_m_coeffs(6, edges6, lam), lam
    print("selftest brute-force m-coeffs OK")


def decode_graph6(code):
    """Decode a short-form graph6 record without external libraries."""
    n = ord(code[0]) - 63
    assert 0 <= n <= 62
    bits = []
    for character in code[1:]:
        value = ord(character) - 63
        assert 0 <= value < 64
        bits.extend((value >> shift) & 1 for shift in range(5, -1, -1))
    pairs = [(i, j) for j in range(1, n) for i in range(j)]
    assert len(bits) >= len(pairs)
    return n, [pair for pair, bit in zip(pairs, bits) if bit]


def check_connected_claw_free(n, edges):
    adj = [set() for _ in range(n)]
    for a, b in edges:
        adj[a].add(b)
        adj[b].add(a)
    reached = {0}
    frontier = [0]
    while frontier:
        vertex = frontier.pop()
        for neighbor in adj[vertex] - reached:
            reached.add(neighbor)
            frontier.append(neighbor)
    assert len(reached) == n
    for vertex in range(n):
        for trio in combinations(adj[vertex], 3):
            assert any(y in adj[x] for x, y in combinations(trio, 2))


if __name__ == "__main__":
    run_selftests()

    # Build root graph H from the PR's prose description:
    # 4-cycle a0 a1 a2 a3; triangles attached at a0 and a2 (opposite);
    # pendant edges at a1 and a3.
    H_edges = [(0,1),(1,2),(2,3),(3,0),          # 4-cycle
               (0,4),(0,5),(4,5),                # triangle at 0
               (2,6),(2,7),(6,7),                # triangle at 2
               (1,8),                            # pendant at 1
               (3,9)]                            # pendant at 3
    assert len(H_edges) == 12
    # line graph
    E = H_edges
    n = len(E)
    G_edges = [(i, j) for i in range(n) for j in range(i + 1, n)
               if set(E[i]) & set(E[j])]
    print(f"L(H): {n} vertices, {len(G_edges)} edges")
    # Line graphs always are claw-free; verify connectivity and claw-freeness
    # directly so this check does not rely on the construction theorem.
    check_connected_claw_free(n, G_edges)
    print("claw-free: verified directly")

    mc = csf_m_coeffs(n, G_edges)
    c = schur_coeffs(n, mc)
    target = (3, 3, 3, 3)
    assert c[target] == -64
    print(f"[s_{target}] X_G = {c[target]}")
    negs = {k: v for k, v in c.items() if v < 0}
    assert negs == {target: Fraction(-64)}
    print("all negative Schur coefficients:", negs)

    # The exhaustive order-12 census found one further isomorphism class.
    # Decode it independently of the census program and recompute all 77
    # Schur coefficients from scratch.
    second_code = "K?`CR@`bAbRB"
    n2, second_edges = decode_graph6(second_code)
    assert n2 == 12 and len(second_edges) == 21
    check_connected_claw_free(n2, second_edges)
    second_coefficients = schur_coeffs(n2, csf_m_coeffs(n2, second_edges))
    second_negs = {
        shape: value
        for shape, value in second_coefficients.items()
        if value < 0
    }
    assert second_negs == {target: Fraction(-40)}
    print(f"graph6 {second_code}: 12 vertices, 21 edges; connected and claw-free")
    print("all negative Schur coefficients:", second_negs)
    print("PASS")
