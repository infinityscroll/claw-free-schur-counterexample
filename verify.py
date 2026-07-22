#!/usr/bin/env python3
"""Exact compact verification of the 12-vertex counterexample.

Only the Python standard library is required.  The calculation follows the
stable-partition/Kostka proof in README.md.
"""

from __future__ import annotations

import itertools
import math
from collections import Counter


ROOT_EDGES = (
    "ab", "bc", "cd", "ad",  # the four-cycle
    "au", "av", "uv",        # triangle at a
    "cx", "cy", "xy",        # triangle at c
    "bl", "dm",               # leaf edges
)

# This order gives the advertised graph6 string.
LINE_VERTICES = (
    "av", "bc", "xy", "dm", "au", "bl",
    "uv", "cy", "ad", "cx", "ab", "cd",
)


def line_graph_edges() -> tuple[tuple[int, int], ...]:
    result = []
    for i, left in enumerate(LINE_VERTICES):
        for j, right in enumerate(LINE_VERTICES[i + 1 :], i + 1):
            if set(left) & set(right):
                result.append((i, j))
    return tuple(result)


def graph6(n: int, edges: tuple[tuple[int, int], ...]) -> str:
    """Encode a graph of order at most 62 using the short graph6 format."""
    edge_set = {tuple(sorted(edge)) for edge in edges}
    bits = [int((i, j) in edge_set) for j in range(1, n) for i in range(j)]
    bits.extend([0] * (-len(bits) % 6))
    payload = []
    for start in range(0, len(bits), 6):
        value = 0
        for bit in bits[start : start + 6]:
            value = 2 * value + bit
        payload.append(chr(value + 63))
    return chr(n + 63) + "".join(payload)


def adjacency(n: int, edges: tuple[tuple[int, int], ...]) -> list[set[int]]:
    result = [set() for _ in range(n)]
    for u, v in edges:
        result[u].add(v)
        result[v].add(u)
    return result


def is_connected(adj: list[set[int]]) -> bool:
    reached = {0}
    todo = [0]
    while todo:
        vertex = todo.pop()
        for neighbor in adj[vertex] - reached:
            reached.add(neighbor)
            todo.append(neighbor)
    return len(reached) == len(adj)


def claws(adj: list[set[int]]) -> list[tuple[int, int, int, int]]:
    result = []
    for center, neighbors in enumerate(adj):
        for leaves in itertools.combinations(sorted(neighbors), 3):
            if all(v not in adj[u] for u, v in itertools.combinations(leaves, 2)):
                result.append((center, *leaves))
    return result


def independence_number(adj: list[set[int]]) -> int:
    for size in range(len(adj), -1, -1):
        for subset in itertools.combinations(range(len(adj)), size):
            if all(v not in adj[u] for u, v in itertools.combinations(subset, 2)):
                return size
    raise AssertionError("unreachable")


def normalized_edge_colorings() -> Counter[tuple[int, ...]]:
    """Count stable partitions by normalizing the four colors at vertex a."""
    colors = set(range(4))
    totals: Counter[tuple[int, ...]] = Counter()
    admissible = 0
    for bc, cd, cx, cy in itertools.permutations(range(4)):
        if bc == 0 or cd == 1:
            continue
        admissible += 1
        choices = (
            colors - {0, bc},      # bl
            colors - {1, cd},      # dm
            colors - {2, 3},       # uv
            colors - {cx, cy},     # xy
        )
        for bl, dm, uv, xy in itertools.product(*(sorted(s) for s in choices)):
            coloring = {
                "ab": 0, "ad": 1, "au": 2, "av": 3,
                "bc": bc, "cd": cd, "cx": cx, "cy": cy,
                "bl": bl, "dm": dm, "uv": uv, "xy": xy,
            }
            # Verify properness directly, independently of the parametrization.
            for left, right in itertools.combinations(ROOT_EDGES, 2):
                if set(left) & set(right):
                    assert coloring[left] != coloring[right]
            sizes = tuple(sorted(Counter(coloring.values()).values(), reverse=True))
            totals[sizes] += 1
    assert admissible == 14
    return totals


def kostka(shape: tuple[int, ...], weight: tuple[int, ...]) -> int:
    """Enumerate SSYT of the given shape and content."""
    cells = [(row, col) for row, width in enumerate(shape) for col in range(width)]
    table = [[-1] * width for width in shape]
    remaining = list(weight)

    def visit(position: int) -> int:
        if position == len(cells):
            return 1
        row, col = cells[position]
        lower = 0
        if col:
            lower = max(lower, table[row][col - 1])
        if row and col < len(table[row - 1]):
            lower = max(lower, table[row - 1][col] + 1)
        total = 0
        for value in range(lower, len(remaining)):
            if not remaining[value]:
                continue
            remaining[value] -= 1
            table[row][col] = value
            total += visit(position + 1)
            remaining[value] += 1
        table[row][col] = -1
        return total

    return visit(0)


def main() -> None:
    edges = line_graph_edges()
    adj = adjacency(12, edges)
    assert len(edges) == 22
    assert graph6(12, edges) == "K?`CRAWWUXIM"
    assert is_connected(adj)
    assert claws(adj) == []
    assert independence_number(adj) == 4

    stable = normalized_edge_colorings()
    expected = Counter({(4, 4, 2, 2): 32, (4, 3, 3, 2): 160, (3, 3, 3, 3): 32})
    assert stable == expected

    monomial_4422 = stable[(4, 4, 2, 2)] * math.factorial(2) ** 2
    monomial_4332 = stable[(4, 3, 3, 2)] * math.factorial(2)
    monomial_3333 = stable[(3, 3, 3, 3)] * math.factorial(4)
    assert (monomial_4422, monomial_4332, monomial_3333) == (128, 320, 768)

    k_4422_4332 = kostka((4, 4, 2, 2), (4, 3, 3, 2))
    k_4422_3333 = kostka((4, 4, 2, 2), (3, 3, 3, 3))
    k_4332_3333 = kostka((4, 3, 3, 2), (3, 3, 3, 3))
    assert (k_4422_4332, k_4422_3333, k_4332_3333) == (1, 2, 3)

    schur_4422 = monomial_4422
    schur_4332 = monomial_4332 - k_4422_4332 * schur_4422
    schur_3333 = (
        monomial_3333
        - k_4422_3333 * schur_4422
        - k_4332_3333 * schur_4332
    )
    assert (schur_4422, schur_4332, schur_3333) == (128, 192, -64)

    print("graph6: K?`CRAWWUXIM")
    print("vertices: 12; edges: 22; connected: yes; claw-free: yes")
    print("stable partitions: 4422=32, 4332=160, 3333=32")
    print("monomial coefficients: m_4422=128, m_4332=320, m_3333=768")
    print("Kostka numbers: K_4422,4332=1, K_4422,3333=2, K_4332,3333=3")
    print("Schur coefficients: s_4422=128, s_4332=192, s_3333=-64")
    print("PASS")


if __name__ == "__main__":
    main()
