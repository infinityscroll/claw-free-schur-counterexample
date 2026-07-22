#!/usr/bin/env python3
"""Independent exact verification by power sums and Frobenius characters.

This program deliberately does not import either graph verifier.  It
enumerates all edge subsets in Stanley's power-sum formula for both
minimum-order counterexamples and obtains the target irreducible character
values from the Frobenius alternant formula.
"""

from __future__ import annotations

import argparse
import itertools
import math
from collections import Counter


FIRST_EDGES = (
    (0, 4), (0, 6), (0, 8), (0, 10),
    (1, 5), (1, 7), (1, 9), (1, 10), (1, 11),
    (2, 7), (2, 9),
    (3, 8), (3, 11),
    (4, 6), (4, 8), (4, 10),
    (5, 10),
    (7, 9), (7, 11),
    (8, 10), (8, 11),
    (9, 11),
)

SECOND_EDGES = (
    (0, 4), (0, 6), (0, 9), (0, 11),
    (1, 5), (1, 7), (1, 8), (1, 10),
    (2, 7), (2, 8),
    (3, 10), (3, 11),
    (4, 6), (4, 9), (4, 11),
    (5, 9),
    (7, 8), (7, 10),
    (8, 10),
    (9, 11),
    (10, 11),
)


def permutation_sign(permutation: tuple[int, ...]) -> int:
    inversions = sum(
        permutation[i] > permutation[j]
        for i in range(len(permutation))
        for j in range(i + 1, len(permutation))
    )
    return -1 if inversions % 2 else 1


def assignments_with_sums(parts: tuple[int, ...], targets: tuple[int, ...]) -> int:
    """Coefficient of x^targets in product_j sum_i x_i^parts[j]."""
    states: Counter[tuple[int, ...]] = Counter({(0,) * len(targets): 1})
    for part in parts:
        following: Counter[tuple[int, ...]] = Counter()
        for sums, multiplicity in states.items():
            for index in range(len(targets)):
                if sums[index] + part <= targets[index]:
                    new_sums = list(sums)
                    new_sums[index] += part
                    following[tuple(new_sums)] += multiplicity
        states = following
    return states[targets]


def frobenius_character(shape: tuple[int, ...], cycle_type: tuple[int, ...]) -> int:
    """Compute chi^shape(cycle_type) from the Frobenius alternant formula."""
    rows = len(shape)
    target = tuple(shape[i] + rows - 1 - i for i in range(rows))
    total = 0
    for permutation in itertools.permutations(range(rows)):
        alternant = tuple(rows - 1 - permutation[i] for i in range(rows))
        remainder = tuple(target[i] - alternant[i] for i in range(rows))
        if min(remainder) < 0:
            continue
        total += permutation_sign(permutation) * assignments_with_sums(cycle_type, remainder)
    return total


def partitions(total: int, upper: int | None = None):
    if total == 0:
        yield ()
        return
    upper = total if upper is None else min(upper, total)
    for first in range(upper, 0, -1):
        for rest in partitions(total - first, first):
            yield (first,) + rest


def centralizer_size(cycle_type: tuple[int, ...]) -> int:
    multiplicities = Counter(cycle_type)
    return math.prod(
        part ** count * math.factorial(count)
        for part, count in multiplicities.items()
    )


def signed_component_counts(
    edges: tuple[tuple[int, int], ...],
) -> Counter[tuple[int, ...]]:
    totals: Counter[tuple[int, ...]] = Counter()
    for mask in range(1 << len(edges)):
        parent = list(range(12))

        def find(vertex: int) -> int:
            while parent[vertex] != vertex:
                parent[vertex] = parent[parent[vertex]]
                vertex = parent[vertex]
            return vertex

        for index, (u, v) in enumerate(edges):
            if not (mask >> index) & 1:
                continue
            ru, rv = find(u), find(v)
            if ru != rv:
                parent[rv] = ru
        sizes = Counter(find(vertex) for vertex in range(12))
        cycle_type = tuple(sorted(sizes.values(), reverse=True))
        totals[cycle_type] += -1 if mask.bit_count() % 2 else 1
    return totals


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--table", action="store_true", help="print every power-sum contribution")
    args = parser.parse_args()

    target = (3, 3, 3, 3)
    assert frobenius_character(target, (1,) * 12) == 462
    character_norm_numerator = sum(
        frobenius_character(target, cycle_type) ** 2
        * (math.factorial(12) // centralizer_size(cycle_type))
        for cycle_type in partitions(12)
    )
    assert character_norm_numerator == math.factorial(12)
    cases = (
        ("K?`CRAWWUXIM", FIRST_EDGES, -64, 5376),
        ("K?`CR@`bAbRB", SECOND_EDGES, -40, 7680),
    )
    for graph6, edges, expected, expected_chromatic_at_four in cases:
        power = signed_component_counts(edges)
        contributions = []
        for cycle_type in sorted(power, reverse=True):
            character = frobenius_character(target, cycle_type)
            contributions.append(
                (
                    cycle_type,
                    power[cycle_type],
                    character,
                    power[cycle_type] * character,
                )
            )
        coefficient = sum(row[3] for row in contributions)
        assert len(power) == 77  # every partition of 12 occurs
        assert coefficient == expected
        chromatic_at_four = sum(
            signed * 4 ** len(cycle_type)
            for cycle_type, signed in power.items()
        )
        assert chromatic_at_four == expected_chromatic_at_four
        if args.table:
            print(f"graph6={graph6}")
            print("cycle_type\tsigned_subsets\tcharacter\tcontribution")
            for cycle_type, signed, character, contribution in contributions:
                print(
                    f"{','.join(map(str, cycle_type))}\t{signed}\t"
                    f"{character}\t{contribution}"
                )
        print(
            f"graph6={graph6}; edges={len(edges)}; "
            f"edge_subsets={1 << len(edges)}; component_types={len(power)}; "
            f"P_G(4)={chromatic_at_four}; [s_(3,3,3,3)]={coefficient}"
        )
    print("character checks: dimension=462; row norm=1")
    print("PASS")


if __name__ == "__main__":
    main()
