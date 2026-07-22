#!/usr/bin/env python3
"""Independent exact verification by power sums and Frobenius characters.

This program deliberately does not import verify.py.  It enumerates all 2^22
edge subsets in Stanley's power-sum formula and obtains the target irreducible
character values from the Frobenius alternant formula.
"""

from __future__ import annotations

import argparse
import itertools
import math
from collections import Counter


EDGES = (
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


def signed_component_counts() -> Counter[tuple[int, ...]]:
    totals: Counter[tuple[int, ...]] = Counter()
    for mask in range(1 << len(EDGES)):
        parent = list(range(12))

        def find(vertex: int) -> int:
            while parent[vertex] != vertex:
                parent[vertex] = parent[parent[vertex]]
                vertex = parent[vertex]
            return vertex

        for index, (u, v) in enumerate(EDGES):
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

    power = signed_component_counts()
    target = (3, 3, 3, 3)
    contributions = []
    for cycle_type in sorted(power, reverse=True):
        character = frobenius_character(target, cycle_type)
        contributions.append((cycle_type, power[cycle_type], character, power[cycle_type] * character))

    coefficient = sum(row[3] for row in contributions)
    assert len(power) == 77  # every partition of 12 occurs
    assert frobenius_character(target, (1,) * 12) == 462
    character_norm_numerator = sum(
        frobenius_character(target, cycle_type) ** 2
        * (math.factorial(12) // centralizer_size(cycle_type))
        for cycle_type in partitions(12)
    )
    assert character_norm_numerator == math.factorial(12)
    chromatic_values = {
        colors: sum(signed * colors ** len(cycle_type) for cycle_type, signed in power.items())
        for colors in range(8)
    }
    assert chromatic_values == {
        0: 0, 1: 0, 2: 0, 3: 0, 4: 5376,
        5: 758160, 6: 23224320, 7: 325500000,
    }
    assert coefficient == -64

    if args.table:
        print("cycle_type\tsigned_subsets\tcharacter\tcontribution")
        for cycle_type, signed, character, contribution in contributions:
            print(f"{','.join(map(str, cycle_type))}\t{signed}\t{character}\t{contribution}")
    print("edge subsets: 4194304")
    print("power-sum component types: 77")
    print("character checks: dimension=462; row norm=1")
    print("chromatic polynomial check: P_G(4)=5376")
    print("[s_(3,3,3,3)] X_G = -64")
    print("PASS")


if __name__ == "__main__":
    main()
