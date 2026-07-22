# Minimum-order classification

Every connected claw-free graph on at most 11 vertices has Schur-positive
chromatic symmetric function (all Schur coefficients computed exactly):

| order | connected graphs | claw-free | non-Schur-positive |
|---|---|---|---|
| ≤ 8 | 12,113 | 1,145 | 0 |
| 9  | 261,080 | 4,494 | 0 |
| 10 | 11,716,571 | 26,389 | 0 |
| 11 | 1,006,700,565 | 184,749 | 0 |
| 12 | 164,059,830,476 | 1,728,404 | 2 |
| **total through 12** | **165,078,520,805** | **1,945,181** | **2** |

The previously quoted number 215,632 is the subtotal for orders 9 through
11; including orders 1 through 8 gives the correct total 216,777.

The two order-12 isomorphism classes have graph6 records
``K?`CRAWWUXIM`` and ``K?`CR@`bAbRB``, with respective negative
coefficients `[s_(3,3,3,3)] = -64` and `-40`. The twelve order-12 shard logs
sum to the known number 164,059,830,476 of connected unlabeled graphs, so
this is an exhaustive minimum-order classification, not a sample search.

Tools: `clawsweep.c` (per-graph: claw filter, stable-partition DP for
monomial coefficients, exact integer Kostka back-substitution; Kostka data
from `kostka_data.txt`), self-tested on the counterexample itself, complete
graphs, and full order-8/9 runs. `verify_claw_schur.py` is a third
independent verifier for both counterexamples and every Schur coefficient,
with its own self-tests. `verify_power_sum.py` independently checks both
negative coefficients by Stanley's edge-subset formula. Disconnected graphs
of order at most 12 reduce to smaller connected components, and Schur
positivity is closed under products.
