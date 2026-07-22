# Minimality of the 12-vertex counterexample

Every connected claw-free graph on at most 11 vertices has Schur-positive
chromatic symmetric function (all Schur coefficients computed exactly):

| order | connected graphs | claw-free | non-Schur-positive |
|---|---|---|---|
| ≤ 8 | 12,113 | 1,145 | 0 |
| 9  | 261,080 | 4,494 | 0 |
| 10 | 11,716,571 | 26,389 | 0 |
| 11 | 1,006,700,565 | 184,749 | 0 |
| **total** | **1,018,690,329** | **216,777** | **0** |

The previously quoted number 215,632 is the subtotal for orders 9 through
11; including orders 1 through 8 gives the correct total 216,777.

Tools: `clawsweep.c` (per-graph: claw filter, stable-partition DP for
monomial coefficients, exact integer Kostka back-substitution; Kostka data
from `kostka_data.txt`), self-tested on the counterexample itself, complete
graphs, and full order-8/9 runs. `verify_claw_schur.py` is a third
independent verifier for the counterexample coefficient with its own
self-tests. Disconnected graphs reduce to components (Schur positivity is
closed under products). Order-12 uniqueness computation in progress.
