# A 12-vertex counterexample to claw-free Schur positivity

This repository gives a connected claw-free graph $G$ whose chromatic
symmetric function is not Schur-positive.  In fact,


$$[s_{(3,3,3,3)}]X_G=-64.$$


Thus $G$ is a counterexample to Stanley's claw-free conjecture (credited to
Gasharov) that every claw-free graph has a Schur-positive chromatic symmetric
function.

## The graph

Start with the cycle $a b c d a$.  Attach a triangle $a u v a$ at
$a$, a triangle $c x y c$ at $c$, a leaf $\ell$ at $b$, and a
leaf $m$ at $d$.  Let $H$ be the resulting graph and set
$G=L(H)$, its line graph.

The twelve vertices of $G$, in the order

```text
av, bc, xy, dm, au, bl, uv, cy, ad, cx, ab, cd
```

have graph6 encoding

```text
K?`CRAWWUXIM
```

and edge set

```text
04 06 08 0A 15 17 19 1A 1B 27 29 38 3B 46 48 4A 5A 79 7B 8A 8B 9B
```

where `A` and `B` denote vertices 10 and 11.

Every line graph is claw-free: three edges of $H$ meeting a fixed edge
cannot be pairwise disjoint, since two of them meet the same endpoint.  The
graph is connected because $H$ is connected.

## Coefficient calculation

A stable set of $G=L(H)$ is a matching of $H$.  The matching number of
$H$ is four.  Indeed, $\{uv,xy,b\ell,dm\}$ is a four-edge matching.  A
perfect matching would have to contain both leaf edges $b\ell$ and $dm$;
the two remaining attached triangles would then each have three unmatched
vertices, which is impossible.

It remains to count partitions of $E(H)$ into four matchings.  At the
degree-four vertex $a$, normalize the colors on
$(ab,ad,au,av)$ to $(1,2,3,4)$.  The colors on
$(bc,cd,cx,cy)$ form a permutation $p$ with $p_1\ne1$ and
$p_2\ne2$, leaving fourteen possibilities.  For every such $p$, each of
$b\ell,dm,uv,xy$ has two choices.  Directly sorting the resulting color
class sizes gives

| permutations at $c$ | $(4,4,2,2)$ | $(4,3,3,2)$ | $(3,3,3,3)$ |
|---|---:|---:|---:|
| `2134`, `2143` | 4 | 8 | 4 |
| the other 12 admissible permutations | 2 each | 12 each | 2 each |
| total | 32 | 160 | 32 |

Normalization at $a$ gives one representative of every unlabeled stable
partition.  Hence, using


$$[m_\lambda]X_G=N_\lambda\prod_j m_j(\lambda)!,$$


the relevant monomial coefficients are


$$[m_{4422}]X_G=32(2!)(2!)=128,\qquad
 [m_{4332}]X_G=160(2!)=320,\qquad
 [m_{3333}]X_G=32(4!)=768.$$


Since the largest stable set has size four, dominance-unitriangularity of the
Kostka matrix makes every Schur coefficient indexed by a partition with first
part greater than four vanish.  Among the remaining partitions that dominate
$(3,3,3,3)$, the only possibilities are
$(4,4,4)$, $(4,4,3,1)$, $(4,4,2,2)$, $(4,3,3,2)$, and
$(3,3,3,3)$.  There are no colorings of the first two types.  The only
Kostka numbers needed below are


$$K_{4422,4332}=1,\qquad K_{4422,3333}=2,\qquad
 K_{4332,3333}=3.$$


Since $s_\nu=\sum_\lambda K_{\nu\lambda}m_\lambda$, triangular inversion
now gives


$$[s_{4422}]X_G=128,$$



$$[s_{4332}]X_G=320-128=192,$$


and finally


$$[s_{3333}]X_G=768-2(128)-3(192)=-64.$$


## Verification

The short verifier uses only the Python standard library.  It constructs both
graphs, checks the graph6 record and claw-freeness, enumerates the normalized
edge colorings, enumerates the relevant semistandard Young tableaux, and
repeats the triangular inversion.

```bash
python3 verify.py
```

An independent verifier uses Stanley's power-sum inclusion--exclusion formula


$$X_G=\sum_{A\subseteq E(G)}(-1)^{|A|}p_{\lambda(V,A)}$$


together with the Frobenius character formula.  It enumerates all $2^{22}$
edge subsets and again obtains $-64$.

```bash
python3 verify_power_sum.py
```

## Minimum-order classification

An exhaustive `geng` census gives the following exact result.

- Every claw-free graph on at most 11 vertices is Schur-positive.
- Among the 1,728,404 connected claw-free graphs on 12 vertices, exactly two
  isomorphism classes are not Schur-positive.

Their graph6 records and negative coefficients are

```text
K?`CRAWWUXIM  [s_(3,3,3,3)] = -64
K?`CR@`bAbRB  [s_(3,3,3,3)] = -40
```

The second graph is also a line graph. One root graph is obtained from a
5-cycle by attaching triangles at two distance-two cycle vertices and a
pendant edge at the intervening cycle vertex. The complete shard totals and
witness records are in [`minimality/RESULTS.md`](minimality/RESULTS.md).
Disconnected graphs of order at most 12 reduce to smaller connected
components, since Schur positivity is closed under products. Thus these are
exactly the two counterexamples of minimum order.

## References

- R. P. Stanley, *A symmetric function generalization of the chromatic
  polynomial of a graph*, Advances in Mathematics 111 (1995), 166--194.
- V. Gasharov, *On Stanley's chromatic symmetric function and clawfree
  graphs*, Discrete Mathematics 205 (1999), 229--234.
- E. Shelburne and S. van Willigenburg, *Schur-positivity for generalized
  nets*, Enumerative Combinatorics and Applications 5:2 (2025), Article
  S2R14, [arXiv:2409.00943](https://arxiv.org/abs/2409.00943).

## License

The verification code is released under the MIT License.  The mathematical
text is released under CC BY 4.0.
