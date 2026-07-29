# Provenance and chronology

The counterexample was first made public in this repository. The
authoritative timestamps are recorded by third parties and can be verified
independently:

| Event | UTC time | Verifiable via |
|---|---|---|
| Initial commit `7747970833d8d4ce2dadba8f89a184cb52d814a9` ("Give a counterexample to claw-free Schur positivity") | 2026-07-22 00:39:51 | `git log`; tagged `v1.0.0` |
| Repository created (GitHub server-side) | 2026-07-22 00:39:59 | `GET /repos/infinityscroll/claw-free-schur-counterexample` |
| Release v1.0.0 published (GitHub server-side) | 2026-07-22 00:40:15 | `GET /repos/infinityscroll/claw-free-schur-counterexample/releases` |
| MathOverflow question citing this repository | 2026-07-23 15:50 | the thread itself |
| Matherne-Morales arXiv:2607.21508 v1 submitted (independent work) | 2026-07-23 16:54:51 | arXiv submission history |

The `v1.0.0` tag still points to the initial commit, which is an ancestor of
`main`: the original history has never been rewritten. GitHub's `created_at`
and `published_at` fields are set server-side and cannot be backdated by a
client. Matherne and Morales have publicly confirmed that their work is
independent; the two pairs of counterexamples are isomorphic, as they must
be, since exactly two counterexamples of minimum order exist (see
`minimality/RESULTS.md`).
