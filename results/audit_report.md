# Audit of the first session's code and results (second session, 2026-10-12)

Scope: every script, test, log and result file left by the first session (and the original
archive in `threebody_codex_project/`) was read; the numerical claims were re-derived with
code written in this session that shares no integration, shooting or Floquet code with the
first session. Verdicts: **confirmed**, **confirmed with caveat**, **incorrect/misleading**,
or **gap**.

## 1. Reproducibility and provenance

| Item | Verdict | Evidence |
|---|---|---|
| Test suite (15 tests) | confirmed | re-run: 15 passed (`pytest tests`, before any change) |
| Original archive unchanged | confirmed | all 24 SHA256 entries of `threebody_codex_project/SHA256SUMS.txt` match |
| Raw MP files referenced by JSON records | confirmed | 106/106 raw CSVs exist and match their recorded SHA256 |
| Source hashes recorded in MP JSONs | confirmed with caveat | recorded hashes are of *earlier* versions of `taylor_mpfr.cpp` and `high_precision.py`; diffs against `results/source_snapshots/` show only input validation, a non-finite check, `--seed-file/--target` options and a later spectrum refinement — integration numerics unchanged. Snapshots for every version used are now present |
| Environment | confirmed with caveat | the offline environment lacked Numba/pandas/mpmath; this session installed `mpmath 1.4.1`, `gmpy2 2.3.2`, `numba 0.68.0` (pinned in `requirements*.txt`) without changing numpy/scipy |

## 2. Binary64 work (Phase A)

| Claim | Verdict | Evidence |
|---|---|---|
| 23-point sweep "reproduced exactly" | confirmed with caveat | true, but each archived point was its own seed, so Newton exited immediately and the same integration was repeated. `dense_sweep_double.py` re-derives all 23 points by continuation from the single mu = 0.613 seed: max \|du\|,\|dv\| = 4.0e-13, \|dT\| = 1.6e-11, \|dD\| = 3.0e-8, \|d rho\| = 6.4e-9 |
| "4 multipliers farthest from 1" selects the nontrivial quartet | confirmed | explicit symplectic-quotient spectrum vs heuristic on 110 masses in [0.580, 0.6345]: max difference 5.7e-8; nontrivial multipliers stay >= 1.356 from 1 |
| Failed archived convergence check diagnosed and replaced | confirmed | `pancake_convergence_check.log` (24 integrations, consistent signs); the exact historical trigger remains unknown (as the first session stated) |
| Binary64 stability window ~[0.59904, 0.61784] | confirmed, but **incomplete** | the branch has further stable stretches outside the scanned [0.585, 0.630] (section 5) |

## 3. Multiprecision work (Phases B/C)

| Claim | Verdict | Evidence |
|---|---|---|
| Seed closure ~5.14e-40 at mu = 0.613 | confirmed | independent GBS integrator reproduces all 12 return components to <= 5.5e-49; the value is a property of the 47-digit seed, not of the integrator. The seed was never Newton-corrected; corrected now (closure 1.4e-62), the parameters differ from the seed by 4.6e-41 (u), 2.8e-41 (v), 1.3e-39 (T) |
| Closures < 1e-28 for all 24 records | confirmed | GBS re-solution of all 12 masses: closures 6e-63 ... 3e-62; parameter agreement with the 85-digit Taylor boundary records 7e-65 (u), 6e-63 (T) |
| D, quartet, monodromy | confirmed | GBS + reduction-free trace formula vs Taylor + Decimal quotient: \|dD\| = 5e-59 (bracket endpoints, 85 digits), \|dPhi\| = 1e-60; all Krein signs agree |
| "Largest D change between precision runs ~8.3e-46" as convergence evidence | misleading | both runs of the 8-point sets stopped at the same Newton iterate, so the comparison cannot see the parameter error; the independent check shows the actual D error is up to ~8e-43 (7e-38 at mu = 0.613, uncorrected seed). Irrelevant for the signs |
| "Independently initialized, time-reversed seed also closes" (defect 1.25e-37) | **misleading** | (q0, -v0) equals the forward initial state rotated by pi with bodies 1 and 3 relabelled. The reversed run is a bitwise mirror: identical 1900-step time grid, identical closure to all 66 printed digits, every node equal to the mirrored forward node (difference 0). It tests the reversibility identity of the computed monodromy, not truncation error |
| Brackets [0.59904378, 0.59904379], [0.61784426, 0.61784428] | confirmed | refined to mu_L = 0.59904378420889985122072886046206841927541372566178042..., mu_R = 0.61784426820906317872736671060041799670207268916883745867616... (four pipelines; see `stability_transitions.md`) |
| Opposite Krein signatures at both brackets | confirmed | independent full-space Krein forms |
| "Exact-transition multiplicities unresolved" | resolved now | sigma_min(Phi - lambda_bar I) = 2.9e-61 / 2.5e-62 vs sigma_2 = 0.020 / 0.092 at mu*: geometric multiplicity 1, algebraic 2; splitting ~ delta^0.5000, rho - 1 ~ delta^0.5000, Krein forms ~ delta^0.501 over delta = 1e-4 ... 1e-40 |

## 4. Literature (Phase D)

| Claim | Verdict | Evidence |
|---|---|---|
| Draft PDF could not be fetched; audit incomplete | resolved now | draft, catalogue and orbit files downloaded; the draft asserts linear stability of one representative only and reports no transitions (`literature_audit.md`) |
| Independent external check | new | our values agree with the draft's 42-digit Table 1 to 3e-38 / 2e-38 / 1.5e-36 and its 19-digit multipliers to 4e-20; the catalogue `stab` values 1.1637 / 1.0 / 2.3815 are reproduced |

## 5. Scientific gaps found during the audit

1. **More than one stable stretch.** Pseudo-arclength continuation of the whole
   AABAbbaabaBB segment (`isola_scan_double.py`) finds three spectrally stable stretches,
   separated by four Krein collisions, plus a tangent (s = +2) crossing at the fold
   mu_max ~ 0.6346043 (matching the draft's 0.63460425). A short stable stretch also exists on
   the (ABab)^3 segment past its mu_min fold.
2. **Half-period reflection symmetry** (implicit in the atlas's "Euler-half-twist" sector).
   Every orbit satisfies X(T/2) = sigma X(0) for a reflection sigma (defect <= 8e-63 at 80
   digits), so M = N^2 with N = sigma DPhi_{T/2}. Hence
   s_M = s_N^2 - 2 >= -2: at the apparent "period-doubling" touches of the multiplier -1 near
   mu ~ 0.5695 and 0.5481 an N-pair passes +-i (a 1:4 resonance of N) and the spectrum of M
   stays on the unit circle -- no change of spectral stability. Nonlinear bifurcations there
   (e.g. symmetry-breaking period-2T orbits) are not excluded and were not studied. The
   near-collision near mu ~ 0.5619 is a pass-through of equal-Krein-sign pairs.
3. **Precision claims.** The first session's text treated the 85-digit/65-digit comparison as
   convergence evidence and the time-reversal run as independent; both are corrected above.

## 6. Code review notes (no numerical errors found)

* `tools/taylor_mpfr.cpp`: power recurrence for s^alpha, force and variational series,
  last-two-term step control and the exact final step were checked line by line; correct.
  The local error estimate is heuristic (as documented).
* `high_precision.py`: Newton uses the full 12-row return Jacobian; quotient construction
  (COM/momentum, f^T Omega, r^T Omega, Euclidean gauges, projection along f and r) is
  correct; Krein forms use the reduced symplectic form consistently.
* `threebody_expand.py`: forces and analytic Jacobian correct (also verified by tests); the
  `stats()` "closure" is a re-integration at rtol 1.5e-13, not the corrector residual.
* Minor: `finalize_results.boundary_report` hard-codes the bracket strings; the Makefile's
  85-digit boundary run intentionally uses binary64 interpolation seeds (independent seeds)
  but this was undocumented; `make test` relied on saved artifacts only.

## 7. External review of this session's work (findings fixed)

| Finding | Fix | Regression test |
|---|---|---|
| The half-period symmetry (`M = N^2`) was presented as excluding period doubling / bifurcations at the -1 touches; it only excludes a *spectral* stability change — 1:4 resonances of the half-period map can still create symmetry-breaking period-2T orbits (Gonchenko et al., Chaos 28, 043123 (2018)) | claims restricted to the spectral statement in README, paper, status, methods, audit, literature audit and `half_period_symmetry.py`; nonlinear analysis listed as an open item | — (text) |
| Parallel Make targets used a bare `wait`, so `make refine` could succeed (and summarize stale results) after worker failures | `tools/run_parallel.sh` waits on each PID and propagates failures; reruns move old outputs to `*.json.prev`; output directories configurable | `tests/test_review_fixes.py` (real `make refine` with failing fake workers) |
| `asymmetric_double.py` dropped failed Newton points, so a partial scan looked valid and negative-D remnants could mimic instability | failures and a coverage record are kept per row; crossings never bracket gaps; `asymmetric_summary.py` excludes incomplete rows from fits, closing estimate and mirror check; Phase E re-run | synthetic-failure tests in `tests/test_review_fixes.py`; artifact test requires `all_rows_complete` |
| `summarize_transitions.py` could count Taylor-vs-Taylor agreement as independent confirmation | confirmed digits now come only from Taylor-vs-GBS pairs, using a strict common prefix that excludes each run's last printed digit; K3 is quoted to 39 (not 40) digits under this rule | reviewer's scenario reproduced as a unit test |

Outputs superseded by the fixes (`asymmetric/main_window_pre_review.json`, its summary, and
`half_period_symmetry_pre_review.json`; identical numbers, but no coverage records or an
over-strong description) were kept in the local working copy and are not part of the repository. The
re-run Phase E scan is complete in every row (61/61 grid points, all crossings refined, no
failures) and reproduces the closing estimate, edge curvatures and mirror check exactly.
