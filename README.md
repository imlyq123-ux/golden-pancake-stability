# Golden Pancake three-body orbit: stability along the isola

Reproducible numerical study of the linear (spectral) stability of the **already known**
Golden Pancake Drama family (Three Body Orbits atlas, 2026) — planar Newtonian three-body
problem, `G = 1`, masses `(1, mu, 1)` — and of its continuation to unequal outer masses.

**Status (2026-10-12, second session): numerically verified with two independent
multiprecision pipelines and against the atlas's own data; not interval-certified.** The
orbit family is not new; what is new is the location and type of its stability transitions.
Read [`results/research_status.md`](results/research_status.md) (findings and limits),
[`results/audit_report.md`](results/audit_report.md) (audit of the first session),
[`results/methodology.md`](results/methodology.md) and
[`results/literature_audit.md`](results/literature_audit.md). Research brief: [`AGENTS.md`](AGENTS.md).

## Main results

**Spectral transitions** (multiprecision secant refinement; values are *truncated* to the digits
on which the independent Taylor and GBS implementations agree, excluding each run's last printed
digit; Taylor-vs-Taylor agreement is not counted; full values and the verified `mu* -/+ 1e-30` sign
brackets are in [`results/stability_transitions.md`](results/stability_transitions.md)):

| | type (linear) | mu* | pipelines |
|---|---|---|---|
| K0 | Krein collision | 0.580094414288468996548897180964879640274... | GBS-80 vs Taylor-50 |
| K1 = mu_L | Krein collision, non-semisimple | 0.599043784208899851220728860462068419275413725661780425744347... | GBS-80 vs Taylor-90 (Taylor-50/70 vs GBS-80: 39/53 digits) |
| K2 = mu_R | Krein collision, non-semisimple | 0.617844268209063178727366710600417996702072689168837458676169... | GBS-80 vs Taylor-90 (Taylor-50/70 vs GBS-80: 39/58 digits) |
| K3 | Krein collision | 0.634197066930554588552605424098150004643... | GBS-80 vs Taylor-50 |
| fold | tangent (s = +2) at mu_max | ≈ 0.6346043 (binary64; atlas: 0.63460425) | scan |
| T1 | s = +2 crossing, (ABab)^3 segment | 0.475901444550864118752069986423460087896... | GBS-80 vs Taylor-50 |

K1 and K2 agree between the independent 80-digit GBS pipeline and the 90-digit Taylor
pipeline to about 60 significant digits (relative differences 1.1e-60 and 2e-61).

**Stable stretches of the AABAbbaabaBB segment** (all four nontrivial multipliers on the unit
circle): `mu in [<=0.5393 (scan end), K0]`, the main window `[K1, K2]` containing the atlas
representative, and a narrow near-fold stretch `[K3, mu_max]`. The (ABab)^3 segment is
unstable except for a short stretch just past its `mu_min` fold (`mu ~ 0.4759 ... >= 0.4797`).
Figure: [`figures/isola_stability_overview.png`](figures/isola_stability_overview.png).

**Transition type.** At K1 and K2 the colliding multipliers have opposite Krein signatures
and split like `delta^0.5000` (stable side) / leave the circle like `delta^0.5000` (unstable
side) over `delta = |mu - mu*| = 1e-4 ... 1e-40`; at `mu*` the matrix `Phi - lambda I` has one
singular value `~1e-61` and the next `~0.02-0.09` — a 2x2 Jordan block (the linear signature of a
Hamiltonian–Hopf bifurcation). Figure: [`figures/transition_scaling.png`](figures/transition_scaling.png).

**Half-period reflection symmetry** (implicit in the atlas's "Euler-half-twist" sector; made
explicit here). Every orbit satisfies `X(T/2) = sigma X(0)` for a
reflection `sigma` (defect `<= 8e-63` at 80 digits), so the monodromy is a square,
`M = N^2`, and `M`, `N` are spectrally stable together. Near mu ≈ 0.5695 and 0.5481 a pair of
multipliers touches −1 without crossing (an `N`-pair passes ±i, a 1:4 resonance of `N`), so the
spectral stability does not change there. This is a spectral statement only: 1:4 resonances can
still create period-2T orbits that break the half-period symmetry (nonlinear bifurcations,
[Gonchenko et al. 2018](https://arxiv.org/abs/1711.09447)); these were not investigated.

**Unequal outer masses (Phase E, binary64).** For masses `(1+delta, mu, 1-delta)` the orbits
stay exactly periodic (`|theta| <= 4e-11`), the main window shrinks as
`mu_L + 126.2 delta^2`, `mu_R - 39.1 delta^2`, and closes at `|delta| ≈ 0.0086` near `mu ≈ 0.611`.
Figure: [`figures/asymmetric_stability_region.png`](figures/asymmetric_stability_region.png).

**Verification.** The new GBS-extrapolation integrator + mpmath Newton/Floquet code shares
no code with the first session's Taylor/Decimal pipeline; on the 85-digit bracket orbits they
agree to `|dD| ~ 5e-59`, `|dPhi| ~ 1e-60`. Independently of both, our orbit at the atlas's
42-digit mass agrees with the atlas draft's Table 1 to `3e-38` (u), `2e-38` (v), `1.5e-36` (T)
and with its 19-digit multipliers to `4e-20`; the catalogue `stab` values 1.1637 / 1.0 / 2.3815
are reproduced. See [`results/independent_crosscheck.md`](results/independent_crosscheck.md).

## Installation (tested on Linux x86-64, Python 3.12.3, GCC 13.3, MPFR 4.2.1)

```bash
python3 -m venv .venv
.venv/bin/python -m pip install -r requirements.lock.txt
make deps        # only if /usr/include/mpfr.h is missing: fetches libmpfr-dev/libgmp-dev headers into .deps/ (no root)
make all         # builds tools/taylor_mpfr and tools/gbs_mpfr against the system libmpfr.so.6
make test
```

`make deps` uses `apt-get download` (Debian/Ubuntu). On other systems install the MPFR and GMP
development packages instead. `tools/taylor_mpfr` can also build with the minimal
`tools/mpfr_abi.hpp` fallback (LP64 only).

## Reproducing the results

| Target | What it does | Approx. wall time in this session (16 cores) |
|---|---|---|
| `make phase-a` | binary64 smoke test, convergence check, 23-point sweep (first session; overwrites its outputs) | not re-run here |
| `make dense-sweep` | single-seed continuation, 110 masses, quotient vs heuristic | ~8 min (NumPy backend) |
| `make isola-scan` | pseudo-arclength scans of both isola segments (parallel) + summary | ~15 min (Numba) |
| `make phase-bc` | first-session Taylor/Decimal runs at 65/85 digits | first session: ~30 min |
| `make independent` | GBS/mpmath re-solution of 15 masses, cross-check report, half-period symmetry | ~15 min |
| `make refine` | secant refinement of K0–K3 and T1 with several pipelines | ~2.5 h (dominated by 90-digit runs) |
| `make transitions` | square-root scaling and Jordan test around K1, K2 | ~12 min |
| `make phase-e` | unequal outer masses, 2D scan + summary | ~4 min |
| `make figures` | transition table and figures from saved results | <1 min |
| `make test` | all tests (physics, Jacobians, integrators, Floquet code, artifacts) | ~1 min |

Every script has `--help`.
Raw multiprecision outputs (`results/mp_raw`, `results/gbs_raw*`, `results/boundary_refinement/raw`,
seed-correction CSVs; about 1.4 GB) are **not included in the repository** — every reported number
is in the processed JSON records, which also store SHA256 hashes of the raw files (all hashes are
listed in `results/artifact_manifest.json`). The make targets regenerate them. Source snapshots of
every code version used are in `results/source_snapshots/`.

## Layout

| Path | Contents |
|---|---|
| `threebody_codex_project/` (+ `.zip`) | original archive, unchanged (SHA256 manifest checked by tests) |
| `threebody_expand.py`, `pancake_*.py`, `make_pancake_results.py` | binary64 tools (first session, portable) |
| `high_precision.py`, `tools/taylor_mpfr.cpp`, `finalize_results.py`, `verify_time_reversal.py` | first-session MP pipeline |
| `tools/gbs_mpfr.cpp`, `independent_mp.py`, `crosscheck_report.py` | independent MP pipeline (this session) |
| `refine_boundaries.py`, `summarize_transitions.py`, `transition_analysis.py`, `half_period_symmetry.py` | transition refinement and classification |
| `dense_sweep_double.py`, `isola_scan_double.py`, `isola_scan_summary.py` | binary64 continuation and stability scans |
| `asymmetric_double.py`, `asymmetric_summary.py` | Phase E (unequal outer masses) |
| `make_figures.py` | figures of this session |
| `results/` | records, logs, reports; `figures/`; `tests/`; `paper/` (manuscript draft) |

## Limitations

No interval-arithmetic proof (existence, closure or spectrum); local error control is
heuristic and confidence comes from agreement of independent methods and precisions. Binary64
scans can miss features smaller than their step; they stop before the binary-collision
limits of the isola. Spectral stability only (no nonlinear/KAM statement). Phase E is binary64.
Novelty was assessed against the atlas draft, catalogue and the literature listed in
`results/literature_audit.md`; worldwide priority cannot be proven.
