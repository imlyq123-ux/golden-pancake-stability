# Research status — second session (2026-10-12)

**Bottom line.** The first session's central numbers are correct and are now confirmed by an
independent multiprecision implementation and by the atlas authors' own published data. The
stability picture is richer than the original two-boundary window: the Golden Pancake isola
has three spectrally stable stretches on its AABAbbaabaBB segment (plus a short one on the
(ABab)^3 segment), separated by non-semisimple Krein collisions and a fold; all orbits carry
a half-period reflection symmetry that makes the monodromy a square; and the main window
closes under a small outer-mass asymmetry. Nothing is interval-certified.

## 1. Audit of the first session

Details in [`audit_report.md`](audit_report.md). Confirmed: all closures, D values, quartets,
Krein signs and the two 1e-8 brackets. Corrected: the "independent time-reversed" check is a
bitwise mirror of the forward run (no independent evidence on truncation error); the 65/85-digit
"D change ~8e-46" compared runs stopped at the same Newton iterate (true D error up to ~8e-43,
7e-38 for the uncorrected 0.613 seed); recorded source hashes refer to earlier (numerically
equivalent) code versions; the literature audit was incomplete (draft not read).

## 2. Independent verification

* **New pipeline**: Gragg–Bulirsch–Stoer extrapolation integrator in MPFR (`tools/gbs_mpfr.cpp`)
  + mpmath Newton + reduction-free Floquet analysis (`independent_mp.py`). Seeds from binary64
  continuation only.
* **Against the first session** (15 masses, 80 digits, closures 6e-63 ... 3e-62): on the 85-digit
  bracket orbits `|du| = 7e-65`, `|dT| = 6e-63`, `|dD| = 5e-59`, `|dPhi| = 1e-60`; all Krein signs
  agree. The supplied seed's return vector (5.14e-40) is reproduced component by component to
  5.5e-49 — it is a property of the 47-digit seed.
* **Against the atlas** (independent 160-bit computation by the family's discoverers): at the
  draft's 42-digit mass our `(u, v, T)` differ from its Table 1 by 3.3e-38, 1.9e-38, 1.5e-36
  (draft closure 1.5e-36) and the four multipliers by 4.3e-20 (19 printed decimals); catalogue
  `stab` 1.1637 / 1.0 / 2.3815 vs ours 1.1636633 / 1.0 / 2.3814702 for the three atlas orbits on
  this branch. See [`independent_crosscheck.md`](independent_crosscheck.md).
* **Binary64 re-derivation**: one continuation chain from the single seed reproduces the archived
  23-point sweep (|du|, |dv| <= 4e-13, |dT| <= 1.6e-11, |dD| <= 3e-8); the "farthest-from-1"
  heuristic picks the correct quartet everywhere scanned (max difference 6e-8 vs the quotient).

**Corrected mu = 0.613 seed.** The supplied 47-digit seed is correct to about 40 digits only.
Newton-corrected independently by both pipelines (GBS-80 closure 1.4e-62; Taylor-85 closure
6.7e-76, `results/seed_correction/`), the parameters agree to 4.7e-64, 1.8e-64, 1.4e-63:

```
u = 0.15148154107900583396835825988003861711921369614127787182422...
v = 0.14030301957068408729274415253719345009221621399465293205890...
T = 8.8466083333323785595199692494278067335242240953689621171902...
E = -1.5442767428070372689259964122887886618...,  S = T|E|^(3/2) = 16.977138689603261341576613470326...
```

The sampled minimum pair separation is 0.0127596 (accepted-step nodes; an upper estimate).

## 3. Transitions (multiprecision)

Full values, pipeline-by-pipeline estimates and verified sign brackets (`mu* -/+ 1e-30`, all
with sign changes) are in [`stability_transitions.md`](stability_transitions.md) /
`stability_transitions.json`; raw records in `boundary_refinement/`.

| | kind | mu* (truncated to Taylor-vs-GBS agreement; Taylor-vs-Taylor not counted) | dD/dmu |
|---|---|---|---|
| K0 | Krein | 0.580094414288468996548897180964879640274... (39 digits, GBS-80 vs Taylor-50) | -22.25 |
| K1 = mu_L | Krein | 0.599043784208899851220728860462068419275413725661780425744347... (first 60 digits identical in GBS-80 and Taylor-90; relative difference 1.1e-60) | +74.874 |
| K2 = mu_R | Krein | 0.617844268209063178727366710600417996702072689168837458676169... (60 digits, GBS-80 vs Taylor-90) | -305.88 |
| K3 | Krein | 0.634197066930554588552605424098150004643... (39 digits, GBS-80 vs Taylor-50) | +5839.8 |
| T1 ((ABab)^3) | s = +2 | 0.475901444550864118752069986423460087896... (39 digits, GBS-80 vs Taylor-50) | +7.83 (q(+1)) |

Transversality: the transition functions cross zero with the nonzero slopes shown.

## 4. Type of the main transitions (K1, K2)

`transition_analysis.json`, `figures/transition_scaling.png`. Over `delta = |mu - mu*|` from
1e-4 to 1e-40 (80-digit GBS orbits): stable-side splitting `~ delta^0.5000`, unstable-side
`rho - 1 ~ delta^0.5000` (K1) / `delta^0.5001` (K2), normalized Krein forms `~ delta^0.501`,
`1 - |<z1,z2>| ~ delta^1.000`; opposite Krein signatures on every stable-side sample. At `mu*`,
`Phi - lambda_bar I` has singular values 2.9e-61 / 0.020 (K1) and 2.5e-62 / 0.092 (K2): geometric
multiplicity one, algebraic two. This is the linear signature of a Hamiltonian–Hopf
bifurcation (non-semisimple 1:1 resonance with a complex quartet on one side). Nonlinear
normal-form coefficients (super/subcritical) were not computed.

## 5. Stability along the whole isola (binary64 pseudo-arclength scans)

`isola_scan_summary.json`, `figures/isola_stability_overview.png`.

AABAbbaabaBB segment (272 points, mu in [0.5393, 0.6346]):
* stable for mu from the scan end 0.5393 up to K0 (no transition detected in between; the
  segment continues toward the collision limit C- ~ 0.5066, not reached);
* complex quartet (K0, K1); stable [K1, K2] (main window, contains the atlas representative);
  complex quartet (K2, K3); stable [K3, mu_max]; fold at mu_max ~ 0.6346043 (parabolic fit of
  scan points; atlas: 0.6346042508) with a pair through +1 (s = +2); beyond the fold one real
  pair (unstable). The far-side scan stops at mu ~ 0.634486 (corrector limit).
* Near mu ~ 0.5695 and 0.5481 a pair sum touches -2 without crossing, and near mu ~ 0.5619 two
  pairs pass through each other (same Krein sign). Section 6 shows why the spectrum stays on the
  unit circle at these points (no change of spectral stability); possible nonlinear
  bifurcations there were not studied.

(ABab)^3 segment (478 points, mu in [0.4689, 0.6161]): unstable (one hyperbolic pair, rho up to
~27) except a short stretch after the mu_min fold, from T1 (s = +2 crossing at mu ~ 0.4759) to
the scan end at mu ~ 0.4797 (toward C-). Consistent with the atlas `stab` values of the four
(ABab)^3 representatives (1.38 ... 23.4).

## 6. Half-period reflection symmetry

The atlas computes this family in an "Euler-half-twist symmetric sector", so a half-period
symmetry is implicitly known; what follows makes it explicit and derives its spectral
consequences. `half_period_symmetry.json`. For all 15 independently computed orbits `X(T/2) = sigma X(0)`
with a reflection `sigma` (defects 1e-63 ... 8e-63 at 80 digits; axis 131.9 ... 151.7 degrees),
and `|N^2 - M| <= 3.3e-59` with `N = sigma DPhi_{T/2}`. Consequences: `s_M = s_N^2 - 2` (identity
verified to <= 1e-56), so real pair sums of `M` can only touch -2 (the touches above are
`N`-multipliers passing `+-i` -- 1:4 resonances of `N` -- with `s_N = 0.0038` and `-0.0028` at the
sampled points). This excludes a change of *spectral* stability there, not nonlinear
bifurcations: 1:4 resonances generically create period-4 orbits of the resonant map, including
non-symmetric ones (Gonchenko, Gonchenko, Ovsyannikov & Vieiro, Chaos 28, 043123 (2018), arXiv:1711.09447), i.e. possibly period-2T orbits breaking the half-period symmetry;
these were not investigated. Further, `|lambda_M| = |lambda_N|^2`, and `D_N` changes sign at the same bracket endpoints as `D_M`
(e.g. -6.35e-8 / +8.74e-8 across the old K1 bracket), so K1, K2 (and K0) are genuine Krein
collisions of the half-period map. T1 is an `s_N = +2` crossing.

## 7. Phase E: outer masses (1 + delta, mu, 1 - delta) (binary64, exploratory)

`asymmetric/main_window*.json`, `figures/asymmetric_stability_region.png`. A general collinear
section with the rotation angle as an unknown (no symmetric ansatz) shows that the orbits stay
exactly periodic (`|theta| <= 4e-11`) — as the reflection symmetry predicts for any masses. The
main window shrinks: `mu_L(delta) = mu_L + 126.2 delta^2 + ...`, `mu_R(delta) = mu_R - 39.1 delta^2 + ...`
(Richardson-extrapolated difference quotients), and closes at `|delta| ~ 0.0086` near `mu ~ 0.611`
(between grid values 0.0085 and 0.0090). `delta = 0` reproduces the multiprecision K1, K2 to
3e-11 / 4e-11, and `delta -> -delta` (relabelling + rotation + time reversal) gives identical
transitions to 3e-11. Not yet verified in multiprecision.

## 8. Novelty (see `literature_audit.md`)

The family, its isola and the stability of one representative are the atlas's (2026). The atlas
draft explicitly asserts the stability of that one orbit only; no inspected source gives the
transition masses, their type, the additional stable stretches, the spectral consequences of
the half-period symmetry (the symmetry itself underlies the atlas's "Euler-half-twist" sector)
or the asymmetric-mass region. A cautious draft manuscript is in `../paper/` for the authors'
review; it credits the atlas and claims only numerical evidence.

## 9. Limitations and next steps

1. Interval-arithmetic (computer-assisted) proof of existence and of the spectral claims at
   least at K1/K2 brackets — none attempted.
2. Hamiltonian–Hopf normal-form coefficient (nonlinear type) at K1, K2.
3. Multiprecision verification of the Phase E curve and of the fold location; extend Phase E to
   the other stable stretches and to non-zero angular momentum.
4. Regularized continuation through the binary-collision limits to close the isola map; the
   scans stop at r_min ~ 1.6e-3 ... 2.2e-3.
5. Full text of Vasiljević, Raonić & Dmitrašinović (2023) not read (paywalled).
6. Nonlinear analysis at the 1:4 resonances of the half-period map (mu ~ 0.5695, 0.5481): look
   for bifurcating period-2T orbits that break the half-period symmetry.
