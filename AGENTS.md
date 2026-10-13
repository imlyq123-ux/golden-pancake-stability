# Codex research handoff — Golden Pancake three-body orbit

**Status:** October 12, 2026. **This is an exploratory numerical research project, not a validated discovery or an arXiv-ready manuscript.**

## Mission

Continue the planar Newtonian three-body study of the **already known** Golden Pancake Drama isola with masses `(1, mu, 1)`. Reproduce and audit the existing double-precision mass continuation and Floquet calculations, then independently establish (or refute) the apparent boundaries of a nontrivial spectral-stability window near `mu ≈ 0.59904` and `mu ≈ 0.61784`. Only after verification, assess scientific novelty against the original atlas paper and relevant literature; if there is a real contribution, develop a short reproducible numerical manuscript.

Do the actual calculations and save scripts, data, logs, plots, and uncertainty estimates. Do not extrapolate precision from solver settings, inherit unsupported earlier claims, or present a numerical sample on a known family as a new orbit family.

## 1. Mathematical specification

Planar Newtonian gravity with `G = 1`:

`q_i'' = sum_{j != i} m_j (q_j - q_i) / |q_j - q_i|^3`.

Masses `m = (1, mu, 1)`. The Euler-phase shooting section is

```
q1(0) = (-1, 0)
q2(0) = ( 0, 0)
q3(0) = ( 1, 0)
v1(0) = v3(0) = (u, v)
v2(0) = -2*(u, v)/mu
```

For the full state use `X = (q1x,q1y,q2x,q2y,q3x,q3y,v1x,v1y,v2x,v2y,v3x,v3y)`. Solve for `(u,v,T)` such that **all twelve labeled position/velocity components** satisfy `X(T) - X(0) = 0`; do not substitute a shape-only or relative-period condition. This ansatz gives zero total momentum and angular momentum and fixes the initial center of mass.

High-precision seed at **mu = 0.613**, copied from the previous numerical session (must be independently verified from code in this project):

```
mu = 0.613
u  = 0.15148154107900583396835825988003861711916800727
v  = 0.14030301957068408729274415253719345009218867155
T  = 8.846608333332378559519969249427806733522940068
```

Construct the middle velocity from the relation `v2 = -2*(u,v)/mu` **at the working precision**, instead of copying a rounded decimal for `v2`. Report energy `E`, scale-invariant period `S = T*|E|^(3/2)`, minimum pairwise distance, conservation diagnostics and all twelve closure residuals.

## 2. Existing files in this folder

- `threebody_expand.py` — exploratory SciPy `DOP853` double-precision integration, Newton-type least-squares shooting, analytic gravitational Jacobian, complete `12 x 12` monodromy integration; standalone command-line interface.
- `pancake_stability_boundaries.py`, `pancake_upper_refine.py` — exploratory root bracketing using a Floquet discriminant.
- `pancake_convergence_check.py`, `pancake_convergence.log` — attempted tighter tolerance checks; **the captured log shows an integration failure near mu=0.599. These checks are NOT currently complete.**
- `pancake_continuation_A.json`, `_B.json`, `_C.json`, `pancake_transition_coarse.json`, `pancake_upper_boundary.json` — intermediate exploratory orbit records.
- `pancake_mass_sweep.csv` — 23 mass samples from about `0.585` to `0.630`, including `u,v,T,E,S,rmin,closure,D,max_lambda` and four selected multipliers.
- `pancake_conservation_tests.csv` — conservation diagnostics at 5 sampled masses.
- `pancake_stability_sweep.png`, `pancake_discriminant_sweep.png`, `pancake_nearest_approach.png` — plots from the exploratory double-precision data.
- `pancake_stability_note.tex`, `.pdf` — 3-page *working note*, NOT a publication-ready preprint.
- `golden_pancake_mu_0p613_trajectory.csv`, `_orbit.png` — earlier high-precision trajectory positions and plot; **the trajectory CSV stores positions only, not velocities or a reproducible multiprecision integrator.**
- `make_pancake_results.py` — plot/CSV generator (currently uses an absolute `/mnt/data` directory and hard-coded estimated boundaries).
- `README_pancake_research.txt` — earlier project README; treat its numerical numbers as provisional.

The **claimed** `~5.14e-40` full-state closure of the single mu=0.613 seed came from the preceding conversation. The executable high-precision Taylor/MP code and its full final-state log are **not included** here. Treat this as an **unverified historical claim**, not independently reproduced evidence. The supplied trajectory is insufficient to certify full 12-component return.

The expected `pancake_boundaries.json` file is missing from the prior archive; regenerate rather than assuming it exists. The continuation and plotting scripts contain absolute `/mnt/data` paths; make them project-relative before claiming cross-machine reproducibility.

## 3. Current exploratory findings — not certified

At double precision, shooting/integration residuals are typically `~4e-12` to `~2.2e-11` in full state. Some close passages have pair separations as small as `~0.008` on the sampled branch, so numerical conditioning is demanding. The tentative spectral-stability interval is

```
mu_lower ≈ 0.5990438
mu_upper ≈ 0.6178443
```

Suggested bracketing checks:

```
lower: mu = 0.5990, 0.5991
upper: mu = 0.6178, 0.6179
```

Approximate reported largest nontrivial Floquet moduli:

```
mu = 0.597: 1.1052   (unstable candidate)
mu = 0.605: 1.0000   (spectrally stable candidate)
mu = 0.613: 1.0000   (spectrally stable candidate)
mu = 0.618: 1.0613   (unstable candidate)
mu = 0.630: 2.3604   (unstable candidate)
```

The exploratory code selects **four multipliers farthest from 1** from the 12-dimensional monodromy. Eight theoretically neutral multipliers are associated with symmetries/integrals, but their numerical splitting can be significant. Confirm that this heuristic actually identifies the correct reduced quartet throughout continuation; near a degeneracy use a principled symmetry reduction or continuation/tracking of eigenvectors and eigenvalues.

For the four nontrivial multipliers, the exploratory discriminant is

`D = 2*sum(lambda_j^2) + 8 - (sum(lambda_j))^2`.

If they form exact reciprocal pairs `{a,1/a,b,1/b}`, this equals `[(a+1/a)-(b+1/b)]^2`. A crossing of `D=0` is **not by itself** proof of Hamiltonian–Hopf bifurcation: additionally verify unit-circle collision, off-circle complex quartet, algebraic/geometric multiplicity, and Krein signature, with appropriate symplectic reduction.

## 4. Required execution plan (in priority order)

### Phase A — reproducibility and code audit (do this first)

1. Create a reproducible Python environment (Python >= 3.11, NumPy, SciPy, Numba, pandas, Matplotlib; add pytest and high-precision backend). Save exact dependencies and versions.
2. Remove hard-coded `/mnt/data` paths, add CLI options for output directory, and retain the supplied baseline artifacts unchanged. Organize new outputs under `results/`, `figures/`, `tests/`, and `paper/`. Do not silently replace prior data.
3. Run a simple smoke test from the project root:

   ```bash
   OPENBLAS_NUM_THREADS=1 python threebody_expand.py \
       --mass-list 0.613,0.612,0.611,0.610,0.608,0.605 \
       --stability --out results/branch_recomputed_double.json
   ```

4. Create tests for pairwise-force antisymmetry/total force, the analytic `12x12` Jacobian versus finite differences or complex-step (where valid), Newton Jacobian versus finite differences, symmetry constraints, energy and angular-momentum conservation, and full labeled closure. Record integration `success`, status, function evaluations, solver convergence, runtime and precision in outputs.
5. **Fix and rerun the failed boundary convergence check** (`pancake_convergence.log` is a real exception, `Required step size is less than spacing between numbers.`). Diagnose whether the Newton trial steps leave the correct family, create near-collision ill-conditioning, or reveal other issues. Repeat tolerance/step-size comparisons on both sides of *both* candidate boundaries. Do not report convergence checks as passed until a successful log exists.
6. Independently reproduce the 23-point double-precision CSV and all three figures. Compare hashes or numerical diffs; clearly explain any disagreements.

### Phase B — independent multiprecision orbit verification

1. Implement or use a genuine arbitrary-precision integrator, preferably adaptive high-order Taylor integration / a compiled MPFR backend, with 60–100 decimal digits and documented local/global tolerance control. An *independent* high-order integration method is strongly preferred for cross-verification. Do not merely cast double-precision results to high precision.
2. Newton-correct `(u,v,T)` using multiprecision arithmetic for `mu=0.613`, both bracketing endpoints, and additional points on each side. Measure the **actual** 12-component terminal residual in a separate stricter full-period integration. Target `< 1e-28` componentwise when computationally feasible. Print all residual components and high-precision parameters, not just the largest error.
3. Quantify convergence with at least two Taylor orders/tolerances, working precisions and/or independent integrator methods; monitor energy, angular momentum, center of mass, momentum and nearest separation. If the target is not reached, report the achievable verified precision without inventing more digits.
4. Store machine-readable `results/high_precision_orbits.json` with equations/conventions, exact decimal inputs, precision, solver version, correction tolerances, return vector, conservation error, minimum distance, period and elapsed time. Add scripts/tests that recreate the JSON.

### Phase C — Floquet stability and transition verification

1. Integrate the variational system at true multiprecision. The present state uses **velocities**, not canonical momenta: if checking symplecticity use `Omega = [[0, M_mass], [-M_mass, 0]]` with `M_mass = diag(m1,m1,m2,m2,m3,m3)`, or transform from `v_i` to canonical `p_i = m_i*v_i` before applying the standard symplectic matrix. Check `Phi(T).T * Omega * Phi(T) - Omega` and expected reciprocal/conjugate multiplier symmetry.
2. Isolate the 4 nontrivial multipliers robustly, preferably by removing the translation, momentum, time/energy and rotation/angular-momentum neutral directions with an explicitly documented reduction. Track modes continuously in `mu`. Avoid using double-precision eigensolvers on supposedly high-precision matrices.
3. With a verified `D(mu)` and monodromy, bracket and refine both transitions. Use independent initial seeds and rising precision (e.g. 50, 70, 90 decimal digits) and report **uncertainty brackets**, not unsupported 10-digit point estimates.
4. At each candidate boundary inspect whether the colliding modes have opposite Krein signatures and the correct quartet departure; describe the result as a **candidate spectral transition** unless full evidence supports a stronger label. Distinguish spectral stability from nonlinear/Lyapunov stability.
5. Save `results/floquet_sweep_high_precision.csv`, `results/stability_boundaries.json`, and figures showing all four nontrivial multipliers in the complex plane, `D(mu)`, spectral radius and uncertainty bands.

### Phase D — literature and novelty audit

Read and compare all relevant text/tables/appendices/supplements from:

- Existing family: https://www.threebodyorbits.com/family/golden-pancake-drama
- Atlas computations: https://www.threebodyorbits.com/about
- Original research draft: https://www.threebodyorbits.com/papers/golden-pancake-drama-2026.pdf
- Vasiljević, Raonić & Dmitrašinović, *An island of linearly stable non-hierarchical unequal mass periodic three-body orbits*, *New Astronomy* **100**, 101969 (2023).

The atlas currently publicly lists only 7 representatives and describes the isola but **that does not prove the authors have never analyzed its full stability window**. Determine whether exact boundary masses, spectra, bifurcation types, or related stability intervals have already been published. Do not claim novelty without checking the original draft and adjacent literature. Cite prior authors appropriately.

### Phase E — stronger new research if boundary study is insufficient

Explore a genuinely higher-dimensional parameter study with asymmetric outer masses, such as `(1+delta,mu,1-delta)` for small `delta`. **Do not use the symmetric 3-parameter Euler-phase shooting ansatz unchanged for unequal outer masses**: derive a valid full periodic shooting/phase-fixing continuation setup with center-of-mass and rotation handling. Locate a 2D stability region or bifurcation curve, quantify its uncertainty, and compare with known literature. This is optional and should follow Phases A–D, not replace them.

## 5. Deliverables / acceptance criteria

- `README.md` with tested installation and run commands on a clean machine.
- `requirements.txt` or locked environment; portable code; tests with explicit pass/fail results.
- Raw and processed output data, each with numerical precision, solver settings and reproducibility provenance; no untracked manual edits.
- High-precision orbit(s) with independently reproduced **full** phase-space closure, conservation diagnostics and a plot regenerated from the verified integration.
- Verified lower/upper stability transition brackets and a documented reliable treatment of the nontrivial Floquet quartet (or an explicit explanation of why these remain unresolved).
- `paper/paper.tex` plus figures **only once** novelty and verification are established. If not established, write `results/research_status.md` stating limitations and next experiments; do not produce an assertive discovery manuscript.

Recommended first Codex action: **read `threebody_expand.py`, `pancake_stability_note.tex`, and `pancake_convergence.log`; create `results/`; run the Phase A smoke test; explain/fix the failed convergence check, then verify numerical results before continuing.**
