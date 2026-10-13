# Numerical methods and interpretation

Conventions: planar Newtonian three-body problem, `G = 1`, masses `(1, mu, 1)` (Phase E:
`(1+delta, mu, 1-delta)`), state `X = (q1x, q1y, q2x, q2y, q3x, q3y, v1x, v1y, v2x, v2y, v3x, v3y)`
in the velocity convention. Euler section: `q = (-1,0), (0,0), (1,0)`, `v1 = v3 = (u, v)`,
`v2 = -2 (u, v)/mu`, formed at working precision from exact decimal inputs. Every closure is the
full labelled 12-component return `X(T) - X(0)`; no body permutation, rotation or shape-only
condition is used. The archive `threebody_codex_project/` is preserved unchanged.

## 1. Binary64 tools

* `threebody_expand.py`: SciPy DOP853, analytic force Jacobian, full 12x12 variational
  equations; damped least-squares Newton (first session). With Numba installed (this session)
  the right-hand sides are compiled; results agree with the NumPy backend to roundoff.
* `dense_sweep_double.py`: predictor-corrector continuation in steps of 0.0005 from the single
  `mu = 0.613` seed (no archived per-mass seeds), nontrivial quartet from an explicit symplectic
  quotient (orthonormal null-space basis by SVD; projection along the flow and rotation
  generators), compared with the legacy "4 multipliers farthest from 1" heuristic.
* `isola_scan_double.py`: pseudo-arclength continuation in `(u, v, T, mu)` with the 12x4
  shooting Jacobian. The mass column comes from the inhomogeneous variational equation
  `dS/dt = A(t) S + dF/dmu`, `S(0) = dX(0)/dmu` (`dF/dmu` is nonzero only in the accelerations
  of bodies 1 and 3), validated against central differences (relative 3e-8). The scaled
  Jacobian has rank 3 with a clean one-dimensional null space (branch tangent). Corrector
  tolerance `2e-10` (closure floor of DOP853 at `rtol = 3e-13` near small `r_min`); the scans
  stop where even this cannot be met (`r_min ~ 1.6e-3 ... 2.2e-3`, before the collision limits).
  Points are classified by the pair sums `s = lambda + 1/lambda`: stable iff both are real and in
  `(-2, 2)`; transitions are D sign changes (Krein), `s = -2` (period doubling) or `s = +2`
  crossings (tangent), located by linear interpolation (exploratory only).

## 2. First-session multiprecision pipeline (MPFR Taylor + Decimal)

`tools/taylor_mpfr.cpp`: adaptive Taylor series of normalized coefficients, `s^alpha`
recurrence `g[k] = sum_{l=1..k} ((alpha+1) l - k) s[l] g[k-l] / (k s[0])`, last-two-term
step control (safety 0.65), exact final step; 2 shooting columns or the 12x12 monodromy.
`high_precision.py`: Decimal least-squares Newton on all 12 return components, separate stricter
verification run, Floquet quotient (COM/momentum rows, `f^T Omega`, `r^T Omega`, Euclidean
gauges, projection along `f, r`), reduced characteristic polynomial roots and Krein forms
`kappa = Re(i z^H Omega_r z)`. Reviewed line by line in this session; no numerical errors found.

## 3. Independent multiprecision pipeline (this session)

Written without reusing any first-session integration, shooting or Floquet code.

* **Integrator** `tools/gbs_mpfr.cpp`: Gragg-Bulirsch-Stoer extrapolation — modified midpoint
  rule with Gragg smoothing, step-number sequence `n_j = 2j`, Aitken-Neville extrapolation in
  `h^2`, `k` stages (order `2k`), error estimate `|T_kk - T_k,k-1| <= tol (1 + |y_i|)`, step factor
  `0.94 (0.65/err)^(1/(2k-1))` clamped to [0.2, 4], exact final step; working precision
  `digits + 20` decimal digits; official MPFR 4.2.1 headers. It needs only forces and the force
  Jacobian, so it is algorithmically independent of the Taylor recurrences.
* **Shooting** `independent_mp.py`: mpmath least-squares Newton on the 12 return components with
  the 2-column variational solution and the vector field at `X(T)`; seeds only from binary64
  continuation of the originally supplied sweep (or from binary64 scans), never from the
  Taylor/Decimal records. Typical settings: 80 digits, Newton `tol = 1e-62` with 20 stages,
  verification `tol = 1e-66` with 24 stages, target closure `1e-55`.
* **Reduction-free Floquet analysis.** For the exact monodromy the eight trivial multipliers
  equal 1 (COM/momentum, time/energy, rotation/angular-momentum blocks), so the power sums
  `p_k = tr(Phi^k) - 8` (k = 1..4) give the nontrivial quartic by Newton's identities and
  `D = 2 p_2 + 8 - p_1^2` without any reduction or eigenvalue selection. Palindromic defects
  `|e3 - e1|`, `|e4 - 1|` measure consistency. All 12 eigenvalues and right eigenvectors come
  from `mpmath.eig`; Krein forms use full-space eigenvectors with
  `Omega = [[0, M], [-M, 0]]`, `M = diag(1,1,mu,mu,1,1)` — eigenvectors of nontrivial multipliers
  are Omega-orthogonal to the generalized 1-eigenspace, so their signs equal the quotient signs.
* **External comparison** `crosscheck_report.py`: against the Taylor records and against the
  atlas draft (42-digit Table 1, 19-digit multipliers) and catalogue (`stab`, `E`, `r_min`, ICs).

## 4. Transition refinement and classification

* `refine_boundaries.py`: secant iteration on a transition function — `D` (Krein collision),
  `q(-1) = prod(1 + lambda_i)` (period doubling) or `q(+1) = prod(1 - lambda_i)` (tangent) —
  where every evaluation is a fresh Newton-corrected orbit plus a separate monodromy
  integration. Pipelines: `taylor50/70/90` (Taylor + Decimal quotient at 50/70/90 digits,
  verification tolerances `1e-44/1e-60/1e-75`) and `gbs80` (GBS + mpmath). After convergence,
  `F` is evaluated at `mu* -/+ 1e-30`; a sign change gives an *empirically verified* bracket
  (not an interval enclosure). `summarize_transitions.py` counts *confirmed* digits only between
  independent implementations (a Taylor run versus the GBS run): the common decimal prefix,
  excluding each run's last printed digit (which may be rounded). Taylor-versus-Taylor agreement
  is reported separately as a precision-consistency check and never adds confirmed digits.
* Parallel batches (`make refine`, `make isola-scan`) run through `tools/run_parallel.sh`, which
  waits on every process and fails the target if any job fails, so summaries never run on a
  partial batch; reruns move a previous output aside (`*.json.prev`) before starting.
* `transition_analysis.py`: orbits at `mu* +/- 10^-k` (k = 4 ... 40). For a non-semisimple
  collision one expects `|lambda1 - lambda2| ~ delta^(1/2)` and normalized Krein forms
  `~ delta^(1/2)` on the stable side, `rho - 1 ~ delta^(1/2)` on the unstable side,
  `1 - |<z1, z2>| ~ delta`, and at `mu*` one tiny and one O(1) singular value of
  `Phi - lambda_bar I` (geometric multiplicity 1). A semisimple collision would give linear
  splitting and two small singular values.

## 5. Half-period reflection symmetry (`half_period_symmetry.py`)

Every orbit on the isola satisfies `X(T/2) = sigma X(0)` with `sigma` a reflection of all
positions and velocities across a line through the centre of mass (axis angle depends on mu;
at `mu = 0.613` it is 139.40 degrees). Reflections commute with the Newtonian flow for any
masses, so `X(t + T/2) = sigma X(t)` and `M = N^2` with `N = sigma DPhi_{T/2}(X(0))`.
Consequences: `s_M = s_N^2 - 2`, so a real `s_M` can only *touch* `-2` (when `N`-multipliers pass
`+-i`, a 1:4 resonance of `N`), so the spectrum of `M` cannot leave the unit circle through -1
inside the symmetric family. This is spectral only: 1:4 resonances can create symmetry-breaking
period-2T orbits (nonlinear bifurcations; Gonchenko, Gonchenko, Ovsyannikov & Vieiro, Chaos 28, 043123 (2018), arXiv:1711.09447), which were not studied.
`|lambda_M| = |lambda_N|^2`, so `M`
and `N` are spectrally stable together; `N` has trivial multipliers `+1` (time/energy,
translation along the axis) and `-1` (rotation/angular momentum, translation across the axis),
four each, so its nontrivial power sums are `tr(N^k) - 4 - 4(-1)^k`. The script integrates to
`T/2` with the GBS integrator, measures the symmetry defect and `|N^2 - M|`, and reports `D_N`
and the identity defect of `s_M = s_N^2 - 2`.

## 6. Phase E: unequal outer masses (`asymmetric_double.py`, binary64)

The symmetric ansatz is replaced by a general collinear section on the x-axis
(`q1 = (x1, 0)`, `q2 = (a, 0)`, `q3 = (x1 + 2, 0)`, `x1 = -(mu a + 2 m3)/(m1 + m3)`), velocities
`v1 = (u1, w1)`, `v3 = (u3, w3)` with `w3` fixed by `L = 0` and `v2` by zero momentum; unknowns
`(a, u1, u3, w1, T, theta)`; equations `X(T) - R_theta X(0) = 0` (12 components, rank 6). The
rotation angle `theta` is solved for, not assumed; it stays at `|theta| <= 4e-11`, as expected
because the reflection symmetry (and with it `L = 0` and exact periodicity) persists for any
masses. Stability uses `R_{-theta} Phi` on the quotient with the mass-weighted form.
Derivatives of `X(0)` use complex-step differentiation. The `delta -> -delta` relabelling
symmetry is checked numerically (transitions agree to 3e-11). Newton failures are kept in each
scan row together with a coverage record; a row is *complete* only if every grid point was
computed and every crossing refined, and only complete rows enter the edge fits, the
window-closing estimate and the mirror check (a failed or missing point must never be read as
instability).

## 7. What the numbers mean (limits)

* All brackets and digits are *empirical*: two independent integrators and Newton/Floquet
  codes agree; local error control is heuristic; nothing is interval-certified; existence of the
  periodic orbits is not proved.
* Binary64 scans (isola, Phase E) are exploratory; their transition locations are accurate to
  roughly `1e-8 ... 1e-10` in mu (validated where multiprecision values exist).
* Sampled minimum distances are upper estimates of the true minimum.
* Spectral (linear) stability only; no nonlinear, KAM or Lyapunov statement. The
  Hamiltonian–Hopf label refers to the linear signature (non-semisimple Krein collision with a
  complex quartet on one side); normal-form coefficients were not computed.
