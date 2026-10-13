# Research status — numerical verification, not publication-ready

Phase A reproduction and the planned symmetric-mass high-precision orbit /
Floquet checks are complete for the saved samples. Phases D/E remain incomplete
or deferred. No existence proof, nonlinear stability result, exact bifurcation
classification, or novelty claim is made.

## Phase A

- The original supplied files match their SHA256 manifest and remain unchanged.
- The six-point smoke test completed successfully. The 23-point binary64 sweep
  reproduced all compared numerical columns exactly in this environment, using
  the original archived seeds and fresh integrations. This is not an independent
  seed search; the already small residuals did not require parameter changes.
- All three original plot types were regenerated from recalculated data.
  PNG hashes intentionally differ because annotations now describe provisional
  candidates and measured sign brackets instead of unsupported precise roots.
- The original working note's claim of completed tolerance checks was unsupported
  by its archived exception. The replacement check completed **24** successful
  integrations: four endpoint masses, three tolerances, and two maximum steps,
  with consistent signs at each endpoint. Detailed solver status/nfev are saved.
- The exact archived nearest-seed exception did not recur. A larger 0.613→0.599
  jump reproduced it after an off-branch near-collision trial reached separation
  ≈2.12e-7. Both successful nearest seeds and failed trial history are preserved.
  Damping, trial rejection and small mass steps address this demonstrated failure
  mechanism. The exact historical trigger remains uncertain.
- A clean isolated Python environment installed successfully from cached wheels.
  Exact Python package versions, compiler and MPFR version are recorded. Numba
  and pandas were unavailable offline; NumPy and stdlib CSV are documented
  substitutes, and MPFR/Decimal provide the arbitrary-precision backend.

## Phases B/C

Eight mass samples and four narrow-bracket endpoints were computed at each
of **65 and 85 decimal digits** (24 corrected/verified orbit records). Correction
and verification use different Taylor orders and tolerances; every verification
is a separate full-period integration. The seed at mu=0.613 retains the supplied
decimal parameters because its measured residual already meets the target.

The largest actual full-state componentwise return across these runs is
**5.14368E-40**, below 1e-28. All twelve components and complete terminal
states are printed in the machine-readable JSON records. Maximum sampled
energy drift is 1.55040E-56. Maximum full mass-weighted symplectic defect
is 6.70341E-52. Sampled nearest distances are upper estimates of
the true minimum and do not certify absence of collisions between nodes.

Across 12 matched mass samples, the largest change in D between precision/order
runs is 8.28769E-46; the largest optimally matched quartet change is
4.98233E-44. These are empirical convergence measurements, not rigorous
error enclosures or evidence that every printed digit is significant.

The full monodromy is reduced explicitly by COM/momentum constraints and a
time/energy, rotation/angular-momentum quotient. Actual reduced polynomial
roots are computed in Decimal complex arithmetic, with eigenvector residuals,
reciprocal/conjugate symmetry, determinant and symplectic checks. Saved mode
tracking matches nearby eigenvalues; individual labels at exact collisions
remain ambiguous.

| Candidate transition | Numerical mass bracket | Width | Endpoint D values |
|---|---|---|---|
| lower | [0.59904378, 0.59904379] | 1E-8 | -3.15137E-7, +4.33603E-7 |
| upper | [0.61784426, 0.61784428] | 2E-8 | +2.51101E-6, -3.60665E-6 |

The positive-imaginary unit-circle modes on the stable side have **opposite
Krein signatures** at both brackets; off-circle complex quartets appear on
the other side. These support candidate spectral transitions. Exact-transition
algebraic/geometric multiplicities are unresolved. Numerical endpoint signs
do not constitute interval-arithmetic enclosures or certify all intermediate
masses as spectrally stable.

An independently initialized, time-reversed seed also closes the full state.
Its full monodromy agrees with R*inverse(Phi)*R to
1.25253E-37, with R=diag(I6,-I6). This is
an additional reversibility check, not an independent MP implementation.

## Tests and saved artifacts

15 final tests run; 0 failures/errors. See `tests_all.xml` and `tests_all.log`.

See `high_precision_orbits.json`, `high_precision_orbits_85.json`,
`boundary_orbits_65.json`, `boundary_orbits_85.json`,
`floquet_sweep_high_precision.csv`, `stability_boundaries.json`,
`high_precision_convergence.json`, and `time_reversal_verification.json`.
Raw accepted-step states, velocities and variational columns are in `mp_raw/`.
The verified full-state trajectory and orbit/conservation/Floquet/boundary
figures are saved. Source snapshots and file hashes preserve provenance.
The initial reciprocal-quartic analysis is retained separately; finalized
spectra refine roots on the actual polynomial instead of imposing symmetry.

## Remaining work

1. Cross-check with an independent arbitrary-precision implementation and/or
   interval arithmetic, particularly the monodromy and transition enclosures.
   Two Taylor settings share an implementation; binary64 DOP853 is the existing
   independent method but cannot verify 28-digit accuracy by itself.
2. Resolve the exact-collision multiplicities and refine Krein/Jordan analysis
   before upgrading the candidate bifurcation label.
3. Obtain the complete original atlas draft and the 2023 article/supplements.
   The draft PDF could not be fetched; only the atlas pages and publisher
   preview were inspected. The literature/novelty audit is incomplete.
4. Only after that audit consider a manuscript or asymmetric-mass continuation
   with a newly derived full shooting/phase setup. No manuscript is prepared.

Reproduce with the root `README.md`; numerical conventions and solver limits
are documented in `methodology.md`. Literature attribution/access limits are
in `literature_audit.md`.
