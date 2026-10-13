EXPLORATORY GOLDEN PANCAKE STABILITY STUDY — October 12, 2026

Purpose: reproducible double-precision exploratory continuation and
Floquet stability analysis for the already-known Golden Pancake Drama
unequal-mass isola. NOT a new-family discovery or arXiv-ready result.

Main findings:
  Approximate boundary masses: mu_lower ~ 0.5990438, mu_upper ~ 0.6178443.
  These are floating-point numerical estimates, not validated digits.
  23 nearby orbit solutions span mu=0.585 to 0.630.

FILE GUIDE
pancake_stability_note.pdf: 3-page technical report
pancake_stability_note.tex: LaTeX source for note
threebody_expand.py: Newton continuation and analytic full monodromy
pancake_stability_boundaries.py: refined spectral discriminant boundary search
pancake_upper_refine.py: upper crossing numerical refinement
pancake_convergence_check.py: tolerance sensitivity checks
make_pancake_results.py: assemble mass-sweep files and plots
pancake_mass_sweep.csv: 23 numerical samples (mu,u,v,T,E,diagnostics,multipliers)
pancake_conservation_tests.csv: five sampled energy and angular momentum checks
pancake_stability_sweep.png: largest nontrivial Floquet modulus
pancake_discriminant_sweep.png: spectral discriminant
pancake_nearest_approach.png: minimum pair distance
pancake_continuation_*.json: numerical continuation raw results
pancake_transition_coarse.json: Floquet samples near upper crossing
pancake_upper_boundary.json: detailed upper transition refinement
pancake_boundaries.json: partial lower and upper numerical bracketing

REPRODUCIBILITY
Requires Python 3.11+, numpy, scipy, numba, pandas and matplotlib.
From the unpacked folder (after Numba compilation on first run):
  OPENBLAS_NUM_THREADS=1 python threebody_expand.py \
    --mass-list 0.613,0.612,0.611,0.610,0.608,0.605 --stability \
    --out checked_branch.json
Each continuation call starts from the previous parameter root; do not
jump by large mass increments without intermediate values.

MODEL
Masses: (1,mu,1), gravitational constant G=1.
Initial q1=(-1,0),q2=(0,0),q3=(1,0),
v1=v3=(u,v), v2=(-2u/mu,-2v/mu).

PRECISION AND LIMITATIONS
This EXPANDED sweep is double precision with full phase-space return
residuals approximately 4e-12 to 2.2e-11. The previous high-precision
single mu=0.613 point from the conversation is separate; it does not
verify the full stability sweep or the transition masses.

At Hamiltonian-Hopf threshold, the nontrivial spectral quartet merges;
root locations should be re-computed with 50-100 digit integrations
and symmetry-reduced symplectic monodromy. Original atlas research note
full text was not accessible for a content-complete novelty comparison.

Original published family:
https://www.threebodyorbits.com/family/golden-pancake-drama
Atlas methodology:
https://www.threebodyorbits.com/about
Original research note (linked but not accessed as PDF):
https://www.threebodyorbits.com/papers/golden-pancake-drama-2026.pdf
