# Three-Body Golden Pancake Stability — Codex starter project

Open this directory in Codex and ask: **“Read AGENTS.md and execute Phase A. Do not claim a research result until reproducibility and precision have been verified.”**

- **AGENTS.md**: full research brief, previous preliminary results, known validation gaps, prioritized implementation plan and acceptance criteria.
- **threebody_expand.py**: executable original exploratory double-precision shooting/variational code.
- **pancake_stability_note.pdf**: preliminary 3-page technical note (not for submission).
- Other `pancake_*` files: original numerical data, scripts, plots and failure log.
- `golden_pancake_mu_0p613_*`: earlier single-orbit trajectory and plot; not a complete multiprecision reproducibility record.

The initial code is intentionally preserved as supplied; some helper scripts are nonportable and must be refactored. See AGENTS.md before using reported digits in a paper.

The central scientific objective is to confirm or refute two tentative nontrivial spectral-stability boundaries on a previously published orbit family—not to claim the discovery of a new family.
