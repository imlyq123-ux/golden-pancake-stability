# Manuscript draft (for author review — not submitted)

`paper.tex` / `paper.pdf`: a 4-page draft, *"Stability windows of the Golden Pancake three-body
isola: Krein collisions, a half-period reflection symmetry and outer-mass sensitivity"*.
Authors: Yuanqing Liu, Han Liu and Anqing Xuan.

Why a draft exists now: the verification and novelty criteria of `../AGENTS.md` were met in the
second session — the transition masses are reproduced by two independent multiprecision
pipelines (about 60 agreeing digits for the main window edges), the representative orbit agrees
with the atlas's own independent data, and none of the inspected sources (atlas draft and
catalogue, adjacent literature; `../results/literature_audit.md`) reports the stability
transitions. The draft claims numerical evidence only, credits the atlas for the family, and
lists its limitations (no interval proof; binary64 scans and Phase E; collision limits not
reached; spectral stability only).

Before any submission: confirm novelty with the atlas authors, read the full text of
Vasiljević, Raonić & Dmitrašinović (2023), add affiliations, and consider the open items in
`../results/research_status.md` (section 9).

Build: `tectonic paper.tex` (or any LaTeX with graphicx, booktabs, hyperref, microtype,
lmodern); figures are read from `../figures/`. The superseded first-session working note is in
`../threebody_codex_project/pancake_stability_note.*`.
