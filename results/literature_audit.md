# Literature and novelty audit

Updated 2026-10-12 (second session). The first session could not download the atlas draft
and recorded the audit as incomplete. In this session the draft PDF, the atlas catalogue
and the atlas orbit files were downloaded directly (SHA256 values and the transcribed
numbers are in [`literature/atlas_reference_data.json`](literature/atlas_reference_data.json)).
No paywalled full text was obtained for the 2023 *New Astronomy* article; its metadata and
reference list (Crossref) and its public data repository were inspected.

## 1. The family and what its discoverers report

**Three Body Orbits (2026), "golden pancake drama: A linearly stable three-body orbit on an
unequal-mass isola", draft numerical research note, 9 September 2026**
(<https://www.threebodyorbits.com/papers/golden-pancake-drama-2026.pdf>, 3 pages).

* Uses exactly our model and Euler-phase section (their Eq. 2), masses `(1, mu, 1)`.
* Reports one representative at `mu = 0.612356088827994835162485140041506525569967` with
  42-digit `(u, v, T)` (their Table 1), full labeled closure `1.52e-36` from 160-bit Taylor
  arithmetic, and four nontrivial multipliers to 19 decimals (their Eq. 4).
* Pseudo-arclength continuation gives a closed isola in the Euler-half-twist, `L = 0`,
  `(1, mu, 1)` sector with `mu_min ~ 0.4689179813307963`, `mu_max ~ 0.6346042508062549`.
  One segment has the primitive word `AABAbbaabaBB` (the segment containing the stable
  orbit — our branch), the other `(ABab)^3`; they join through binary-collision limits near
  `mu ~ 0.50655` and `0.62944`.
* **Stability scope:** the caption of their Fig. 1 states *"Only the marked orbit's linear
  stability is asserted."* The note gives no stability interval, transition masses, Krein
  signatures or bifurcation types along the isola.

**Atlas catalogue** (<https://data.threebodyorbits.com/catalogue.json>, generated stamp
1790531492) lists seven `pancake` representatives with summary data. The field `stab`
equals the largest nontrivial multiplier modulus (we reproduce it, see below):

| key | mu | T0 | stab | segment (our identification) |
|---|---|---|---|---|
| pancake2026_isola_300 | 0.5876435519022455 | 8.920021894707649 | 1.1637 | AABAbbaabaBB (our branch) |
| pancake2026_golden_pancake_drama | 0.6123560888279949 | 8.845713110998389 | 1.0 | AABAbbaabaBB (our branch) |
| pancake2026_isola_020 | 0.630314929124448 | 9.157983596689066 | 2.3815 | AABAbbaabaBB (our branch) |
| pancake2026_isola_100 | 0.6041898278825628 | 10.898987931782164 | 23.44 | (ABab)^3 |
| pancake2026_isola_120 | 0.5603560388149639 | 11.449696087028208 | 8.1655 | (ABab)^3 |
| pancake2026_isola_140 | 0.5052499206446759 | 11.443137022249925 | 2.5586 | (ABab)^3 |
| pancake2026_isola_164 | 0.4689192990913388 | 10.51730512875858 | 1.3792 | (ABab)^3 |

The three AABAbbaabaBB samples straddle our stability window (unstable, stable, unstable)
but do not locate any boundary. The atlas "About" page lists research notes up to
27 September 2026; none concerns the Golden Pancake stability window.

## 2. External numerical agreement (independent computations)

Our independent GBS/mpmath pipeline (`independent_mp.py`) was continued from our own
branch, never seeded from atlas values, and compared afterwards
([`independent_crosscheck.md`](independent_crosscheck.md)):

* Draft Table 1 at its 42-digit mass: `|du| = 3.3e-38`, `|dv| = 1.9e-38`, `|dT| = 1.5e-36`
  (the draft's own closure is `1.5e-36`); multipliers agree to `4.3e-20` (19 decimals
  printed); `E` and `T*` agree to all printed digits.
* Catalogue `stab`: ours `1.1636633` vs `1.1637`, `1.0` vs `1.0`, `2.3814702` vs `2.3815`;
  energies and minimum distances agree to the printed digits; published binary64 initial
  conditions agree to their stated closure (`closure0` up to `1.3e-11`).

This places the atlas orbit on the branch studied here and gives an external check of
the integrators and of the Floquet computation at one stable point.

## 3. Adjacent literature (prior work to cite; none treats this family)

* Vasiljević, Raonić & Dmitrašinović, *An island of linearly stable non-hierarchical unequal
  mass periodic three-body orbits — A harbinger of circumbinary exoplanets?*, New Astron.
  **100**, 101969 (2023), doi:10.1016/j.newast.2022.101969. Metadata/reference list
  checked via Crossref; data and analysis notebooks at
  <https://github.com/trsavi/Three-body-orbits-correlations>. It is a statistical analysis of
  the stability classes already reported in published catalogues (zero angular momentum,
  symmetric collinear configuration, central mass `m3`), finding a compact stable domain in
  the `(L_f, m3/m)` plane. It predates the 2026 Golden Pancake family and does not
  continue individual orbits in mass to locate transitions. Its stability terminology
  (elliptic / marginal / hyperbolic / loxodromic) is compatible with ours. **Full text not
  read** (paywalled); conclusions above rest on the abstract, README and metadata.
* Li, Jing & Liao, PASJ **70**, 64 (2018) and Li, Li & Liao, Sci. China PMA **64**, 219511
  (2021): unequal-mass catalogues in the same symmetric collinear sector, including mass
  continuation of a large stable family. The atlas draft reports that Golden Pancake has no
  counterpart in these catalogues.
* Galán, Muñoz-Almaraz, Freire, Doedel & Vanderbauwhede, PRL **88**, 241101 (2002):
  stability and bifurcations of the figure-eight under mass continuation (fold connecting
  the eight with Simó's H orbit). Methodological precedent for stability changes along
  mass continuation.
* Simó, Contemp. Math. **292** (2002); Kapela & Simó, Nonlinearity **20**, 1241 (2007):
  dynamical properties and computer-assisted proof of linear stability of the figure-eight —
  precedent for rigorous (interval) stability statements, which this project does not make.
* Dmitrašinović, Hudomal, Shibayama & Sugita, J. Phys. A **51**, 315101 (2018); Hristov,
  Hristova & Tanikawa, New Astron. **125** (2026): linear stability surveys of zero
  angular-momentum orbits (equal or fixed masses).
* Krein-collision / Hamiltonian–Hopf background: Howard & MacKay, J. Math. Phys. **28**,
  1036 (1987); van der Meer, *The Hamiltonian Hopf Bifurcation*, LNM 1160 (1985).

Searches (web, 2026-10-12) for "golden pancake" stability and for Krein / Hamiltonian–Hopf
transitions on unequal-mass three-body periodic orbits found no study of this family's
stability window.

## 4. Novelty assessment

* **Not new:** the orbit family, its isola, its topology, the stable representative and its
  multipliers (atlas, 2026); the general approach of tracking linear stability along mass
  continuation (e.g. Galán et al. 2002; Li, Li & Liao 2021).
* **Not found in any inspected source:** the location of the stability transitions along
  the Golden Pancake isola, their type (Krein collisions of opposite-signature pairs with a
  square-root splitting, i.e. the linear signature of a Hamiltonian–Hopf bifurcation), the
  additional stable stretches found in this session, high-precision values for the
  transition masses, and the outer-mass (1+delta, mu, 1-delta) stability region.
* **Known in part:** the atlas computes the family in an "Euler-half-twist symmetric sector"
  with half-period, symmetry-constrained shooting, so a half-period symmetry of these orbits is
  implicitly known. We identify it concretely as X(T/2) = sigma X(0) with a reflection sigma
  and work out its spectral consequences (M = N^2; multipliers touch but cannot cross -1 within
  the family -- a spectral statement, not an exclusion of nonlinear bifurcations);
  we do not claim the symmetry itself as new.
* **Limits of this assessment:** worldwide priority cannot be proven; the 2023 article was
  not read in full; the atlas authors may hold unpublished continuation data (their note
  says the project archive retains continuation paths). Any manuscript must credit the
  atlas for the family, cite the draft, and phrase the contribution as a quantitative
  stability delimitation of a known family.
