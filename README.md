# methane_inversion — Paper 2: joint CH₄–C₂H₆–CO source inversion for Denver and the DJ Basin

Code for *Which methane sources are identifiable from airborne observations? An ethane-constrained Bayesian
inversion over Denver and the Denver–Julesburg Basin* (submitted to ES&T Air). A low-dimensional Bayesian
inversion of the 22 AMMBEC Twin Otter flights of June–July 2024. HRRR-driven STILT footprints map
sector-resolved prior emissions (EPA GHGI 2020, GRA2PES v1.1, GRA2PES v2.0β) onto 634 level-leg segments. A
Stan model estimates one scale factor per source component, uses ethane to separate fossil from biogenic
methane, and (as a secondary test) uses CO to constrain transport. The paper's contribution is the
identifiability analysis around that inversion: what the observations determine, what ethane adds, and what
remains controlled by the prior or by transport assumptions.

This repo sits next to Paper 1's `methane_AMMBEC/` and never modifies it. It sources Paper 1's ICARTT
reader, leg detection and enhancement code from `../methane_AMMBEC/R/`.

## Reproducing the paper

```
Rscript run_all.R      # 01–03: observations, HRRR download, STILT footprints (hours to days; skips finished work)
Rscript run_paper.R    # 04–09: priors, Jacobians, inversions at five mixed-layer scalings, OSSEs,
                       #        transport and background sensitivities, aircraft MLH, results tables (~2 h)
```

`run_paper.R` writes `inversion_outputs/results/table1_main_results.csv … table6_prior_dependence.csv` and
`results_summary.md`; every number in the manuscript's Tables 1–6 comes from those files. Inputs (ICARTT
files, HRRR, GRA2PES and EPA grids) live outside the repo; paths are in `config.R`. Environment:
`environment.md` (R, CmdStan and package versions of the final run) and `renv.lock`.

Git tags mark the code state behind each manuscript version (`v0.3-paper2`, …).

## Pipeline

| step | script | what it does |
|---|---|---|
| 0 | `scripts/00_setup.R` | installs cmdstanr + CmdStan, fetches STILT binaries, checks inputs |
| 1 | `scripts/01_observations.R` | all 22 flights, all level legs in Denver + DJB below 3 km, cut into ≈60 s segments → `obs_segments.csv` (634: 115 urban, 148 edge, 371 basin) |
| 2 | `scripts/02_met.R` | HRRR ARL files (46), download + cut to Colorado |
| 3 | `scripts/03_footprints.R` | one STILT receptor per segment, 1000 particles, 10 h back → `footprints/<id>/foot.rds` |
| 4 | `scripts/04_priors.R` | EPA GHGI / GRA2PES v1.1 / v2.0β split into the inversion components; region totals with area-weighted masks (EPA box = 2.80 t/h, as in Paper 1) |
| 5 | `scripts/05_jacobian.R` | footprints × priors → H (ppb per component) at the chosen `METHANE_ZISCALE`; forward-model check |
| 5b | `scripts/05b_transport_checks.R` | particle-number test; lidar vs HRRR mixed-layer height |
| 5c | `scripts/05c_blh_aircraft.R` | mixed-layer height from the aircraft's own profiles (Heffter criterion, censoring); the constraint that sets the central ZISCALE = 0.8 |
| 6 | `scripts/06_osse.R` | OSSE with realistic errors: coverage of totals and fossil shares, CH₄ vs +C₂H₆ vs +CO |
| 6b | `scripts/06b_identifiability_regimes.R` | synthetic experiment: ethane contrast × noise × sampling density |
| 7 | `scripts/07_invert.R` | the real inversion, all priors × tracer sets, at one transport setting; identifiability and prior dependence |
| 7b | `scripts/07b_background_sensitivity.R` | nine background treatments at the central transport |
| 8 | `scripts/08_metro_wwtp.R` | Metro Water Recovery secondary analysis (SI S9); `METHANE_METRO_COMPONENT=wwtp` prior variant |
| 9 | `scripts/09_results_table.R` | Tables 1–6 and `results_summary.md` from all runs |

Transport sweep: `METHANE_ZISCALE` ∈ {0.37, 0.5, 0.8, 1.0, 1.2}, each a full 05 → 07 chain in
`inversion_outputs/runs/np1000_h10[_ziX.XX]/`. 0.8 is the central run (aircraft profiles imply ≥ 0.78),
0.8–1.0 the supported range, 0.37 and 1.2 bracketing failure cases (Section 3.4 of the paper).
`METHANE_CORES=8` avoids a macOS fork error in the STILT step.

## Components (config.R)

Up to nine per prior: `dads`, `tower_road` (landfills, 3 km), `metro_complex` (Metro Water Recovery / Suncor,
5 km), `waste` (the rest), `og_basin` / `og_urban` (oil and gas split by the Paper 1 box + 8 km), `postmeter`
(residential gas), `ag`, `other` (total minus sectors). Components absent from a prior (v1.1 has no waste,
landfill or wastewater methane inside the box) or unseen by the footprints are dropped from that fit.

Ethane ratios (normal priors truncated at zero) from Paper 1: basin 0.080 ± 0.015, urban gas 0.110 ± 0.012,
other 0.02 ± 0.02; waste and agriculture fixed at 0. The fossil share is defined as
(og_basin + og_urban + postmeter) / total (`FOSSIL_COMPONENTS` in `R/inversion.R`).

## Model (stan/inversion.stan; settings in R/inversion.R::inv_settings, "v2")

For segment *i* on flight *f*, leg *l*:

- CH₄ = b_f + τ_f Σ_k α_k H_ik + δ_l + ε
- C₂H₆ = b'_f + τ_f Σ_k α_k r_k H_ik + δ'_l + ε'
- CO = b''_f + τ_f γ Hco_i + δ''_l + ε''

α_k ~ lognormal(0, 1); τ_f ~ lognormal(0, 0.5), shared by all species; γ ~ lognormal(0, 0.15); backgrounds
b_f centred on the flight's 5th percentile with sd 5 / 0.3 / 5 ppb; leg offsets δ_l with scale 3 / 0.2 / 3 ppb.
Errors are Student-t (ν = 4) with scale √(s² + (φ·ŷ_i)²), where ŷ_i = τ_f Σ_k α_k H_ik is the **modeled**
emission enhancement at segment *i* and φ is estimated. Model "v1" (Gaussian errors, 20 ppb background sd) is
kept for reference (`METHANE_INV_MODEL=v1`); it assigned about half the urban enhancement to backgrounds.

## Identifiability (R/identifiability.R)

From every fit (06, 06b, 07), in log α space, using the posterior covariance of the MCMC draws as a
second-moment Gaussian approximation:

- error reduction and averaging-kernel diagonal per component;
- DOFS = K − tr(S_prior⁻¹ S_post);
- information gain H(prior) − H(posterior) in bits, and its change when a tracer is added;
- posterior correlation between the fossil and biogenic totals, the component correlation matrix, the sd of the
  fossil share;
- a regime per component: well identified (error reduction ≥ 0.5 **and** |corr| < 0.5 with every other
  component), not identified (error reduction < 0.2), partially identified otherwise.

07 writes `posterior_prior_dependence_v2.csv` and 09 writes `table6_prior_dependence.csv`: D = across-prior
spread of posterior medians / spread of the priors (0 = data-determined, 1 = prior-determined), paired with
R = across-prior spread / posterior 1σ.

## Status (26 Sep 2026)

Final run for manuscript v0.4: all 634 receptors, five transport settings, three priors, three tracer sets,
OSSE (06, 06b), nine background cases (07b), aircraft MLH (05c), results tables (09). All real-data fits:
0 divergent transitions, R̂ ≤ 1.02 (the ZISCALE 0.37 bracketing case reports low E-BFMI in one or two chains;
it is a documented failure mode, not a headline run). Headline checks: EPA box prior 2.80 t/h reproduces
Paper 1; 0.8 and 1.0 runs reproduce bit-for-bit on rerun; box total prior dependence D = 0.25, box fossil-share
D = 1.1.

## Known limitations

- GRA2PES sector files are July 2023 weekdays; the flights are July 2024. CO uses Saturday/Sunday fields
  where they exist (v1.1).
- A 10-h footprint within DOMAIN; older or out-of-domain air is folded into the backgrounds.
- The OSSE shows the posterior intervals on *totals* are under-covered by about a factor of two; the paper
  reports ensemble ranges (priors × transport range × background) for totals, and posterior intervals only
  for fossil shares, which are calibrated.
