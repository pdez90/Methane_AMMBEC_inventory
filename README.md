# methane_inversion — Paper 2: joint CH₄–C₂H₆–CO source inversion for Denver and the DJ Basin

A low-dimensional Bayesian inversion of the 2024 AMMBEC Twin Otter flights. HRRR-driven STILT footprints
map sector-resolved prior emissions (GRA2PES v1.1, v2.0beta, EPA GHGI) onto every level-leg segment.
A Stan model then estimates one scale factor per source component, uses ethane to separate fossil
from biogenic methane, and uses CO to constrain transport.

This repo sits next to Paper 1's `methane_AMMBEC/` and never modifies it. It sources Paper 1's ICARTT
reader, leg detection and enhancement code from `../methane_AMMBEC/R/`.

## Pipeline

| step | script | what it does | where it runs |
|---|---|---|---|
| 0 | `scripts/00_setup.R` | installs cmdstanr + CmdStan, fetches STILT binaries, checks inputs | Mac |
| 1 | `scripts/01_observations.R` | all 22 flights, all level legs in Denver + DJB below 3 km, cut into ~60 s segments (CH₄, C₂H₆, CO) → `obs_segments.csv` | anywhere |
| 2 | `scripts/02_met.R` | HRRR ARL files the footprints need (46; the 15 from Paper 1 are reused), download + cut to Colorado | Mac (needs NOAA) |
| 3 | `scripts/03_footprints.R` | one STILT receptor per segment, 300 particles, 10 h back → `footprints/<id>/foot.rds` | Mac (local engine) or any Linux box (`run_stilt_linux.sh`) |
| 4 | `scripts/04_priors.R` | GRA2PES v1.1 / v2.0beta (hourly, 4 km) and EPA (0.1°) split into inversion components; CO fields | Mac (needs the GRA2PES trees) |
| 5 | `scripts/05_jacobian.R` | footprints × priors → H (ppb per component) and prior CO; forward-model check vs observations | anywhere |
| 6 | `scripts/06_osse.R` | synthetic-data tests: what these flights can and cannot identify, with CH₄ / +C₂H₆ / +CO | anywhere |
| 6b | `scripts/06b_identifiability_regimes.R` | synthetic experiment on what controls identifiability: ethane contrast between basin and urban gas × ethane noise × sampling density, CH₄ vs CH₄ + C₂H₆ | anywhere |
| 7 | `scripts/07_invert.R` | the real inversion, all priors × tracer sets (run after reading the OSSE) | anywhere |

`Rscript run_all.R` runs 01–06 in order. Every step skips work that is already done.

## Components (config.R)

`dads`, `tower_road`, `metro_complex`: the three GRA2PES v2.0beta waste hotspot cells from Paper 1. They get
their own factors because Paper 1 found their inventory errors have opposite signs. The other components are
`waste` (the rest), `og_basin` / `og_urban` (GRA2PES OG split by the Paper 1 box + 8 km), `postmeter` (RES),
`ag` and `other` (total minus sectors).

Ethane ratios have priors from Paper 1: basin 0.080 ± 0.015, urban gas 0.110 ± 0.012, other 0.02 ± 0.02.
Waste and agriculture are fixed at 0.

## Model (stan/inversion.stan)

For segment *i* on flight *f*:

- CH₄ = b_CH₄,f + τ_f Σ_k α_k H_ik + leg offset + ε
- C₂H₆ = b_C₂H₆,f + τ_f Σ_k α_k r_k H_ik + leg offset + ε
- CO = b_CO,f + τ_f γ Hco_i + leg offset + ε

The terms:

- α_k: component scale factors, lognormal(0, 1)
- r_k: ethane ratios
- τ_f: flight transport factor, lognormal(0, 0.5), shared by all three species
- γ: CO inventory scale, lognormal(0, 0.15)
- b_f: flight backgrounds
- ε: error with sd √(s² + (φ·enhancement)²)

CO pins τ_f, which CH₄ alone cannot separate from α. Posterior box emissions are α_k × the prior box totals,
and the fossil share is (og_basin + og_urban + postmeter) / total.

## Status (24 Sep 2026)

- **Tested:**
  - 01: 634 segments, 22 flights: 115 urban, 148 edge, 371 basin.
  - 04 on real GRA2PES subsets: box totals reproduce Paper 1 exactly (v1.1 1.69 t/h; v2.0beta 16.54 t/h, waste 12.20, OG 3.67).
  - 03 end to end on 6 real receptors: STILT's Linux build with real HRRR on the Mac.
  - 05, 06, 07 on Paper 1's 100 urban receptors.
- **First OSSE (100 urban receptors only, 2 replicates):**
  - The landfill, Metro-complex and basin-OG factors are well constrained.
  - Ethane raises basin-OG error reduction from 0.58 to 0.72 and urban-OG from 0.30 to 0.48.
  - Post-meter and "other" are not identifiable.
  - Urban OG is biased low. Its diffuse signal trades off against the flight backgrounds, whose prior is
    centred on each flight's 5th percentile. With urban plume receptors only, that percentile already
    contains urban methane. Check this with the full 634-segment set, which includes basin, edge and
    upwind legs; if it persists, centre the background prior on upwind or free-troposphere segments.
- **Not yet run:**
  - 02 for the other 31 HRRR files (~105 GB of downloads, ~5 h).
  - 03 for all 634 segments.
  - The OSSE on the full set.
  - No real-data result.

## Known limitations to address

- GRA2PES sector files are July 2023 weekdays; the flights are July 2024. CO uses Saturday/Sunday fields where they exist (v1.1).
- A 10-h footprint within DOMAIN. Older or out-of-domain air is folded into the backgrounds.
- Priors for backgrounds and error scales are in `R/inversion.R::inv_settings()`; test them with the OSSE before the real run.

## Identifiability (added 25 Sep 2026)

The paper's question is framed as "what can CH₄ alone identify, and what does ethane add?", not only "what are
Denver's emissions?". `R/identifiability.R` computes these from every fit (06, 06b, 07), in log α space:

- error reduction and averaging-kernel diagonal per component;
- DOFS = K − tr(S_prior⁻¹ S_post);
- information gain H(prior) − H(posterior) in bits, and dIG = the gain from each tracer relative to CH₄ alone;
- posterior correlation between the fossil and biogenic box totals (strongly negative means the total is known but
  the split is not), the full component correlation matrix, and the sd of the fossil share;
- a regime for each component: well identified (error reduction ≥ 0.5 and |corr| < 0.5 with every other component),
  partially identified, or not identified (error reduction < 0.2).

07 also writes `posterior_prior_dependence.csv`: how much of the ~10-fold spread between the three inventory priors
survives in the posterior (0 = the data decide, 1 = the prior decides).
