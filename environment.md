# Computational environment (Paper 2)

Recorded 2026-09-26 18:21 MDT on Priyankas-MBP-2 (aarch64-apple-darwin20).

- R 4.5.0
- CmdStan cmdstan not found ()
- macOS / OS: Darwin 27.0.0
- Cores used for STILT/Stan (METHANE_CORES): unset

## R packages

- cmdstanr not installed
- posterior 1.6.1
- ncdf4 1.24
- parallel 4.5.0
- renv not installed

## Transport

- STILT (HYSPLIT-based, Fasoli et al. 2018) run through `R/stilt.R`; 1000 particles, 10 h backward
- HRRR 3 km analyses converted to ARL by `scripts/02_met.R`


## Seeds

- Stan: seed = 1 in `R/inversion.R::fit_inversion()`; 4 chains x 1000 warmup + 1000 sampling iterations, adapt_delta 0.95
- OSSE replicates: seeds in `scripts/06_osse.R`

## Reproducing the manuscript numbers

`Rscript run_all.R` (01–03: observations, HRRR, footprints), then `Rscript run_paper.R` (04–09).
`git tag` marks the code state of each manuscript version.
