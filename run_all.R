# run_all.R — the Paper 2 pipeline in order --------------------------------------------------------
#   cd ~/Methane_AMMBEC/methane_inversion && Rscript run_all.R [from_step]
# Steps 02 (downloads) and 03 (STILT) are long; each skips work already done, so re-running is cheap.
steps <- c("01_observations.R", "02_met.R", "03_footprints.R", "04_priors.R", "05_jacobian.R", "06_osse.R", "06b_identifiability_regimes.R")
# 07_invert.R (real-data inversion) is deliberately not in the default run: read the OSSE first.
from <- as.integer(commandArgs(TRUE)[1]); if (is.na(from)) from <- 1L
for (s in steps[as.integer(substr(steps, 1, 2)) >= from]) {
  cat("\n==", s, "==\n")
  st <- system2(file.path(R.home("bin"), "Rscript"), file.path("scripts", s))
  if (st != 0) stop(s, " failed (exit ", st, ")")
}
