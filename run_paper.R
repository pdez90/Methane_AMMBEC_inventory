# run_paper.R — every analysis behind the Paper 2 manuscript, from the priors to the results tables ------
#   cd ~/Methane_AMMBEC/methane_inversion && Rscript run_paper.R
# Assumes 01–03 (observations, HRRR, STILT footprints) have been run: see run_all.R. Everything here is
# deterministic given the footprints and the Stan seed, and reruns in ~2 h on an 8-core Mac.
#
#   04   priors (area-weighted region masks; EPA box = 2.80 t/h as in Paper 1)
#   05/07 Jacobian and real-data inversion at each mixed-layer scaling (METHANE_ZISCALE); 0.8 is the
#        central run, 0.8–1.0 the supported range, 0.37 / 0.5 / 1.2 the bracketing cases (Section 3.4)
#   05b  transport checks (particle number, lidar vs HRRR mixed-layer height)
#   05c  aircraft-profile mixed-layer height (the constraint that sets the central ZISCALE)
#   06   OSSE (realistic errors) and 06b identifiability regimes, at the central transport
#   07b  background sensitivity (nine cases) at the central transport
#   12   ethane endmember sensitivity (ten cases)
#   13   posterior predictive checks, PSIS-LOO / WAIC model comparison, OSSE coverage
#   14   leave-one-flight-out jackknife (METHANE_RUN_JACKKNIFE=FALSE to skip)
#   09   Tables 1–6 and results_summary.md
#   11   three-way split (fossil / biogenic / other) and the missing-sector test (v1.1 + waste)
#   10   Figures 1–6 (inversion_outputs/figures/paper2)
#   08   Metro Water Recovery secondary analysis (SI S9); not needed for the main tables
ZI     <- c("0.37", "0.5", "0.8", "1.0", "1.2")
ZI_MAIN <- "0.8"
Sys.setenv(METHANE_CORES = Sys.getenv("METHANE_CORES", "8"))
rs <- function(script, env = character()) {
  cat("\n==", script, if (length(env)) paste(names(env), env, sep = "=", collapse = " "), "==\n")
  a <- strsplit(script, " ")[[1]]
  st <- system2(file.path(R.home("bin"), "Rscript"), c(file.path("scripts", a[1]), a[-1]), env = paste0(names(env), "=", env))
  if (st != 0) stop(script, " failed (exit ", st, ")")
}
rs("04_priors.R")
for (zi in ZI) { e <- c(METHANE_ZISCALE = zi); rs("05_jacobian.R", e); rs("07_invert.R", e) }
rs("05b_transport_checks.R", c(METHANE_ZISCALE = ZI_MAIN))
rs("05c_blh_aircraft.R")
rs("06_osse.R gra2pes_v2.0beta 20", c(METHANE_ZISCALE = ZI_MAIN))   # 20 replicates for coverage
rs("06b_identifiability_regimes.R", c(METHANE_ZISCALE = ZI_MAIN))
rs("07b_background_sensitivity.R", c(METHANE_ZISCALE = ZI_MAIN))
rs("12_endmember_sensitivity.R", c(METHANE_ZISCALE = ZI_MAIN))
rs("13_predictive_checks.R", c(METHANE_ZISCALE = ZI_MAIN))
if (isTRUE(as.logical(Sys.getenv("METHANE_RUN_JACKKNIFE", "TRUE")))) rs("14_jackknife.R", c(METHANE_ZISCALE = ZI_MAIN))   # ~1 h
rs("09_results_table.R", c(METHANE_MAIN_ZI = ZI_MAIN, METHANE_RANGE_ZI = "0.8,1.0"))
if (isTRUE(as.logical(Sys.getenv("METHANE_RUN_WWTP", "FALSE")))) rs("08_metro_wwtp.R", c(METHANE_ZISCALE = ZI_MAIN))
# Missing-sector test: GRA2PES v1.1 given a waste sector (two magnitudes), central transport only
for (v in c("epa", "v2")) { e <- c(METHANE_V11_WASTE = v, METHANE_ZISCALE = ZI_MAIN)
  rs("04_priors.R v1.1", e); rs("05_jacobian.R", e); rs("07_invert.R", e) }
rs("11_missing_sector.R", c(METHANE_MAIN_ZI = ZI_MAIN))
rs("10_figures.R", c(METHANE_MAIN_ZI = ZI_MAIN, METHANE_RANGE_ZI = "0.8,1.0"))
cat("\nAll Paper 2 analyses complete. Tables are in inversion_outputs/results/.\n")
