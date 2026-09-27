# record_environment.R — write environment.md: the R, CmdStan, package and platform versions used --------
#   Rscript record_environment.R        (run on the machine that produced the manuscript numbers)
# Complements renv.lock (renv::snapshot()) with a human-readable record for the paper's Data and Code
# Availability statement.
pk <- c("cmdstanr", "posterior", "ncdf4", "parallel", "renv")
ver <- function(p) tryCatch(as.character(packageVersion(p)), error = function(e) "not installed")
cs <- tryCatch(cmdstanr::cmdstan_version(), error = function(e) "cmdstan not found")
csp <- tryCatch(cmdstanr::cmdstan_path(), error = function(e) "")
stilt <- Sys.getenv("STILT_DIR", "")
lines <- c(
  "# Computational environment (Paper 2)", "",
  sprintf("Recorded %s on %s (%s).", format(Sys.time(), "%Y-%m-%d %H:%M %Z"), Sys.info()[["nodename"]], R.version$platform), "",
  sprintf("- R %s", getRversion()),
  sprintf("- CmdStan %s (%s)", cs, csp),
  sprintf("- macOS / OS: %s", paste(Sys.info()[c("sysname", "release")], collapse = " ")),
  sprintf("- Cores used for STILT/Stan (METHANE_CORES): %s", Sys.getenv("METHANE_CORES", "unset")),
  "", "## R packages", "",
  sprintf("- %s %s", pk, sapply(pk, ver)),
  "", "## Transport", "",
  "- STILT (HYSPLIT-based, Fasoli et al. 2018) run through `R/stilt.R`; 1000 particles, 10 h backward",
  "- HRRR 3 km analyses converted to ARL by `scripts/02_met.R`",
  if (nzchar(stilt)) sprintf("- STILT_DIR: %s", stilt) else "",
  "", "## Seeds", "",
  "- Stan: seed = 1 in `R/inversion.R::fit_inversion()`; 4 chains x 1000 warmup + 1000 sampling iterations, adapt_delta 0.95",
  "- OSSE replicates: seeds in `scripts/06_osse.R`",
  "", "## Reproducing the manuscript numbers", "",
  "`Rscript run_all.R` (01–03: observations, HRRR, footprints), then `Rscript run_paper.R` (04–09).",
  "`git tag` marks the code state of each manuscript version.")
writeLines(lines[nzchar(lines) | lines == ""], "environment.md")
cat(paste(lines, collapse = "\n"), "\n\nwrote environment.md\n")
