# 00_setup.R — one-time checks and installs for the Paper 2 pipeline ------------------------------
#   Rscript scripts/00_setup.R
# 1. R packages: ncdf4 (CRAN), cmdstanr (Stan's r-universe), then CmdStan itself (~10 min compile).
#    macOS needs the Xcode command-line tools for CmdStan (xcode-select --install).
# 2. STILT's hycs_std and xtrct_grid (uataq/stilt bin/macos_x64) into <Methane_AMMBEC>/stilt/bin,
#    with the Intel gfortran runtime beside them (R/stilt_install.R). Apple Silicon needs Rosetta 2.
# ----------------------------------------------------------------------------------------------
proj <- if (file.exists("config.R")) "." else ".."; source(file.path(proj, "config.R"))
ok <- function(x, msg) cat(sprintf("  [%s] %s\n", if (x) " ok " else "MISS", msg))

cat("R packages\n")
if (!requireNamespace("ncdf4", quietly = TRUE)) install.packages("ncdf4", lib = R_LIB, repos = "https://cloud.r-project.org")
ok(requireNamespace("ncdf4", quietly = TRUE), "ncdf4")
if (!requireNamespace("cmdstanr", quietly = TRUE))
  install.packages("cmdstanr", lib = R_LIB, repos = c("https://stan-dev.r-universe.dev", "https://cloud.r-project.org"))
ok(requireNamespace("cmdstanr", quietly = TRUE), "cmdstanr")
if (requireNamespace("cmdstanr", quietly = TRUE)) {
  invisible(use_local_cmdstan())
  v <- tryCatch(cmdstanr::cmdstan_version(error_on_NA = FALSE), error = function(e) NULL)
  if (is.null(v)) {
    cat("  installing CmdStan into ", CMDSTAN_DIR, " (one time, ~10 min)...\n", sep = "")
    dir.create(CMDSTAN_DIR, FALSE, TRUE)
    cmdstanr::install_cmdstan(dir = CMDSTAN_DIR, cores = max(1, N_CORES))
    invisible(use_local_cmdstan())
  }
  v <- tryCatch(cmdstanr::cmdstan_version(error_on_NA = FALSE), error = function(e) NULL)
  ok(!is.null(v), paste("CmdStan", if (is.null(v)) "" else v))
}

cat("STILT\n")
source(file.path(proj, "R", "stilt_install.R"))
if (Sys.info()[["sysname"]] == "Darwin") {
  if (!rosetta_ok()) {
    ok(FALSE, "Rosetta 2 (STILT's macOS binaries are Intel-only)")
    cat("    run this once in Terminal (asks for your Mac password), then run 00_setup.R again:\n",
        "      softwareupdate --install-rosetta --agree-to-license\n")
    quit(save = "no", status = 1)
  }
  ok(TRUE, "Rosetta 2")
  install_stilt_mac(dirname(STILT_EXE))
}
for (e in c(STILT_EXE, XTRCT_EXE)) {
  s <- if (file.exists(e)) stilt_starts(e) else FALSE
  ok(isTRUE(as.logical(s)), paste(basename(e), if (isTRUE(as.logical(s))) "runs" else paste("does not run:", attr(s, "log"))))
}
cat("Data\n")
ok(dir.exists(file.path(DATA_DIR, "Aircraft")), file.path(DATA_DIR, "Aircraft"))
ok(dir.exists(file.path(GRA2PES_DIR, "sectors")), file.path(GRA2PES_DIR, "sectors"))
ok(file.exists(GHGI_FILE), GHGI_FILE)
ok(file.exists(file.path(STILT_BDY, "LANDUSE.ASC")), file.path(STILT_BDY, "LANDUSE.ASC"))
cat("outputs go to ", INV_OUT, "\n", sep = "")
