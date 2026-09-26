# run_mac.R — run the Paper 2 pipeline steps on the Mac, in R, logged to ../logs/ ----------------
#   cd ~/Methane_AMMBEC/methane_inversion
#   Rscript run_mac.R 00 01 02 03 04 05 06        # any subset, in order; stops at the first failure
# Each step writes ../logs/<step>.log (full output) and ../logs/<step>.status (RUNNING / OK / FAIL).
# Inputs are read from /Volumes/Elements/Methane_AMMBEC; everything is written under
# /Users/priyanka/Methane_AMMBEC. 07_invert.R is left out on purpose (read the OSSE first).
# ----------------------------------------------------------------------------------------------
if (!file.exists("config.R")) stop("run this from the methane_inversion folder")
if (!dir.exists("/Volumes/Elements/Methane_AMMBEC")) stop("the Elements drive is not mounted")
steps <- commandArgs(TRUE); if (!length(steps)) steps <- sprintf("%02d", 0:6)
logs <- normalizePath(file.path("..", "logs"), mustWork = FALSE); dir.create(logs, FALSE, TRUE)
# keep the Mac awake while this R session runs (macOS caffeinate, released when R exits)
if (Sys.info()[["sysname"]] == "Darwin")
  system2("caffeinate", c("-i", "-w", Sys.getpid()), wait = FALSE)
rscript <- file.path(R.home("bin"), "Rscript")
for (s in steps) {
  f <- list.files("scripts", pattern = paste0("^", s, "_.*\\.R$"), full.names = TRUE)[1]
  if (is.na(f)) stop("no script for step ", s)
  lg <- file.path(logs, paste0(s, ".log")); stf <- file.path(logs, paste0(s, ".status"))
  writeLines("RUNNING", stf); t0 <- Sys.time()
  cat(sprintf("[%s] %s ... (log: %s)\n", format(t0, "%H:%M"), basename(f), lg))
  rc <- system2(rscript, f, stdout = lg, stderr = lg)   # the step's whole output goes to its log
  rc <- if (is.null(rc)) 0L else rc
  cat(paste("== exit", rc, "started", format(t0), "ended", format(Sys.time())), file = lg, sep = "\n", append = TRUE)
  mins <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
  if (rc == 0) {
    writeLines("OK", stf); cat(sprintf("    ok (%.1f min)\n", mins))
    tl <- tail(readLines(lg, warn = FALSE), 4); cat(paste0("    ", tl[-length(tl)]), sep = "\n")
  } else {
    writeLines(paste("FAIL", rc), stf)
    cat(sprintf("    FAILED (exit %d after %.1f min). Last lines of the log:\n", rc, mins))
    cat(paste0("    ", tail(readLines(lg, warn = FALSE), 15)), sep = "\n")
    quit(save = "no", status = 1)
  }
}
cat("done:", paste(steps, collapse = " "), "\n")
