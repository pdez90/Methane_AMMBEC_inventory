# 03_footprints.R — STILT footprints for every observation segment ------------------------------
# One receptor per segment (mid-point time, position and height above ground), NUMPAR particles,
# N_HOURS back, HRRR subgrids from 02. Each run directory is <FOOT_DIR>/<segment id>/ and ends with
# foot.rds (particle positions and footprint weights, foot > 0 only). Re-running skips receptors
# that already have foot.rds, so it can be stopped and resumed.
#
# Engines (config.R):
#   METHANE_STILT_ENGINE=local     run METHANE_STILT_EXE here, N_CORES at a time (default)
#   METHANE_STILT_ENGINE=external  only write the run directories; run them elsewhere with
#                                  run_stilt_linux.sh (any Linux box or VM that sees the folder),
#                                  then run this script again to collect foot.rds
#
#   Rscript scripts/03_footprints.R [--only <regex on segment id>] [--collect]
#   Out: <FOOT_DIR>/<id>/foot.rds, <TRANSPORT_DIR>/footprint_status.csv
# ----------------------------------------------------------------------------------------------
proj <- if (file.exists("config.R")) "." else ".."; source(file.path(proj, "config.R"))
source(file.path(proj, "R", "stilt.R"))
suppressPackageStartupMessages(library(parallel))

args <- commandArgs(TRUE)
only <- if ("--only" %in% args) args[which(args == "--only") + 1] else NULL
collect_only <- "--collect" %in% args

S <- read.csv(file.path(INV_OUT, "obs_segments.csv"), stringsAsFactors = FALSE)
if (!is.null(only)) S <- S[grepl(only, S$id), ]
R <- data.frame(id = S$id, run_time = as.POSIXct(round(as.numeric(as.POSIXct(S$time_utc, tz = "UTC")) / 60) * 60,
                                                    origin = "1970-01-01", tz = "UTC"),
                lati = round(S$lat, 5), long = round(S$lon, 5), zagl = pmax(5, round(S$agl_m)), stringsAsFactors = FALSE)
dir.create(FOOT_DIR, FALSE, TRUE)

miss <- setdiff(unique(unlist(lapply(R$run_time, hrrr_files_for))), list.files(file.path(MET_DIR, "sub")))
if (length(miss) && !collect_only && STILT_ENGINE == "external")
  message("note: ", length(miss), " HRRR subgrid(s) not visible here; fine if the runs happen on another machine")
if (length(miss) && !collect_only && STILT_ENGINE == "local") stop(length(miss), " HRRR subgrid(s) missing (run scripts/02_met.R): ", paste(head(miss, 5), collapse = " "))

one <- function(i) {
  r <- R[i, ]; wd <- file.path(FOOT_DIR, r$id)
  tryCatch({
    if (file.exists(file.path(wd, "foot.rds"))) return(file.path(wd, "foot.rds"))
    if (collect_only || STILT_ENGINE == "external") {
      if (!file.exists(file.path(wd, "CONTROL"))) stilt_prepare(r, wd)
      return(stilt_collect(r, wd))
    }
    stilt_run_local(r, wd)
  }, error = function(e) paste("ERROR", r$id, conditionMessage(e)))
}
if (STILT_ENGINE == "local" && !collect_only && !file.exists(STILT_EXE))
  stop("STILT executable not found: ", STILT_EXE, " (set METHANE_STILT_EXE, or METHANE_STILT_ENGINE=external)")
cat(sprintf("%d receptors, %d particles, %g h, engine %s, %d cores\n", nrow(R), NUMPAR, N_HOURS,
            if (collect_only) "collect" else STILT_ENGINE, N_CORES))
t0 <- Sys.time()
res <- unlist(mclapply(seq_len(nrow(R)), one, mc.cores = N_CORES, mc.preschedule = FALSE))
st <- data.frame(id = R$id, status = ifelse(grepl("foot\\.rds$", res), "ok", sub(" .*", "", res)), detail = ifelse(grepl("foot\\.rds$", res), "", res))
write.csv(st, file.path(TRANSPORT_DIR, "footprint_status.csv"), row.names = FALSE)
cat(sprintf("done in %.1f min: %s\n", as.numeric(difftime(Sys.time(), t0, units = "mins")),
            paste(names(table(st$status)), table(st$status), collapse = ", ")))
if (any(st$status == "MISSING") && STILT_ENGINE == "external")
  cat("run directories are ready in", FOOT_DIR, "- run them with run_stilt_linux.sh, then:\n  Rscript scripts/03_footprints.R --collect\n")
bad <- st[!st$status %in% c("ok", "MISSING"), ]
if (nrow(bad)) { cat("problems:\n"); print(head(bad, 10), row.names = FALSE) }
