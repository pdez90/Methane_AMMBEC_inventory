# 02_met.R — HRRR 3-km meteorology for every receptor ------------------------------------------
# Lists the NOAA ARL HRRR 6-hour files that the footprints need (N_HOURS back from each segment,
# plus the file after, for interpolation), and makes sure each is present in <MET_DIR>/sub as a
# DOMAIN subgrid: files Paper 1 script 64 already cut are reused; the rest are downloaded
# (~3.4 GB each, resumable) and cut with xtrct_grid, and the full file is deleted.
# Needs internet access to www.ready.noaa.gov or ftp.arl.noaa.gov — run it on the Mac.
#
#   Rscript scripts/02_met.R            # download + subset what is missing
#   Rscript scripts/02_met.R --list     # only list what is needed / missing
#   Out: <MET_DIR>/sub/*_hrrr, <INV_OUT>/met_files.csv
# ----------------------------------------------------------------------------------------------
proj <- if (file.exists("config.R")) "." else ".."; source(file.path(proj, "config.R"))
source(file.path(proj, "R", "stilt.R"))

S <- read.csv(file.path(INV_OUT, "obs_segments.csv"), stringsAsFactors = FALSE)
tt <- as.POSIXct(S$time_utc, tz = "UTC")
need <- sort(unique(unlist(lapply(tt, hrrr_files_for))))
# Subgrids already cut on the Elements drive are linked (read-only), not copied
dir.create(file.path(MET_DIR, "sub"), FALSE, TRUE)
for (f in need) { src <- file.path(MET_READ_DIR, f); dst <- file.path(MET_DIR, "sub", f)
  if (!file.exists(dst) && file.exists(src) && file.size(src) > 1e6) file.symlink(src, dst) }
have <- file.exists(file.path(MET_DIR, "sub", need)) & file.size(file.path(MET_DIR, "sub", need)) > 1e6
have[is.na(have)] <- FALSE
write.csv(data.frame(file = need, present = have), file.path(INV_OUT, "met_files.csv"), row.names = FALSE)
cat(sprintf("%d HRRR files needed, %d already in %s, %d to fetch (~%.0f GB download, ~%.0f MB kept)\n",
            length(need), sum(have), file.path(MET_DIR, "sub"), sum(!have), 3.4 * sum(!have), 65 * sum(!have)))
if ("--list" %in% commandArgs(TRUE)) { if (any(!have)) cat("missing:", need[!have], sep = "\n  "); quit(save = "no") }
# NOAA's server gives each connection ~1 MB/s, so several files download at once
# (METHANE_HRRR_PARALLEL, default 4; ~3.4 GB of disk each while in flight). Interrupting is safe:
# partial downloads resume with curl -C and finished subgrids are skipped.
np <- as.integer(Sys.getenv("METHANE_HRRR_PARALLEL", "4"))
todo <- need[!have]
cat(sprintf("downloading %d at a time; a line is printed as each file finishes\n", np))
t0 <- Sys.time()
res <- parallel::mclapply(seq_along(todo), function(i) {
  f <- todo[i]
  r <- tryCatch({ ensure_hrrr(f, quiet = np > 1); "ok" }, error = function(e) conditionMessage(e))
  cat(sprintf("[%s] %-22s %s (%.0f min since start)\n", format(Sys.time(), "%H:%M"), f, r,
              as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  r
}, mc.cores = np, mc.preschedule = FALSE)
bad <- todo[unlist(res) != "ok"]
if (length(bad)) stop(length(bad), " file(s) failed (run 02_met.R again to retry): ", paste(bad, collapse = " "))
cat("all", length(need), "files ready in", file.path(MET_DIR, "sub"), "\n")
