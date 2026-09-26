# stilt.R — HRRR meteorology and STILT (HYSPLIT STILT-emulation) runs, base R ------------------
# The SETUP.CFG / CONTROL writers and the near-field dilution follow STILT v2 (Fasoli et al. 2018,
# uataq/stilt r/src); they are re-implemented here so that nothing but hycs_std is needed from
# STILT. Validated in Paper 1 script 64 against the uataq Linux build.
# ----------------------------------------------------------------------------------------------

STILT_VARS <- c("time", "indx", "long", "lati", "zagl", "foot", "mlht", "dens", "samt", "sigw", "tlgr")

## ---- HRRR ARL files (NOAA ARL archive, 6-hour files YYYYMMDD_HH-HH_hrrr) ------------------------
hrrr_name <- function(t) {
  h0 <- (as.integer(format(t, "%H", tz = "UTC")) %/% 6) * 6
  sprintf("%s_%02d-%02d_hrrr", format(t, "%Y%m%d", tz = "UTC"), h0, h0 + 5)
}
hrrr_files_for <- function(t, n_hours = N_HOURS) {
  tt <- c(seq(t + n_hours * 3600, t, by = 3600), t, t + 3600)   # +1 h: interpolation at file ends
  unique(vapply(tt, hrrr_name, ""))
}
HRRR_URLS <- c(Sys.getenv("METHANE_HRRR_URL", ""),
               "https://www.ready.noaa.gov/data/archives/hrrr/",
               "ftp://ftp.arl.noaa.gov/archives/hrrr/")
HRRR_URLS <- HRRR_URLS[nzchar(HRRR_URLS)]

.curl_head <- function(u) {
  out <- suppressWarnings(system2("curl", c("-sI", "-L", "--max-time", "30", shQuote(u)), stdout = TRUE, stderr = FALSE))
  len <- suppressWarnings(as.numeric(gsub("\\D", "", tail(grep("^content-length", out, ignore.case = TRUE, value = TRUE), 1))))
  list(ok = is.null(attr(out, "status")) && length(len) == 1 && is.finite(len) && len > 1e8, bytes = if (length(len)) len else NA)
}
arl_nz <- function(path) {          # levels in an ARL file header (after STILT read_met_header.r)
  h <- substring(readChar(path, 166), 51)
  as.integer(as.numeric(unlist(strsplit(gsub("(.{3})", "\\1_", substring(h, 94, 102)), "_")))[3])
}
#' Make sure the HRRR subgrid for `f` exists in <MET_DIR>/sub: download (resumable) and cut it with
#' xtrct_grid to DOMAIN, then delete the ~3.4 GB original unless keep_raw.
ensure_hrrr <- function(f, xtrct = XTRCT_EXE, quiet = FALSE,
                        keep_raw = Sys.getenv("METHANE_HRRR_KEEP_RAW", "0") == "1") {
  sub <- file.path(MET_DIR, "sub", f); raw <- file.path(MET_DIR, "raw", f)
  dir.create(dirname(sub), FALSE, TRUE); dir.create(dirname(raw), FALSE, TRUE)
  if (file.exists(sub) && file.size(sub) > 1e6) { if (!keep_raw && file.exists(raw)) unlink(raw); return(invisible(sub)) }
  got <- FALSE
  # A dropped connection (Wi-Fi blip, NOAA hiccup, TLS error) used to fail the file at once, so a
  # 1-minute outage burned through the whole queue. Now: up to 8 tries, 2 min apart (~15 min).
  for (attempt in 1:8) {
  for (b in HRRR_URLS) {
    h <- .curl_head(paste0(b, f)); if (!h$ok) next
    if (!(file.exists(raw) && file.size(raw) == h$bytes)) {
      message(sprintf("  downloading %s (%.1f GB)", f, h$bytes / 1e9))
      system2("curl", c(if (quiet) "-sS", "-f", "-L", "-C", "-", "--retry", "10", "--retry-all-errors", "--retry-delay", "30", "-o", shQuote(raw), shQuote(paste0(b, f))))
    }
    if (file.exists(raw) && file.size(raw) == h$bytes) { got <- TRUE; break }
  }
  if (got) break
  if (attempt < 8) { message(sprintf("  %s: server not reachable (try %d of 8), waiting 2 min", f, attempt)); Sys.sleep(120) }
  }
  if (!got) stop("could not download ", f, " from ", paste(HRRR_URLS, collapse = ", "))
  wd <- tempfile("xtrct_"); dir.create(wd)
  writeLines(c(paste0(dirname(raw), "/"), basename(raw), paste(DOMAIN["ymn"], DOMAIN["xmn"]),
               paste(DOMAIN["ymx"], DOMAIN["xmx"]), arl_nz(raw)), file.path(wd, "input"))
  st <- system(paste("cd", shQuote(wd), "&&", shQuote(xtrct), "< input > xtrct.log 2>&1"))
  if (st %in% c(9, 137)) stop("xtrct_grid was killed by macOS; run: xattr -dr com.apple.quarantine ", dirname(xtrct),
                              " && codesign --force -s - ", xtrct)
  if (!file.exists(file.path(wd, "extract.bin"))) stop("xtrct_grid failed on ", f, "; see ", file.path(wd, "xtrct.log"))
  file.copy(file.path(wd, "extract.bin"), sub, overwrite = TRUE); unlink(wd, recursive = TRUE)
  if (!keep_raw) unlink(raw)
  invisible(sub)
}

## ---- run files ----------------------------------------------------------------------------------
stilt_write_setup <- function(file, numpar = NUMPAR, zicontrol = ZISCALE != 1) {   # STILT v2 run_stilt.r defaults
  kv <- list(CAPEMIN = -1, CMASS = 0, CONAGE = 48, CPACK = 1, DELT = 1, DXF = 1, DYF = 1, DZF = 0.01,
    EFILE = "''", FRHMAX = 3, FRHS = 1, FRME = 0.1, FRMR = 0, FRTS = 0.1, FRVS = 0.01, HSCALE = 10800,
    ICHEM = 8, IDSP = 2, INITD = 0, IVMAX = length(STILT_VARS), K10M = 1, KAGL = 1, KBLS = 1, KBLT = 5, KDEF = 0,
    KHINP = 0, KHMAX = 9999, KMIX0 = 150, KMIXD = 3, KMSL = 0, KPUFF = 0, KRAND = 4, KRND = 6, KSPL = 1,
    KWET = 1, KZMIX = 0, MAXDIM = 1, MAXPAR = numpar, MGMIN = 10, MHRS = 9999, NBPTYP = 1, NCYCL = 0,
    NDUMP = 0, NINIT = 1, NSTR = 0, NTURB = 0, NUMPAR = numpar, NVER = 0, OUTDT = 0, P10F = 1,
    PINBC = "''", PINPF = "''", POUTF = "''", QCYCLE = 0, RHB = 80, RHT = 60, SPLITF = 1, TKERD = 0.18,
    TKERN = 0.18, TLFRAC = 0.1, TOUT = 0, TRATIO = 0.75, TVMIX = 1,
    VARSIWANT = paste0("'", paste(toupper(STILT_VARS), collapse = "', '"), "'"),
    VEGHT = 0.5, VSCALE = 200, VSCALEU = 200, VSCALES = -1, WBBH = 0, WBWF = 0, WBWR = 0, WVERT = ".FALSE.",
    WINDERRTF = 0, ZICONTROLTF = as.integer(zicontrol))
  writeLines(c("$SETUP", paste0(names(kv), "=", vapply(kv, function(v) format(v, scientific = FALSE), ""), ","), "$END"), file)
}
stilt_write_control <- function(r, met, file, n_hours = N_HOURS) {
  writeLines(c(format(r$run_time, "%y %m %d %H %M", tz = "UTC"), "1", paste(r$lati, r$long, r$zagl),
    n_hours, "0", "25000.0", length(met), as.vector(rbind(paste0(dirname(met), "/"), basename(met))),
    "1", "test", "1", "0.01", "00 00 00 00 00", "1", "0.0 0.0", "0.5 0.5", "30.0 30.0", "./", "cdump",
    "1", "100", "00 00 00 00 00", "00 00 00 00 00", "00 2 00", "1", "0.0 0.0 0.0",
    "0.0 0.0 0.0 0.0 0.0", "0.0 0.0 0.0", "0.0", "0.0"), file)
}
#' Write the run directory for one receptor (CONTROL, SETUP.CFG, ASCDATA.CFG pointing at ./ and
#' copies of the land-use grids). The CONTROL file lists met files by absolute path on this machine;
#' run_stilt_linux.sh rewrites that prefix when the runs happen elsewhere.
stilt_prepare <- function(r, wd) {
  dir.create(wd, FALSE, TRUE)
  stilt_write_setup(file.path(wd, "SETUP.CFG"))
  stilt_write_control(r, file.path(MET_DIR, "sub", hrrr_files_for(r$run_time)), file.path(wd, "CONTROL"))
  # Mixed-layer sensitivity (STILT ziscale): HYSPLIT multiplies its diagnosed mixed-layer height by
  # the factor for each hour back. ZICONTROL = number of hours, then one factor per line.
  if (ZISCALE != 1) writeLines(c(abs(N_HOURS), rep(format(ZISCALE), abs(N_HOURS))), file.path(wd, "ZICONTROL"))
  asc <- c("LANDUSE.ASC", "ROUGLEN.ASC", "TERRAIN.ASC")
  for (f in c("ASCDATA.CFG", asc)) if (file.exists(file.path(STILT_BDY, f))) file.copy(file.path(STILT_BDY, f), wd, overwrite = TRUE)
  if (file.exists(file.path(wd, "ASCDATA.CFG"))) {
    a <- readLines(file.path(wd, "ASCDATA.CFG")); a[length(a)] <- "'./'           directory location of data files"
    writeLines(a, file.path(wd, "ASCDATA.CFG"))
  }
  invisible(wd)
}

## ---- particle output ----------------------------------------------------------------------------
read_particles <- function(pf) {
  con <- if (grepl("\\.gz$", pf)) gzfile(pf) else file(pf); ln <- readLines(con); close(con)
  nf <- lengths(strsplit(trimws(ln), "[[:space:]]+"))
  ln <- ln[nf == length(STILT_VARS)]
  if (!length(ln)) return(NULL)
  x <- read.table(text = ln, header = FALSE, colClasses = "numeric"); names(x) <- STILT_VARS; x
}
plume_dilution <- function(p, r_zagl, veght = 0.5) {    # STILT calc_plume_dilution (hnf_plume = TRUE)
  p <- p[order(p$indx, -p$time), ]
  sig <- p$samt * sqrt(2) * p$sigw * sqrt(p$tlgr * abs(p$time * 60) + p$tlgr^2 * exp(-abs(p$time * 60) / p$tlgr) - 1)
  sig[!is.finite(sig)] <- 0
  plume <- r_zagl + ave(sig, p$indx, FUN = cumsum)
  near <- plume < veght * p$mlht
  p$foot[near] <- (0.02897 / (plume * p$dens) * p$samt * 60)[near]
  p
}
#' Turn a finished run directory into the compact footprint file <wd>/foot.rds:
#' only rows with foot > 0 are kept (time [min, negative], indx, long, lati, foot [ppm per umol m-2 s-1])
#' plus attributes np (number of particles) and mlht_receptor (median HRRR mixed-layer height, m AGL,
#' over the first 15 min). Returns the path, or an error string.
stilt_collect <- function(r, wd) {
  out <- file.path(wd, "foot.rds"); if (file.exists(out)) return(out)
  pf <- c(file.path(wd, "PARTICLE_STILT.DAT.gz"), file.path(wd, "PARTICLE_STILT.DAT"))
  pf <- pf[file.exists(pf)][1]
  if (is.na(pf)) return(paste("MISSING", basename(wd)))
  p <- read_particles(pf); if (is.null(p) || !nrow(p)) return(paste("EMPTY", basename(wd)))
  np <- max(p$indx); p <- plume_dilution(p, r$zagl)
  # mixed-layer height seen by the particles in their first 15 min (receptor-local HRRR mlht, m AGL),
  # for the check against the Doppler-lidar boundary-layer height (05b)
  mlh0 <- p$mlht[p$time >= -15]
  p <- p[p$foot > 0, c("time", "indx", "long", "lati", "foot")]   # indx kept for particle-number convergence (05b)
  attr(p, "np") <- np
  attr(p, "mlht_receptor") <- if (length(mlh0)) stats::median(mlh0) else NA_real_
  saveRDS(p, out)
  unlink(file.path(wd, c("PARTICLE_STILT.DAT", "PARTICLE_STILT.DAT.gz", "PARTICLE.DAT", "cdump",
                         "LANDUSE.ASC", "ROUGLEN.ASC", "TERRAIN.ASC")))
  out
}
stilt_run_local <- function(r, wd, exe = STILT_EXE) {
  if (file.exists(file.path(wd, "foot.rds"))) return(file.path(wd, "foot.rds"))
  if (!any(file.exists(file.path(wd, c("PARTICLE_STILT.DAT.gz"))))) {
    stilt_prepare(r, wd)
    system(paste("cd", shQuote(wd), "&&", shQuote(exe), "> stilt.log 2>&1 < /dev/null"), timeout = 3600)
  }
  stilt_collect(r, wd)
}
