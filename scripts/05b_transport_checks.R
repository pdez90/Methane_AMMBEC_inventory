# 05b_transport_checks.R — the three STILT choices a reviewer will ask about ----------------------
# Everything else in SETUP.CFG is the STILT v2 default (Lin et al. 2003; Fasoli et al. 2018) and is
# cited, not tuned. The settings chosen here are checked:
#
#  (a) PARTICLE NUMBER. Particles in each foot.rds are subsampled (n = 50...np, 40 draws) and the
#      prior-predicted CH4 enhancement recomputed. The relative particle noise of a full run is
#      extrapolated from the subsamples (finite-population corrected) and reported for np = 300 and
#      the current NUMPAR, next to the inversion's model-data mismatch.
#  (b) BACK-TRAJECTORY LENGTH (N_HOURS). Share of each receptor's prior signal picked up in the last
#      hour, and the share of particle positions still inside the flux domain. A small last-hour
#      share means a longer run would change H little.
#  (c) MIXED-LAYER HEIGHT. HRRR mlht seen by the particles at the receptor versus the Doppler-lidar
#      boundary-layer height (Paper 1 script 08, blh_timeseries.csv) for urban segments within 30 min.
#      If HRRR is consistently off, run 03-07 with METHANE_ZISCALE = lidar/HRRR; in any case run the
#      +/-20% test (METHANE_ZISCALE=0.8 and 1.2) and report the change in the posterior.
#
#   Rscript scripts/05b_transport_checks.R [prior=gra2pes_v2.0beta] [n_receptors=80]
#   Out: <TRANSPORT_DIR>/transport_checks.csv, transport_checks.txt, figures/transport_checks.png
# ----------------------------------------------------------------------------------------------
proj <- if (file.exists("config.R")) "." else ".."; source(file.path(proj, "config.R"))
source(file.path(proj, "R", "grids.R"))
a <- commandArgs(TRUE)
PRIOR <- if (length(a) >= 1) a[1] else "gra2pes_v2.0beta"
NREC  <- if (length(a) >= 2) as.integer(a[2]) else 80L
set.seed(42)

S <- read.csv(file.path(INV_OUT, "obs_segments.csv"), stringsAsFactors = FALSE)
S$foot <- file.path(FOOT_DIR, S$id, "foot.rds"); S <- S[file.exists(S$foot), ]
S$run_time <- as.POSIXct(round(as.numeric(as.POSIXct(S$time_utc, tz = "UTC")) / 60) * 60, origin = "1970-01-01", tz = "UTC")
g <- readRDS(file.path(PRIOR_DIR, paste0("prior_", PRIOR, ".rds")))
Ftot <- Reduce(`+`, g$comps)
cat(sprintf("transport %s: %d footprints, prior %s\n", TRANSPORT, nrow(S), PRIOR))

per_particle <- function(i) {        # prior-predicted CH4 (ppb) carried by each particle
  p <- readRDS(S$foot[i]); np <- attr(p, "np"); mlr <- attr(p, "mlht_receptor")
  if (is.null(p$indx)) stop("foot.rds has no particle index; rerun 03 with the current R/stilt.R")
  ci <- grid_index(g, p$long, p$lati); inside <- mean(ci$ok)
  late_pos <- p$time < (N_HOURS + 1) * 60
  inside_late <- if (any(late_pos)) mean(ci$ok[late_pos]) else NA
  q <- p[ci$ok, ]; ix <- ci$ix[ci$ok]; iy <- ci$iy[ci$ok]
  h <- as.integer(format(S$run_time[i] + q$time * 60, "%H", tz = "UTC"))
  cc <- q$foot * Ftot[ix + (iy - 1L) * length(g$x) + h * length(g$x) * length(g$y)] * 1000
  v <- numeric(np); t <- tapply(cc, q$indx, sum); v[as.integer(names(t))] <- t
  late <- sum(cc[q$time < (N_HOURS + 1) * 60]) / max(sum(cc), 1e-12)
  list(v = v, np = np, late = late, inside_late = inside_late, mlht = if (is.null(mlr)) NA else mlr)
}
idx <- sort(sample(seq_len(nrow(S)), min(NREC, nrow(S))))
pp <- parallel::mclapply(seq_len(nrow(S)), function(i) tryCatch(per_particle(i), error = function(e) conditionMessage(e)), mc.cores = N_CORES)
bad <- vapply(pp, is.character, NA); if (any(bad)) stop(pp[[which(bad)[1]]])

## (a) particle convergence
NS <- c(50, 100, 200, 300, 500)
conv <- do.call(rbind, lapply(idx, function(i) {
  v <- pp[[i]]$v; N <- length(v); H <- mean(v)
  if (H <= 0.5) return(NULL)                     # skip receptors with essentially no prior signal
  do.call(rbind, lapply(NS[NS < N], function(n) {
    hs <- replicate(40, mean(sample(v, n)))
    rs <- sd(hs) / H
    sig_p <- rs * sqrt(n) / sqrt(1 - n / N)      # per-particle relative sd
    data.frame(id = S$id[i], region = S$region[i], H = H, n = n, rel_sd_n = rs,
               rel_sd_300 = sig_p / sqrt(300), rel_sd_np = sig_p / sqrt(N), np = N)
  }))
}))
cv <- aggregate(cbind(rel_sd_300, rel_sd_np) ~ id + region + H + np, data = conv, FUN = median)

## (b) back-trajectory length
tl <- data.frame(id = S$id, region = S$region,
                 late_share = vapply(pp, `[[`, 0, "late"), inside_late = vapply(pp, `[[`, 0, "inside_late"),
                 mlht = vapply(pp, `[[`, 0, "mlht"))

## (c) HRRR mixed-layer height vs Doppler lidar (urban segments, daytime)
bl_f <- file.path(DATA_BASE, "outputs", "blh_timeseries.csv"); mlh <- NULL
if (file.exists(bl_f)) {
  bl <- read.csv(bl_f); bl$t <- as.POSIXct(bl$time, tz = "UTC")
  bl <- bl[is.finite(bl$blh_m) & bl$blh_m > 150 & bl$blh_m < 4500, ]
  U <- which(S$region == "urban" & is.finite(tl$mlht))
  mlh <- do.call(rbind, lapply(U, function(i) {
    k <- abs(as.numeric(difftime(bl$t, S$run_time[i], units = "mins"))) <= 30
    if (sum(k) < 1) return(NULL)
    data.frame(id = S$id[i], hrrr_mlht = tl$mlht[i], lidar_blh = median(bl$blh_m[k]), agl = S$agl_m[i])
  }))
} else message("note: ", bl_f, " not found; MLH check skipped")

out <- merge(tl, cv[, c("id", "H", "rel_sd_300", "rel_sd_np")], by = "id", all.x = TRUE)
write.csv(out, file.path(TRANSPORT_DIR, "transport_checks.csv"), row.names = FALSE)
if (!is.null(mlh)) write.csv(mlh, file.path(TRANSPORT_DIR, "transport_mlh_vs_lidar.csv"), row.names = FALSE)

q <- function(v, p = c(.5, .1, .9)) round(quantile(v, p, na.rm = TRUE), 3)
sink(file.path(TRANSPORT_DIR, "transport_checks.txt"), split = TRUE)
cat(sprintf("\n=== transport checks, %s (NUMPAR %d, N_HOURS %d, ZISCALE %g) ===\n", TRANSPORT, NUMPAR, N_HOURS, ZISCALE))
cat(sprintf("(a) particle noise in the prior-predicted enhancement, %d receptors with H > 0.5 ppb:\n", nrow(cv)))
cat(sprintf("    relative sd at np = 300: median %.1f%% (10-90%%: %.1f-%.1f%%)\n", 100*median(cv$rel_sd_300), 100*q(cv$rel_sd_300)[2], 100*q(cv$rel_sd_300)[3]))
cat(sprintf("    relative sd at np = %d: median %.1f%% (10-90%%: %.1f-%.1f%%)\n", cv$np[1], 100*median(cv$rel_sd_np), 100*q(cv$rel_sd_np)[2], 100*q(cv$rel_sd_np)[3]))
cat("    (compare with the inversion's multiplicative model-data mismatch phi ~ 0.6, i.e. 60%)\n")
cat(sprintf("(b) share of the prior signal from the last hour (%d to %d h): median %.1f%%, 90th pct %.1f%%\n",
            N_HOURS + 1, N_HOURS, 100*median(tl$late_share, na.rm = TRUE), 100*quantile(tl$late_share, .9, na.rm = TRUE)))
cat(sprintf("    particle positions still inside the flux domain in the last hour: median %.0f%%\n", 100*median(tl$inside_late, na.rm = TRUE)))
if (!is.null(mlh) && nrow(mlh)) {
  r <- mlh$hrrr_mlht / mlh$lidar_blh
  cat(sprintf("(c) HRRR mlht / lidar BLH, %d urban segments within 30 min: median %.2f (IQR %.2f-%.2f); lidar median %.0f m, HRRR %.0f m\n",
              nrow(mlh), median(r), quantile(r, .25), quantile(r, .75), median(mlh$lidar_blh), median(mlh$hrrr_mlht)))
  cat(if (abs(median(r) - 1) > 0.2) sprintf("    HRRR/lidar differ by more than 20%% (lidar-implied ZISCALE %.2f). This ratio is driven by a few afternoons with\n    HRRR mlht > 4 km; use it as a bracket, not the central run. Set the central ZISCALE from the aircraft profiles (05c).\n", 1 / median(r)) else
      "    within 20%: keep ZISCALE = 1 and report the +/-20% sensitivity runs\n")
} else cat("(c) no HRRR mlht in these footprints (np300 runs predate it) or no lidar match; rerun 03 first\n")
sink()

png(file.path(TRANSPORT_DIR, "figures", "transport_checks.png"), 1800, 600, res = 150); par(mfrow = c(1, 3), mar = c(4.5, 4.5, 2.5, 1))
m <- aggregate(rel_sd_n ~ n, data = conv, FUN = median)
plot(m$n, 100 * m$rel_sd_n, log = "xy", type = "b", pch = 19, xlab = "particles", ylab = "relative sd of H (%)",
     main = "(a) particle convergence", cex.main = 0.9); lines(m$n, 100 * m$rel_sd_n[1] * sqrt(m$n[1] / m$n), lty = 3)
abline(v = c(300, NUMPAR), col = c("gray50", "#d73027"), lty = 2)
hist(100 * tl$late_share, breaks = 30, col = "#4575b4", border = "white", xlab = "% of prior signal from the last hour",
     main = sprintf("(b) back-trajectory length (%d h)", abs(N_HOURS)), cex.main = 0.9)
if (!is.null(mlh) && nrow(mlh)) { plot(mlh$lidar_blh, mlh$hrrr_mlht, pch = 19, col = "#00000066", xlab = "lidar BLH (m)", ylab = "HRRR mlht at receptor (m)",
  main = "(c) mixed-layer height", cex.main = 0.9); abline(0, 1); abline(0, 1.2, lty = 3); abline(0, 0.8, lty = 3)
} else { plot.new(); title("(c) no MLH data yet", cex.main = 0.9) }
invisible(dev.off())
cat("wrote transport_checks.csv/.txt and figures/transport_checks.png to", TRANSPORT_DIR, "\n")
