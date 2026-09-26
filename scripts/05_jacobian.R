# 05_jacobian.R — footprints x prior fields: the forward model H ---------------------------------
# For every segment with a footprint and every prior from 04, the predicted enhancement from each
# CH4 component k and from total CO:
#
#   H[i, k] = (1/Np) * sum_particles sum_t foot_p,t * F_k(cell(x_p,t), hour(t))      (ppb)
#
# foot is STILT's sensitivity (ppm per umol m-2 s-1, with the near-field dilution applied in 03),
# F the prior flux (umol m-2 s-1) in the particle's cell at the particle's UTC hour. Particles are
# placed on the GRA2PES Lambert grid with the spherical LCC in R/grids.R (checked against the
# file coordinates in 04). CO uses the weekday / Saturday / Sunday field of the particle's local
# date where the prior has one. Particles that leave DOMAIN contribute nothing (that air is
# background, which the inversion estimates).
#
# With alpha_k = 1 for every k, rowSums(H) is the prior's forward simulation, so this script also
# writes the prior-versus-observed comparison per segment — a first test of v1.1 against v2.0beta.
#
#   Rscript scripts/05_jacobian.R [prior names ...]      (default: every prior_*.rds in PRIOR_DIR)
#   Out: <RUN_DIR>/jacobian_<prior>.rds, jacobian_<prior>.csv, forward_summary.csv
# ----------------------------------------------------------------------------------------------
proj <- if (file.exists("config.R")) "." else ".."; source(file.path(proj, "config.R"))
source(file.path(proj, "R", "grids.R"))
suppressPackageStartupMessages(library(parallel))

priors <- commandArgs(TRUE)
if (!length(priors)) priors <- sub("^prior_(.*)\\.rds$", "\\1", list.files(PRIOR_DIR, pattern = "^prior_.*\\.rds$"))
if (!length(priors)) stop("no prior_*.rds in ", PRIOR_DIR, " (run scripts/04_priors.R)")
S <- read.csv(file.path(INV_OUT, "obs_segments.csv"), stringsAsFactors = FALSE)
S$foot <- file.path(FOOT_DIR, S$id, "foot.rds"); S <- S[file.exists(S$foot), ]
if (!nrow(S)) stop("no foot.rds under ", FOOT_DIR, " (run scripts/03_footprints.R)")
S$run_time <- as.POSIXct(round(as.numeric(as.POSIXct(S$time_utc, tz = "UTC")) / 60) * 60, origin = "1970-01-01", tz = "UTC")
cat(sprintf("%d segments with footprints\n", nrow(S)))

daytype_of <- function(t) {                     # local (UTC-6) calendar day
  wd <- as.POSIXlt(t - 6 * 3600, tz = "UTC")$wday
  ifelse(wd == 6, "satdy", ifelse(wd == 0, "sundy", "weekdy"))
}
jac_one <- function(i, g) {
  p <- readRDS(S$foot[i]); np <- attr(p, "np"); mlr <- attr(p, "mlht_receptor"); if (is.null(mlr)) mlr <- NA_real_
  ci <- grid_index(g, p$long, p$lati); p <- p[ci$ok, ]; ix <- ci$ix[ci$ok]; iy <- ci$iy[ci$ok]
  nx <- length(g$x); ny <- length(g$y)
  tabs <- S$run_time[i] + p$time * 60
  h <- if (length(g$hours) == 1 && is.na(g$hours)) rep(0L, nrow(p)) else as.integer(format(tabs, "%H", tz = "UTC"))
  lin <- ix + (iy - 1L) * nx + h * nx * ny
  Hk <- vapply(g$comps, function(a) sum(p$foot * a[lin]) / np * 1000, 0)
  hco <- NA_real_
  if (!is.null(g[["co"]])) {
    dt <- daytype_of(tabs); dt[!dt %in% names(g[["co"]])] <- "weekdy"
    hco <- 0; for (d in unique(dt)) { k <- dt == d; hco <- hco + sum(p$foot[k] * g[["co"]][[d]][lin[k]]) }
    hco <- hco / np * 1000
  }
  # share of the prior signal picked up in the last hour of the back-trajectory: small means the
  # footprint has converged in time and N_HOURS is long enough (05b reports the distribution)
  late <- p$time < (N_HOURS + 1) * 60
  contrib <- p$foot * Reduce(`+`, lapply(g$comps, function(a) a[lin])); tot <- sum(contrib)
  late_share <- if (tot > 0) sum(contrib[late]) / tot else NA
  c(Hk, co = hco, foot_sum = sum(p$foot) / np, frac_in_domain = if (length(ci$ok)) mean(ci$ok) else NA,
    late_share = late_share, mlht_receptor = mlr)
}
fs <- list()
for (pr in priors) {
  g <- readRDS(file.path(PRIOR_DIR, paste0("prior_", pr, ".rds")))
  t0 <- Sys.time()
  M <- do.call(rbind, mclapply(seq_len(nrow(S)), jac_one, g = g, mc.cores = N_CORES))
  K <- names(g$comps)
  J <- list(prior = pr, ids = S$id, components = K, H = M[, K, drop = FALSE], Hco = M[, "co"],
            foot_sum = M[, "foot_sum"], frac_in_domain = M[, "frac_in_domain"], late_share = M[, "late_share"],
            mlht_receptor = M[, "mlht_receptor"], transport = TRANSPORT, box_t_hr = g$box_t_hr, region_t_hr = g$region_t_hr,
            has_co = !is.null(g[["co"]]))
  saveRDS(J, file.path(RUN_DIR, paste0("jacobian_", pr, ".rds")))
  out <- data.frame(id = S$id, flight = S$flight, region = S$region, agl_m = round(S$agl_m),
                    ch4_enh_p1 = S$ch4_enh_p1, co_enh_p1 = S$co_enh_p1, round(M, 4), check.names = FALSE)
  out$ch4_prior_total <- round(rowSums(M[, K, drop = FALSE]), 4)
  write.csv(out, file.path(RUN_DIR, paste0("jacobian_", pr, ".csv")), row.names = FALSE)
  low <- S$agl_m < 1000
  sp <- function(a, b, k = TRUE) { ok <- k & is.finite(a) & is.finite(b); if (sum(ok) < 5) NA else round(suppressWarnings(cor(a[ok], b[ok], method = "spearman")), 2) }
  fs[[pr]] <- data.frame(prior = pr, n = nrow(S),
    median_pred_ch4_ppb = round(median(out$ch4_prior_total), 2), median_obs_enh_ppb = round(median(S$ch4_enh_p1), 2),
    rho_ch4 = sp(out$ch4_prior_total, S$ch4_enh_p1), rho_ch4_low = sp(out$ch4_prior_total, S$ch4_enh_p1, low),
    median_pred_co_ppb = round(median(M[, "co"]), 2), rho_co = sp(M[, "co"], S$co_enh_p1),
    rho_co_low = sp(M[, "co"], S$co_enh_p1, low), minutes = round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1))
  cat(sprintf("%-18s %d segments in %.1f min; median predicted CH4 %.1f ppb (observed rolling-baseline enhancement %.1f)\n",
              pr, nrow(S), fs[[pr]]$minutes, fs[[pr]]$median_pred_ch4_ppb, fs[[pr]]$median_obs_enh_ppb))
}
FS <- do.call(rbind, fs); rownames(FS) <- NULL
write.csv(FS, file.path(RUN_DIR, "forward_summary.csv"), row.names = FALSE)
print(FS, row.names = FALSE)
