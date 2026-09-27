# 05c_blh_aircraft.R — mixed-layer height from the aircraft's own vertical profiles ------------------
# 05b found HRRR's mixed-layer height (as the STILT particles see it) about 2.7x the Doppler-lidar
# BLH, and the ZISCALE runs fit the observations better as the mixed layer is made shallower. Before
# adopting a lidar-based scaling for the main run, this checks both against a third, independent
# estimate: the Twin Otter's climbs and descents.
#
# Profiles: stretches where the aircraft climbs or descends steadily (|dz/dt| > VS_MIN m/s from GPS
# altitude, smoothed over 15 s) through at least SPAN_MIN m, starting or ending within BOTTOM_MAX m
# of the ground. Each is binned every BIN_M m (AGL) and three estimates are made:
#   zi_theta  base of the lowest inversion (Heffter 1980): a layer with dtheta/dz >= GRAD_MIN K/m across
#             which theta rises by >= DTHETA K (the primary estimate). Climb temperatures can lag, so
#             descents are the more reliable; both are kept and flagged.
#   zi_h2o    height of the sharpest drop in water vapour (most negative dq/dz)
#   zi_ch4    height of the sharpest drop in CH4
# If theta never jumps below the profile top, the mixed layer is at least as deep as the top: the
# estimate is censored (zi_theta = NA, above_top = TRUE, zi_lower_bound = top). These profiles still
# count: a lidar BLH lower than a censored aircraft bound is evidence the lidar is too low.
# Pressure comes from the pressure altitude (standard atmosphere); temperature from Temp_True.
#
# Each profile is paired with (a) the lidar BLH within +/-30 min (Paper 1 outputs/blh_timeseries.csv)
# and (b) HRRR mlht at the receptors of segments on the same flight within +/-60 min and 30 km
# (jacobian of the unscaled transport run, ZISCALE = 1).
#
#   Rscript scripts/05c_blh_aircraft.R
#   Out: <TRANSPORT_DIR of the ZISCALE=1 run>/blh_aircraft_profiles.csv, blh_aircraft_summary.txt,
#        figures/blh_aircraft.png, figures/blh_aircraft_profiles.png
# ----------------------------------------------------------------------------------------------
proj <- if (file.exists("config.R")) "." else ".."; source(file.path(proj, "config.R"))
VS_MIN <- 1.0; SPAN_MIN <- 400; BOTTOM_MAX <- 500; BIN_M <- 50; DTHETA <- 2.0; GRAD_MIN <- 0.005; MIN_DUR_S <- 60
SURF_MARGIN <- 100   # m: a stable layer whose base is this close to the profile bottom is the surface layer
                     # (takeoff / landing through a superadiabatic-then-stable surface layer), not an inversion base
REF_TRANSPORT <- Sys.getenv("METHANE_BLH_REF_TRANSPORT", sprintf("np%d_h%d", NUMPAR, abs(N_HOURS)))
REF_DIR <- file.path(INV_OUT, "runs", REF_TRANSPORT); dir.create(file.path(REF_DIR, "figures"), FALSE, TRUE)
hav <- function(la1, lo1, la2, lo2) { p <- pi / 180
  a <- sin((la2 - la1) * p / 2)^2 + cos(la1 * p) * cos(la2 * p) * sin((lo2 - lo1) * p / 2)^2
  2 * 6371.0088 * asin(pmin(1, sqrt(a))) }
p_from_palt <- function(z) 1013.25 * (1 - 2.25577e-5 * z)^5.25588          # hPa, standard atmosphere

## ---- reference data: lidar BLH and HRRR mlht at receptors ------------------------------------------
BL <- NULL; bl_f <- file.path(DATA_BASE, "outputs", "blh_timeseries.csv")
if (file.exists(bl_f)) { BL <- read.csv(bl_f); BL$t <- as.POSIXct(BL$time, tz = "UTC")
  BL <- BL[is.finite(BL$blh_m) & BL$blh_m > 150 & BL$blh_m < 4500 & is.finite(BL$t), ] } else message("no lidar file ", bl_f)
jf <- list.files(REF_DIR, pattern = "^jacobian_.*\\.rds$", full.names = TRUE, recursive = FALSE)
HR <- NULL
if (length(jf)) { J <- readRDS(jf[1]); S0 <- read.csv(file.path(INV_OUT, "obs_segments.csv"), stringsAsFactors = FALSE)
  HR <- data.frame(id = J$ids, mlht = J$mlht_receptor)
  HR <- merge(HR, S0[, c("id", "flight", "time_utc", "lat", "lon")], by = "id")
  HR$t <- as.POSIXct(HR$time_utc, tz = "UTC"); HR <- HR[is.finite(HR$mlht), ]
} else message("no jacobian in ", REF_DIR, ": HRRR comparison skipped")

## ---- profiles -----------------------------------------------------------------------------------------
files <- list_flights(DATA_DIR); files <- files[grepl("ARL-Suite", files)]
P <- list(); PROF <- list()
for (p in files) {
  ic <- tryCatch(read_icartt(p), error = function(e) NULL); if (is.null(ic)) next
  d <- ic$data; fl <- sub("AMMBEC-ARL-Suite_TwinOtter_", "", sub("\\.ict$", "", basename(p)))
  need <- c("ALTGPS", "ALTAGL", "Temp_True", "Alt_Pressure", "timestamp")
  if (!all(need %in% names(d))) { message("skip ", fl, ": missing ", paste(setdiff(need, names(d)), collapse = ",")); next }
  d <- d[order(d$timestamp), ]; n <- nrow(d)
  tt <- as.numeric(d$timestamp); z <- d$ALTGPS
  vs <- c(NA, diff(z) / pmax(1, diff(tt)))
  vs <- as.numeric(stats::filter(ifelse(is.finite(vs), vs, 0), rep(1 / 15, 15), sides = 2))
  dirn <- ifelse(is.finite(vs) & vs > VS_MIN, 1L, ifelse(is.finite(vs) & vs < -VS_MIN, -1L, 0L))
  r <- rle(dirn); ends <- cumsum(r$lengths); starts <- ends - r$lengths + 1
  gnd <- median(d$ALTGPS - d$ALTAGL, na.rm = TRUE)
  for (k in which(r$values != 0)) {
    i <- starts[k]:ends[k]; s <- d[i, ]
    agl <- ifelse(is.finite(s$ALTAGL), s$ALTAGL, s$ALTGPS - gnd)
    if (sum(is.finite(agl)) < 20 || diff(tt[range(i)]) < MIN_DUR_S) next
    span <- diff(range(agl, na.rm = TRUE)); bot <- min(agl, na.rm = TRUE); top <- max(agl, na.rm = TRUE)
    if (span < SPAN_MIN || bot > BOTTOM_MAX) next
    Tk <- s$Temp_True + ifelse(median(s$Temp_True, na.rm = TRUE) < 150, 273.15, 0)
    th <- Tk * (1000 / p_from_palt(s$Alt_Pressure))^0.2857
    b <- floor(agl / BIN_M) * BIN_M + BIN_M / 2
    agg <- function(v) if (!is.null(v)) tapply(v, b, function(x) if (sum(is.finite(x)) >= 2) mean(x, na.rm = TRUE) else NA_real_) else NULL
    zb <- as.numeric(names(agg(th))); TH <- as.numeric(agg(th))
    Q <- if ("H2O_ppm" %in% names(s)) as.numeric(agg(s$H2O_ppm)) else rep(NA_real_, length(zb))
    C <- if ("CH4_ppb" %in% names(s)) as.numeric(agg(s$CH4_ppb)) else rep(NA_real_, length(zb))
    ok <- is.finite(TH); if (sum(ok) < 6) next
    zb <- zb[ok]; TH <- TH[ok]; Q <- Q[ok]; C <- C[ok]
    # Heffter (1980) inversion test on the 3-bin running-mean theta: the lowest layer where
    # dtheta/dz >= GRAD_MIN and theta rises by at least DTHETA across the layer; zi = layer base
    ths <- as.numeric(stats::filter(TH, rep(1 / 3, 3), sides = 2)); ths[is.na(ths)] <- TH[is.na(ths)]
    g <- diff(ths) / diff(zb); st <- is.finite(g) & g >= GRAD_MIN
    rr <- rle(st); e <- cumsum(rr$lengths); b0 <- e - rr$lengths + 1; zi_theta <- NA_real_
    for (m in which(rr$values)) { lo <- b0[m]; hi <- e[m] + 1
      if (zb[lo] < bot + SURF_MARGIN) next                       # surface layer, not the mixed-layer top
      if (ths[hi] - ths[lo] >= DTHETA) { zi_theta <- zb[lo]; break } }
    grad_min_z <- function(v) { g <- diff(v) / diff(zb); zm <- (zb[-1] + zb[-length(zb)]) / 2
      k <- is.finite(g) & zm > bot + 100 & zm < top - 50; if (sum(k) < 3) return(NA_real_); zm[k][which.min(g[k])] }
    zi_h2o <- grad_min_z(Q); zi_ch4 <- grad_min_z(C)
    tm <- d$timestamp[round(mean(range(i)))]; lat <- median(s$Latitude, na.rm = TRUE); lon <- median(s$Longitude, na.rm = TRUE)
    lid <- if (!is.null(BL)) { kk <- abs(as.numeric(difftime(BL$t, tm, units = "mins"))) <= 30; if (any(kk)) median(BL$blh_m[kk]) else NA_real_ } else NA_real_
    hr <- if (!is.null(HR)) { kk <- HR$flight == fl & abs(as.numeric(difftime(HR$t, tm, units = "mins"))) <= 60 &
      hav(HR$lat, HR$lon, lat, lon) <= 30; if (any(kk)) median(HR$mlht[kk]) else NA_real_ } else NA_real_
    P[[length(P) + 1]] <- data.frame(flight = fl, profile = length(P) + 1, t_utc = format(tm, "%Y-%m-%d %H:%M", tz = "UTC"),
      hour_mdt = as.integer(format(tm - 6 * 3600, "%H", tz = "UTC")), direction = ifelse(r$values[k] > 0, "climb", "descent"),
      lat = round(lat, 3), lon = round(lon, 3), bottom_agl = round(bot), top_agl = round(top),
      zi_theta = zi_theta, above_top = is.na(zi_theta), zi_lower_bound = ifelse(is.na(zi_theta), round(top), NA),
      zi_h2o = zi_h2o, zi_ch4 = zi_ch4, lidar_blh = round(lid), hrrr_mlht = round(hr), stringsAsFactors = FALSE)
    PROF[[length(PROF) + 1]] <- list(id = sprintf("%s #%d %s", fl, length(P), format(tm, "%H:%M", tz = "UTC")),
                                     z = zb, th = TH, q = Q, zi = zi_theta, lid = lid, hr = hr)
  }
  message(sprintf("%-16s %d profiles so far", fl, length(P)))
}
if (!length(P)) stop("no vertical profiles found (relax VS_MIN / SPAN_MIN / BOTTOM_MAX)")
PR <- do.call(rbind, P)
write.csv(PR, file.path(REF_DIR, "blh_aircraft_profiles.csv"), row.names = FALSE)

## ---- comparison -------------------------------------------------------------------------------------
sink(file.path(REF_DIR, "blh_aircraft_summary.txt"), split = TRUE)
day <- PR$hour_mdt >= 10 & PR$hour_mdt <= 17
cat(sprintf("aircraft profiles: %d (%d daytime 10-17 MDT); inversion found below the top in %d, censored (mixed layer deeper than the profile) in %d\n",
            nrow(PR), sum(day), sum(!PR$above_top), sum(PR$above_top)))
cat(sprintf("profile tops: median %.0f m AGL (range %.0f-%.0f)\n", median(PR$top_agl), min(PR$top_agl), max(PR$top_agl)))
q <- function(x) sprintf("%.2f (IQR %.2f-%.2f, n=%d)", median(x), quantile(x, .25), quantile(x, .75), length(x))
D <- PR[day, ]
a <- D[!D$above_top & is.finite(D$lidar_blh), ]; if (nrow(a)) cat("lidar / aircraft theta BLH:  ", q(a$lidar_blh / a$zi_theta), "\n")
a <- D[!D$above_top & is.finite(D$hrrr_mlht), ]; if (nrow(a)) cat("HRRR  / aircraft theta BLH:  ", q(a$hrrr_mlht / a$zi_theta), "\n")
a <- D[is.finite(D$lidar_blh) & is.finite(D$hrrr_mlht), ]; if (nrow(a)) cat("HRRR  / lidar (same times):  ", q(a$hrrr_mlht / a$lidar_blh), "\n")
a <- D[D$above_top & is.finite(D$lidar_blh), ]
if (nrow(a)) cat(sprintf("censored profiles: lidar BLH below the aircraft lower bound in %d of %d (evidence the lidar reads low there)\n",
                         sum(a$lidar_blh < a$zi_lower_bound), nrow(a)))
a <- D[D$above_top & is.finite(D$hrrr_mlht), ]
if (nrow(a)) cat(sprintf("censored profiles: HRRR mlht below the aircraft lower bound in %d of %d\n", sum(a$hrrr_mlht < a$zi_lower_bound), nrow(a)))
a <- D[!D$above_top & is.finite(D$zi_h2o), ]; if (nrow(a)) cat("H2O-gradient / theta BLH:    ", q(a$zi_h2o / a$zi_theta), "  (consistency of the two aircraft methods)\n")
a <- D[!D$above_top & is.finite(D$hrrr_mlht), ]
if (nrow(a) >= 3) { z <- median(a$zi_theta / a$hrrr_mlht)
  cat(sprintf("\nimplied ZISCALE (aircraft / HRRR, daytime, uncensored): %.2f   compare 05b lidar-based 0.37\n", z))
  cat("  Censored profiles mean the aircraft value is a lower bound: the true ZISCALE is at least this large where they dominate.\n") }
# Censoring-aware estimate: treat aircraft/HRRR as right-censored at bound/HRRR for profiles whose mixed-layer
# top lay above the profile, and estimate the distribution of the ratio with the Kaplan-Meier product-limit
# estimator (daytime profiles with an HRRR match). This is what the uncensored median cannot give: the
# population median and the probability that the ratio exceeds each candidate ZISCALE.
r_all <- D[is.finite(D$hrrr_mlht) & (is.finite(D$zi_theta) | is.finite(D$zi_lower_bound)), ]
if (nrow(r_all) >= 5) {
  ratio <- ifelse(r_all$above_top, r_all$zi_lower_bound, r_all$zi_theta) / r_all$hrrr_mlht; event <- !r_all$above_top
  o <- order(ratio); ratio <- ratio[o]; event <- event[o]; S <- 1; surv <- numeric(length(ratio)); at_risk <- length(ratio)
  for (i in seq_along(ratio)) { if (event[i]) S <- S * (1 - 1 / at_risk); surv[i] <- S; at_risk <- at_risk - 1 }
  km_q <- function(p) { i <- which(surv <= 1 - p)[1]; if (is.na(i)) NA else ratio[i] }
  km_S <- function(x) { i <- which(ratio <= x); if (length(i)) surv[max(i)] else 1 }
  cat(sprintf("\nKaplan-Meier (censoring-aware) aircraft / HRRR ratio, daytime, n = %d (%d resolved, %d censored):\n", nrow(r_all), sum(event), sum(!event)))
  cat(sprintf("  median %.2f, IQR %.2f-%.2f\n", km_q(0.5), km_q(0.25), km_q(0.75)))
  for (x in c(0.37, 0.5, 0.65, 0.8, 1.0, 1.2)) cat(sprintf("  P(ratio > %.2f) = %.2f\n", x, km_S(x)))
  # Bootstrap (profiles resampled with replacement, 2000 draws) for the uncertainty of the KM median, quartiles and
  # exceedance fractions. These describe the DISTRIBUTION of the profile ratio; they are not a test of a single
  # scaling factor. KM assumes the censoring value (profile top / HRRR) is independent of the true ratio, which
  # holds only approximately here because both share the HRRR denominator: stated as an assumption in the paper.
  km_fit <- function(rt, ev) { o <- order(rt); rt <- rt[o]; ev <- ev[o]; S <- 1; sv <- numeric(length(rt)); ar <- length(rt)
    for (i in seq_along(rt)) { if (ev[i]) S <- S * (1 - 1 / ar); sv[i] <- S; ar <- ar - 1 }
    list(q = function(p) { i <- which(sv <= 1 - p)[1]; if (is.na(i)) NA else rt[i] }, S = function(x) { i <- which(rt <= x); if (length(i)) sv[max(i)] else 1 }) }
  set.seed(1); B <- 2000; xs <- c(0.37, 0.5, 0.65, 0.8, 1.0, 1.2)
  bs <- t(replicate(B, { i <- sample(length(ratio), replace = TRUE); k <- km_fit(ratio[i], event[i]); c(k$q(0.5), k$q(0.25), k$q(0.75), sapply(xs, k$S)) }))
  ci <- apply(bs, 2, quantile, c(0.025, 0.975), na.rm = TRUE)
  cat(sprintf("  bootstrap 95%% CI: median %.2f-%.2f, q25 %.2f-%.2f, q75 %.2f-%.2f\n", ci[1, 1], ci[2, 1], ci[1, 2], ci[2, 2], ci[1, 3], ci[2, 3]))
  for (j in seq_along(xs)) cat(sprintf("  P(ratio > %.2f): 95%% CI %.2f-%.2f\n", xs[j], ci[1, 3 + j], ci[2, 3 + j]))
  write.csv(data.frame(ziscale = xs, p_ratio_exceeds = sapply(xs, km_S), p_lo95 = ci[1, 4:9], p_hi95 = ci[2, 4:9],
                       km_median = km_q(0.5), km_median_lo95 = ci[1, 1], km_median_hi95 = ci[2, 1], km_q25 = km_q(0.25), km_q75 = km_q(0.75),
                       n = nrow(r_all), n_resolved = sum(event)),
            file.path(REF_DIR, "blh_aircraft_km.csv"), row.names = FALSE)
}
sink()
print(PR[, c("flight", "t_utc", "direction", "bottom_agl", "top_agl", "zi_theta", "zi_lower_bound", "zi_h2o", "lidar_blh", "hrrr_mlht")], row.names = FALSE)

## ---- figures ----------------------------------------------------------------------------------------
png(file.path(REF_DIR, "figures", "blh_aircraft.png"), 1800, 650, res = 150); par(mfrow = c(1, 3), mar = c(4.5, 4.5, 2.5, 1))
lim <- c(0, max(c(PR$zi_theta, PR$zi_lower_bound, PR$lidar_blh, PR$hrrr_mlht), na.rm = TRUE) * 1.05)
pc <- ifelse(PR$above_top, 2, 19)
plot(pmax(PR$zi_theta, PR$zi_lower_bound, na.rm = TRUE), PR$lidar_blh, pch = pc, xlim = lim, ylim = lim,
     xlab = "aircraft BLH (m AGL; open = lower bound)", ylab = "lidar BLH (m)", main = "Lidar vs aircraft"); abline(0, 1, lty = 2)
plot(pmax(PR$zi_theta, PR$zi_lower_bound, na.rm = TRUE), PR$hrrr_mlht, pch = pc, xlim = lim, ylim = lim,
     xlab = "aircraft BLH (m AGL; open = lower bound)", ylab = "HRRR mlht at receptors (m)", main = "HRRR vs aircraft"); abline(0, 1, lty = 2)
plot(PR$lidar_blh, PR$hrrr_mlht, pch = 19, xlim = lim, ylim = lim, xlab = "lidar BLH (m)", ylab = "HRRR mlht (m)", main = "HRRR vs lidar")
abline(0, 1, lty = 2); abline(0, 2.69, lty = 3, col = "grey50")
dev.off()
np <- min(length(PROF), 16); nc <- 4; nr <- ceiling(np / nc)
png(file.path(REF_DIR, "figures", "blh_aircraft_profiles.png"), 420 * nc, 380 * nr, res = 120); par(mfrow = c(nr, nc), mar = c(4, 4, 2, 1))
for (z in PROF[seq_len(np)]) {
  plot(z$th, z$z, type = "b", pch = 20, cex = 0.6, xlab = "theta (K)", ylab = "m AGL", main = z$id, cex.main = 0.85,
       ylim = c(0, max(c(z$z, z$lid, z$hr), na.rm = TRUE)))
  hv <- c(z$zi, z$lid, z$hr); k <- is.finite(hv)
  if (any(k)) abline(h = hv[k], col = c("black", "blue", "red")[k], lty = (1:3)[k])
}
legend("bottomright", bty = "n", cex = 0.7, lty = 1:3, col = c("black", "blue", "red"), legend = c("aircraft theta BLH", "lidar", "HRRR"))
dev.off()
cat("wrote blh_aircraft_profiles.csv, blh_aircraft_summary.txt, figures/blh_aircraft*.png to", REF_DIR, "\n")
