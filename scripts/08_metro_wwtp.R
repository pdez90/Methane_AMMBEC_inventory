# 08_metro_wwtp.R — what the campaign data say about the Metro Water Recovery plant --------------
# The inversion hardly sees the plant: in every inventory it is small against the aircraft noise
# and it sits 1 km from the Suncor refinery, so its scale factor stays near the prior. This
# script collects the direct evidence at plant scale, separately from the inversion:
#
#  (A) GROUND MOBILE SURVEYS (Paper 1 outputs/mobile_hotspots.csv, CDPHE CAT/EMU 1 Hz CH4).
#      200-m grid around the plant: samples, days, median and 90th-percentile enhancement, share
#      of samples > 100 ppb, distance to each point source; coverage per source (closest
#      approach, samples within 0.5/1/2 km); one row per survey day. These surveys carry no wind
#      and no ethane, so they show WHERE methane is persistently high, not who emits it.
#      (MOOSE MeFTIR, which does carry ethane, never came within 3 km of the plant: checked 26 Sep.)
#  (B) AIRCRAFT 1 Hz TRANSECTS. Every level leg (AGL < 1500 m) that passes downwind of any of the
#      corridor sources (within D_MAX_KM, sample bearing within ALIGN_DEG of the wind) gives a
#      single-transect mass balance (as Paper 1 script 59: L = sum n_air dX u_perp dy, Q = L x BLH,
#      lidar BLH; the rate scales with BLH, so Q at 0.8 and 1.2 x BLH is also given). The crossing
#      is then split laterally: each source is advected along the leg-mean wind to its predicted
#      crossing point on the leg, and every sample is assigned to the nearest predicted crossing.
#      Each window gets its own flux and ethane slope, so the plant (no ethane) and the refinery
#      (fossil) are separated where the geometry allows. "resolved" = the plant and Suncor
#      crossings are further apart than the lateral error from a WD_UNC_DEG wind-direction error.
#      Biogenic rate = (1 - fossil fraction) x Q, fossil fraction = York slope / 0.102 (Paper 1).
#  (C) CARBON MAPPER: CH4 and CO2 plumes within CM_KM of each source.
#  (D) SUMMARY against the plant component of each prior (PRIOR_DIR, METRO_COMPONENT = "wwtp"),
#      and against Paper 1's complex encounters (outputs/v2_pointsource_encounters.csv).
#
# The plant location is POINT_SOURCES in config.R (APPROXIMATE, override with METHANE_METRO_LAT /
# METHANE_METRO_LON). Everything here is diagnostic; none of it feeds the inversion.
#
#   Rscript scripts/08_metro_wwtp.R
#   Out: <INV_OUT>/metro_wwtp/mobile_grid.csv, mobile_coverage.csv, mobile_days.csv,
#        aircraft_crossings.csv, aircraft_vs_paper1.csv, carbonmapper_near.csv, summary.txt,
#        figures/mobile_map.png, figures/aircraft_crossings.png
# -----------------------------------------------------------------------------------------------
source(file.path(if (file.exists("config.R")) "." else "..", "config.R"))
OUT <- file.path(INV_OUT, "metro_wwtp"); FIG <- file.path(OUT, "figures"); dir.create(FIG, FALSE, TRUE)
PS <- POINT_SOURCES; iM <- which(PS$name == "metro_wwtp"); iS <- which(PS$name == "suncor")

MOB_KM <- 5; CELL_KM <- 0.2                      # (A) window and grid
D_MAX_KM <- 12; ALIGN_DEG <- 25; AGL_MAX_M <- 1500 # (B) same criterion as Paper 1 script 59
MIN_PTS <- 10; MIN_ENH_PPB <- 20; WD_UNC_DEG <- 15; BETA0 <- 0.102
CM_KM <- 5                                       # (C)

hav <- function(la1, lo1, la2, lo2) { p <- pi / 180
  a <- sin((la2 - la1) * p / 2)^2 + cos(la1 * p) * cos(la2 * p) * sin((lo2 - lo1) * p / 2)^2
  2 * 6371.0088 * asin(pmin(1, sqrt(a))) }
bearing <- function(la1, lo1, la2, lo2) { p <- pi / 180
  y <- sin((lo2 - lo1) * p) * cos(la2 * p); x <- cos(la1 * p) * sin(la2 * p) - sin(la1 * p) * cos(la2 * p) * cos((lo2 - lo1) * p)
  (atan2(y, x) / p + 360) %% 360 }
ang_diff <- function(a, b) abs((a - b + 180) %% 360 - 180)
to_xy <- function(lat, lon, lat0, lon0) cbind(x = (lon - lon0) * 111.32 * cos(lat0 * pi / 180), y = (lat - lat0) * 110.57)
q <- function(x, p) if (sum(is.finite(x))) unname(quantile(x, p, na.rm = TRUE)) else NA_real_
wd_vmean <- function(wd, ws) { u <- -ws * sin(wd * pi / 180); v <- -ws * cos(wd * pi / 180)
  (atan2(-mean(u), -mean(v)) * 180 / pi + 360) %% 360 }

cat("point sources (config POINT_SOURCES):\n"); print(PS, row.names = FALSE)
cat(sprintf("Metro-Suncor separation %.2f km, Metro-Cherokee %.2f km\n\n",
            hav(PS$lat[iM], PS$lon[iM], PS$lat[iS], PS$lon[iS]),
            hav(PS$lat[iM], PS$lon[iM], PS$lat[PS$name == "cherokee"], PS$lon[PS$name == "cherokee"])))

## ---- (A) ground mobile surveys --------------------------------------------------------------------
mf <- file.path(DATA_BASE, "outputs", "mobile_hotspots.csv"); MG <- MC <- MD <- NULL
if (file.exists(mf)) {
  M <- read.csv(mf, stringsAsFactors = FALSE)
  M$d_metro <- hav(M$Latitude, M$Longitude, PS$lat[iM], PS$lon[iM])
  M <- M[is.finite(M$d_metro) & M$d_metro <= MOB_KM & is.finite(M$CH4_ppmv_enh), ]
  M$enh <- M$CH4_ppmv_enh * 1000; M$day <- substr(M$timestamp, 1, 10)
  cat(sprintf("(A) mobile: %d samples on %d days within %.0f km of the plant (sites: %s)\n", nrow(M),
              length(unique(M$day)), MOB_KM, paste(unique(M$site), collapse = ", ")))
  xy <- to_xy(M$Latitude, M$Longitude, PS$lat[iM], PS$lon[iM])
  M$cell <- paste(floor(xy[, 1] / CELL_KM), floor(xy[, 2] / CELL_KM))
  MG <- do.call(rbind, lapply(split(M, M$cell), function(s) data.frame(
    lat = mean(s$Latitude), lon = mean(s$Longitude), n = nrow(s), n_days = length(unique(s$day)),
    p50_ppb = q(s$enh, 0.5), p90_ppb = q(s$enh, 0.9), max_ppb = max(s$enh), frac_gt100 = mean(s$enh > 100))))
  for (j in seq_len(nrow(PS))) MG[[paste0("km_", PS$name[j])]] <- round(hav(MG$lat, MG$lon, PS$lat[j], PS$lon[j]), 2)
  dk <- as.matrix(MG[, paste0("km_", PS$name)]); MG$nearest <- PS$name[apply(dk, 1, which.min)]
  MG <- MG[order(-MG$p90_ppb), ]; num <- vapply(MG, is.numeric, NA); MG[num] <- lapply(MG[num], round, 4)
  write.csv(MG, file.path(OUT, "mobile_grid.csv"), row.names = FALSE)
  MC <- do.call(rbind, lapply(seq_len(nrow(PS)), function(j) {
    d <- hav(M$Latitude, M$Longitude, PS$lat[j], PS$lon[j])
    data.frame(source = PS$name[j], closest_km = round(min(d), 2), n_0.5km = sum(d <= 0.5), n_1km = sum(d <= 1),
               n_2km = sum(d <= 2), days_1km = length(unique(M$day[d <= 1])),
               p90_1km_ppb = round(q(M$enh[d <= 1], 0.9)), p90_1to2km_ppb = round(q(M$enh[d > 1 & d <= 2], 0.9)))
  }))
  write.csv(MC, file.path(OUT, "mobile_coverage.csv"), row.names = FALSE)
  MD <- do.call(rbind, lapply(split(M, paste(M$day, M$site)), function(s) { k <- which.max(s$enh)
    dmx <- hav(s$Latitude[k], s$Longitude[k], PS$lat, PS$lon)
    data.frame(day = s$day[1], site = s$site[1], n = nrow(s), closest_metro_km = round(min(s$d_metro), 2),
               p90_ppb = round(q(s$enh, 0.9)), max_ppb = round(s$enh[k]), max_lat = s$Latitude[k], max_lon = s$Longitude[k],
               max_nearest = PS$name[which.min(dmx)], max_nearest_km = round(min(dmx), 2)) }))
  MD <- MD[order(MD$day), ]
  write.csv(MD, file.path(OUT, "mobile_days.csv"), row.names = FALSE)
  cat("coverage per source:\n"); print(MC, row.names = FALSE)
  cat("highest 200-m cells (p90):\n"); print(head(MG[, c("lat", "lon", "n", "n_days", "p50_ppb", "p90_ppb", "frac_gt100", "nearest",
                                                       "km_metro_wwtp", "km_suncor")], 12), row.names = FALSE)
  # map
  png(file.path(FIG, "mobile_map.png"), 1500, 1300, res = 160); par(mar = c(4.5, 4.5, 3, 6))
  br <- c(-Inf, 10, 25, 50, 100, 250, 500, 1000, Inf); pal <- hcl.colors(length(br) - 1, "YlOrRd", rev = TRUE)
  cl <- pal[cut(MG$p90_ppb, br)]
  plot(MG$lon, MG$lat, pch = 15, cex = 1.1, col = cl, asp = 1 / cos(39.8 * pi / 180), xlab = "longitude", ylab = "latitude",
       main = sprintf("Ground mobile CH4, 200-m cells within %d km of Metro: p90 enhancement", MOB_KM))
  th <- seq(0, 2 * pi, length = 100)
  for (r in c(1, 2)) lines(PS$lon[iM] + r / (111.32 * cos(PS$lat[iM] * pi / 180)) * cos(th), PS$lat[iM] + r / 110.57 * sin(th), lty = 3)
  points(PS$lon, PS$lat, pch = c(8, 17, 15), cex = 2, col = "black", lwd = 2)
  text(PS$lon, PS$lat, PS$name, pos = 3, cex = 0.8)
  par(xpd = TRUE); legend("right", inset = -0.2, bty = "n", fill = pal, legend = levels(cut(0, br)), title = "ppb", cex = 0.75)
  dev.off()
} else message("(A) ", mf, " not found; mobile analysis skipped")

## ---- (B) aircraft transects -----------------------------------------------------------------------
bl_f <- file.path(DATA_BASE, "outputs", "blh_timeseries.csv"); BL <- NULL
if (file.exists(bl_f)) { BL <- read.csv(bl_f); BL$t <- as.POSIXct(BL$time, tz = "UTC")
  BL <- BL[is.finite(BL$blh_m) & BL$blh_m > 150 & BL$blh_m < 4500 & is.finite(BL$t), ] }
blh_at <- function(t0, t1) { if (is.null(BL)) return(NA_real_)
  k <- BL$t >= t0 - 1800 & BL$t <= t1 + 1800; if (sum(k)) median(BL$blh_m[k]) else NA_real_ }

files <- list_flights(DATA_DIR); files <- files[grepl("ARL-Suite", files)]
cr <- list(); prof <- list()
for (p in files) {
  ic <- tryCatch(read_icartt(p), error = function(e) NULL)
  if (is.null(ic) || !all(c("CH4_ppb", "C2H6_ppb", "WS", "WD") %in% names(ic$data))) next
  fl <- sub("AMMBEC-ARL-Suite_TwinOtter_", "", sub("\\.ict$", "", basename(p)))
  d <- detect_level_legs(ic$data); if (!any(d$leg_id > 0)) next
  d <- add_enhancements(d, "CH4_ppb"); d <- add_enhancements(d, "C2H6_ppb")
  lg <- leg_metrics(d)
  for (g in sort(unique(d$leg_id[d$leg_id > 0]))) {
    s <- d[d$leg_id == g, ]
    s <- s[is.finite(s$WD) & is.finite(s$WS) & is.finite(s$ALTAGL) & is.finite(s$CH4_ppb_enh) & is.finite(s$Latitude), ]
    if (nrow(s) < MIN_PTS || mean(s$ALTAGL) > AGL_MAX_M) next
    if (min(hav(s$Latitude, s$Longitude, PS$lat[iM], PS$lon[iM])) > D_MAX_KM) next
    link <- matrix(FALSE, nrow(s), nrow(PS))
    for (j in seq_len(nrow(PS))) {
      dist <- hav(s$Latitude, s$Longitude, PS$lat[j], PS$lon[j])
      mis <- ang_diff(bearing(PS$lat[j], PS$lon[j], s$Latitude, s$Longitude), (s$WD + 180) %% 360)
      link[, j] <- dist <= D_MAX_KM & mis <= ALIGN_DEG
    }
    if (sum(link[, iM]) < MIN_PTS) next                         # the plant itself must be upwind of the leg
    any_l <- rowSums(link) > 0; seg <- s[min(which(any_l)):max(which(any_l)), ]
    if (max(seg$CH4_ppb_enh) < MIN_ENH_PPB) next
    track <- lg$track_deg[lg$leg_id == g]
    n <- nrow(seg)
    dy <- c(hav(seg$Latitude[-n], seg$Longitude[-n], seg$Latitude[-1], seg$Longitude[-1]) * 1000, 0)
    up <- perp_wind(seg$WS, seg$WD, track)
    nair <- air_molar_density(seg$ALTGPS, if ("Temp_True" %in% names(seg)) seg$Temp_True else NA)
    integ <- nair * pmax(0, seg$CH4_ppb_enh) * 1e-9 * abs(up)
    wts <- (c(0, dy[-n]) + dy) / 2                               # trapezoid weight of each sample (m)
    H <- blh_at(min(seg$timestamp), max(seg$timestamp))
    to_t <- function(L) L * H * MW[["CH4"]] * 3.6 / 1000         # mol m-1 s-1 x m -> t/h
    # geometry: along-wind / cross-wind coordinates about the segment centre
    wd <- wd_vmean(seg$WD, seg$WS); wto <- (wd + 180) %% 360; ux <- sin(wto * pi / 180); uy <- cos(wto * pi / 180)
    lat0 <- mean(seg$Latitude); lon0 <- mean(seg$Longitude)
    sxy <- to_xy(seg$Latitude, seg$Longitude, lat0, lon0); lat_s <- sxy[, 1] * (-uy) + sxy[, 2] * ux
    pxy <- to_xy(PS$lat, PS$lon, lat0, lon0)
    along <- -(pxy[, 1] * ux + pxy[, 2] * uy); lat_p <- pxy[, 1] * (-uy) + pxy[, 2] * ux
    upw <- along > 0 & along <= D_MAX_KM
    inview <- upw & lat_p >= min(lat_s) - 1 & lat_p <= max(lat_s) + 1
    own <- if (any(upw)) apply(abs(outer(lat_s, lat_p[upw], "-")), 1, which.min) else rep(NA_integer_, n)
    own <- which(upw)[own]
    slope_of <- function(k) { k <- k & seg$CH4_ppb_enh > MIN_ENH_PPB & is.finite(seg$C2H6_ppb_enh)
      if (sum(k) < MIN_PTS) return(NA_real_)
      tryCatch(york_slope(seg$CH4_ppb_enh[k], seg$C2H6_ppb_enh[k], 1, 0.2)$slope, error = function(e) NA_real_) }
    ffrac <- function(b) if (is.finite(b)) max(0, min(1, b / BETA0)) else NA_real_
    L_all <- sum(integ * wts, na.rm = TRUE); b_all <- slope_of(rep(TRUE, n)); f_all <- ffrac(b_all)
    pk <- seg$CH4_ppb_enh > MIN_ENH_PPB
    cen <- if (sum(pk) >= 3) sum(lat_s[pk] * seg$CH4_ppb_enh[pk]) / sum(seg$CH4_ppb_enh[pk]) else NA_real_
    sep <- abs(lat_p[iM] - lat_p[iS]); sig <- along[iM] * tan(WD_UNC_DEG * pi / 180)
    row <- data.frame(flight = fl, leg_id = g, t_utc = format(min(seg$timestamp), "%Y-%m-%d %H:%M", tz = "UTC"),
      n_seg = n, seg_km = round(sum(dy) / 1000, 1), agl_m = round(mean(seg$ALTAGL)), wd_deg = round(wd), ws_ms = round(mean(seg$WS), 1),
      u_perp_ms = round(mean(abs(up)), 1), blh_m = round(H), dch4_max_ppb = round(max(seg$CH4_ppb_enh)),
      metro_along_km = round(along[iM], 1), metro_suncor_sep_km = round(sep, 2), lateral_err_km = round(sig, 2),
      resolved = sep > sig, centroid_minus_metro_km = round(cen - lat_p[iM], 2), centroid_minus_suncor_km = round(cen - lat_p[iS], 2),
      Q_t_hr = round(to_t(L_all), 2), Q_lo = round(0.8 * to_t(L_all), 2), Q_hi = round(1.2 * to_t(L_all), 2),
      slope = round(b_all, 4), fossil_frac = round(f_all, 2), Q_bio_t_hr = round((1 - f_all) * to_t(L_all), 2),
      stringsAsFactors = FALSE)
    for (j in seq_len(nrow(PS))) {                                  # per-source lateral windows
      k <- !is.na(own) & own == j; nm <- PS$name[j]
      Lj <- sum((integ * wts)[k], na.rm = TRUE); bj <- slope_of(k); fj <- ffrac(bj)
      row[[paste0("inview_", nm)]] <- inview[j]
      row[[paste0("Q_", nm)]] <- round(to_t(Lj), 2)
      row[[paste0("slope_", nm)]] <- round(bj, 4)
      row[[paste0("Qbio_", nm)]] <- round((1 - fj) * to_t(Lj), 2)
    }
    cr[[length(cr) + 1]] <- row
    prof[[length(prof) + 1]] <- list(id = sprintf("%s L%d", fl, g), lat_s = lat_s, enh = seg$CH4_ppb_enh, c2 = seg$C2H6_ppb_enh,
                                     lat_p = lat_p, upw = upw, resolved = sep > sig)
  }
  message(sprintf("%-16s %d crossings so far", fl, length(cr)))
}
CR <- if (length(cr)) do.call(rbind, cr) else NULL
if (!is.null(CR)) {
  write.csv(CR, file.path(OUT, "aircraft_crossings.csv"), row.names = FALSE)
  cat(sprintf("\n(B) aircraft: %d leg crossings downwind of the plant (%d laterally resolved from Suncor, %d with lidar BLH)\n",
              nrow(CR), sum(CR$resolved), sum(is.finite(CR$blh_m))))
  print(CR[, c("flight", "leg_id", "agl_m", "wd_deg", "blh_m", "dch4_max_ppb", "metro_along_km", "metro_suncor_sep_km",
               "lateral_err_km", "resolved", "centroid_minus_metro_km", "Q_t_hr", "fossil_frac", "Q_bio_t_hr",
               "Q_metro_wwtp", "slope_metro_wwtp", "Qbio_metro_wwtp", "Q_suncor", "slope_suncor")], row.names = FALSE)
  # per-crossing profiles
  np <- min(length(prof), 16); nc <- 4; nr <- ceiling(np / nc)
  png(file.path(FIG, "aircraft_crossings.png"), 450 * nc, 330 * nr, res = 120); par(mfrow = c(nr, nc), mar = c(4, 4, 2.2, 1))
  ordp <- order(-vapply(prof, function(z) max(z$enh), 1))[seq_len(np)]
  for (z in prof[ordp]) {
    o <- order(z$lat_s)
    plot(z$lat_s[o], z$enh[o], type = "l", xlab = "cross-wind distance (km)", ylab = "CH4 enh (ppb)",
         main = paste0(z$id, if (z$resolved) "  [resolved]" else ""), cex.main = 0.9)
    if (any(is.finite(z$c2))) lines(z$lat_s[o], z$c2[o] / BETA0, col = "orange")
    abline(v = z$lat_p[z$upw], col = c("blue", "red", "grey40")[z$upw], lty = 2, lwd = 1.5)
  }
  legend("topright", bty = "n", cex = 0.7, lty = c(1, 1, 2, 2, 2), col = c("black", "orange", "blue", "red", "grey40"),
         legend = c("CH4", "C2H6 / 0.102", "Metro", "Suncor", "Cherokee"))
  dev.off()
  # compare with Paper 1's complex encounters on the same legs
  pe <- file.path(DATA_BASE, "outputs", "v2_pointsource_encounters.csv")
  if (file.exists(pe)) {
    P1 <- read.csv(pe, stringsAsFactors = FALSE); P1 <- P1[grepl("Metro", P1$hotspot), ]
    if (nrow(P1)) {
      CMP <- merge(CR[, c("flight", "leg_id", "Q_t_hr", "Q_bio_t_hr", "Qbio_metro_wwtp", "resolved")],
                   P1[, c("flight", "leg_id", "Q_t_hr", "Q_bio_t_hr", "hotspot_in_fetch")], by = c("flight", "leg_id"),
                   all = TRUE, suffixes = c("_here", "_paper1"))
      write.csv(CMP, file.path(OUT, "aircraft_vs_paper1.csv"), row.names = FALSE)
      cat("\nsame legs in Paper 1 script 59 (complex encounters):\n"); print(CMP, row.names = FALSE)
    }
  }
} else cat("\n(B) aircraft: no leg crossed downwind of the plant under the criterion\n")

## ---- (C) Carbon Mapper ---------------------------------------------------------------------------
cmf <- file.path(DATA_BASE, "outputs", "carbonmapper_denver_plumes.csv"); CMN <- NULL
if (file.exists(cmf)) {
  C <- read.csv(cmf, stringsAsFactors = FALSE)
  CMN <- do.call(rbind, lapply(seq_len(nrow(PS)), function(j) { d <- hav(C$lat, C$lon, PS$lat[j], PS$lon[j]); k <- d <= CM_KM
    if (!any(k)) return(NULL)
    data.frame(source = PS$name[j], km = round(d[k], 2), C[k, intersect(c("plume_id", "gas", "sector", "date", "platform", "kg_hr", "kg_hr_unc"), names(C))]) }))
  cat(sprintf("\n(C) Carbon Mapper: %d plumes in the file; within %d km of a corridor source:\n", nrow(C), CM_KM))
  if (is.null(CMN)) cat("   none\n") else print(CMN, row.names = FALSE)
  cat(sprintf("   CH4 plumes within %d km of the plant: %d\n", CM_KM, if (is.null(CMN)) 0 else sum(CMN$source == "metro_wwtp" & CMN$gas == "CH4")))
  if (!is.null(CMN)) write.csv(CMN, file.path(OUT, "carbonmapper_near.csv"), row.names = FALSE)
} else message("(C) ", cmf, " not found; Carbon Mapper check skipped")

## ---- (D) summary against the priors -------------------------------------------------------------
pri <- list.files(PRIOR_DIR, pattern = "^prior_.*\\.rds$", full.names = TRUE)
PT <- do.call(rbind, lapply(pri, function(f) { g <- readRDS(f); k <- intersect(c("metro_wwtp", "metro_complex"), names(g$box_t_hr))
  if (!length(k)) return(NULL); data.frame(prior = g$meta$name, component = k[1], t_hr = round(g$box_t_hr[[k[1]]], 3)) }))
sink(file.path(OUT, "summary.txt"), split = TRUE)
cat("\n==== Metro Water Recovery: summary ====\n")
cat(sprintf("plant location used: %.4f N %.4f W (APPROXIMATE; config POINT_SOURCES)\n", PS$lat[iM], -PS$lon[iM]))
if (!is.null(PT)) { cat("prior plant component (t CH4/h):\n"); print(PT, row.names = FALSE) } else
  cat("no priors with a metro component in", PRIOR_DIR, "(run 04 first)\n")
if (!is.null(MC)) cat(sprintf("mobile: closest approach to the plant %.2f km, %d samples within 1 km on %d days; p90 within 1 km %s ppb\n",
                              MC$closest_km[iM], MC$n_1km[iM], MC$days_1km[iM], MC$p90_1km_ppb[iM]))
if (!is.null(CR)) {
  bh <- is.finite(CR$Q_t_hr) & is.finite(CR$blh_m); rs <- bh & CR$resolved
  cat(sprintf("aircraft: %d crossings, %d with BLH; median Q (whole crossing) %.2f t/h [%.2f-%.2f], biogenic %.2f t/h\n",
              nrow(CR), sum(bh), median(CR$Q_t_hr[bh]), q(CR$Q_t_hr[bh], 0), q(CR$Q_t_hr[bh], 1), median(CR$Q_bio_t_hr[bh], na.rm = TRUE)))
  if (any(rs)) cat(sprintf("          %d laterally resolved from Suncor: plant-window Q median %.2f t/h, plant-window biogenic %.2f t/h\n",
                           sum(rs), median(CR$Q_metro_wwtp[rs]), median(CR$Qbio_metro_wwtp[rs], na.rm = TRUE)))
  else cat("          no crossing resolves the plant from Suncor at the assumed wind-direction error\n")
}
cat(sprintf("Carbon Mapper CH4 detections within %d km of the plant: %d\n", CM_KM,
            if (is.null(CMN)) 0 else sum(CMN$source == "metro_wwtp" & CMN$gas == "CH4")))
sink()
cat("\nwrote", OUT, "\n")
