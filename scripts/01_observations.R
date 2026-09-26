# 01_observations.R — the observation vector: CH4, C2H6 and CO on every level leg ---------------
# All Twin Otter flights with CH4, C2H6 and CO (the ARL-Suite files), all level legs (Paper 1 leg
# detection), inside OBS_BOX (Denver + the DJ Basin) and below MAX_AGL_M, cut into segments of about SEG_S
# seconds within each leg (equal parts, so no tail is dropped). Unlike Paper 1 nothing is gated on enhancement: low values constrain
# the model as much as plumes do. Each segment is one STILT receptor (its mid-point).
#
# Written per segment: flight, leg, region (urban / edge / basin, Paper 1 definition), time,
# position, height, mean and sd of the three species, the number of valid 1-Hz samples, the
# aircraft wind, and Paper 1's rolling-baseline enhancements (diagnostic only; the inversion
# estimates its own backgrounds).
#
#   Rscript scripts/01_observations.R
#   Out: <INV_OUT>/obs_segments.csv, obs_flights.csv
# ----------------------------------------------------------------------------------------------
proj <- if (file.exists("config.R")) "." else ".."; source(file.path(proj, "config.R"))

vmean_wd <- function(wd, ws) {
  ok <- is.finite(wd) & is.finite(ws); if (!any(ok)) return(c(NA, NA))
  u <- -ws[ok] * sin(wd[ok] * pi / 180); v <- -ws[ok] * cos(wd[ok] * pi / 180)
  c((atan2(-mean(u), -mean(v)) * 180 / pi + 360) %% 360, sqrt(mean(u)^2 + mean(v)^2))
}
files <- list_flights(DATA_DIR); files <- files[grepl("ARL-Suite", files)]
if (!length(files)) stop("no ARL-Suite ICARTT files under ", file.path(DATA_DIR, "Aircraft"))

segs <- list(); fl_rows <- list()
for (p in files) {
  ic <- tryCatch(read_icartt(p), error = function(e) NULL)
  if (is.null(ic) || !all(c("CH4_ppb", "C2H6_ppb", "CO_ppb") %in% names(ic$data))) { message("skip ", basename(p)); next }
  fl <- sub("AMMBEC-ARL-Suite_TwinOtter_", "", sub("\\.ict$", "", basename(p)))
  d <- detect_level_legs(ic$data)
  for (s in c("CH4_ppb", "C2H6_ppb", "CO_ppb")) d <- add_enhancements(d, s)
  inbox <- d$Latitude >= OBS_BOX["ymn"] & d$Latitude <= OBS_BOX["ymx"] & d$Longitude >= OBS_BOX["xmn"] &
           d$Longitude <= OBS_BOX["xmx"] & is.finite(d$ALTAGL) & d$ALTAGL <= MAX_AGL_M
  use <- d$leg_id > 0 & inbox & is.finite(d$CH4_ppb)
  fl_rows[[length(fl_rows) + 1]] <- data.frame(flight = fl, date = as.character(ic$meta$date),
    t_start = format(min(d$timestamp), "%H:%M", tz = "UTC"), t_end = format(max(d$timestamp), "%H:%M", tz = "UTC"),
    n_legs = length(unique(d$leg_id[use])), n_samples = sum(use),
    ch4_p05 = quantile(d$CH4_ppb[use], 0.05, na.rm = TRUE), c2h6_p05 = quantile(d$C2H6_ppb[use], 0.05, na.rm = TRUE),
    co_p05 = quantile(d$CO_ppb[use], 0.05, na.rm = TRUE), stringsAsFactors = FALSE)
  if (!any(use)) next
  for (g in sort(unique(d$leg_id[use]))) {
    s <- d[use & d$leg_id == g, ]
    # split the leg into equal parts of about SEG_S seconds (at least one), so no tail is lost
    npart <- max(1L, round(nrow(s) / SEG_S)); k <- floor((seq_len(nrow(s)) - 1) * npart / nrow(s))
    for (kk in unique(k)) {
      q <- s[k == kk, ]
      if (sum(is.finite(q$CH4_ppb)) < MIN_SEG_N) next
      w <- vmean_wd(q$WD, q$WS)
      mean_or_na <- function(x) if (sum(is.finite(x)) >= MIN_SEG_N / 2) mean(x, na.rm = TRUE) else NA_real_
      sd_or_na   <- function(x) if (sum(is.finite(x)) >= MIN_SEG_N / 2) sd(x, na.rm = TRUE) else NA_real_
      segs[[length(segs) + 1]] <- data.frame(flight = fl, leg_id = g, seg = kk,
        time_utc = format(as.POSIXct(mean(as.numeric(q$timestamp)), origin = "1970-01-01", tz = "UTC"), "%Y-%m-%d %H:%M:%S"),
        lat = mean(q$Latitude), lon = mean(q$Longitude), agl_m = mean(q$ALTAGL, na.rm = TRUE), alt_m = mean(q$ALTGPS, na.rm = TRUE),
        n = sum(is.finite(q$CH4_ppb)),
        ch4 = mean(q$CH4_ppb, na.rm = TRUE), ch4_sd = sd(q$CH4_ppb, na.rm = TRUE),
        c2h6 = mean_or_na(q$C2H6_ppb), c2h6_sd = sd_or_na(q$C2H6_ppb),
        co = mean_or_na(q$CO_ppb), co_sd = sd_or_na(q$CO_ppb),
        ch4_enh_p1 = mean(q$CH4_ppb_enh, na.rm = TRUE), c2h6_enh_p1 = mean(q$C2H6_ppb_enh, na.rm = TRUE),
        co_enh_p1 = mean(q$CO_ppb_enh, na.rm = TRUE),
        wd_deg = w[1], ws_ms = w[2], stringsAsFactors = FALSE)
    }
  }
  message(sprintf("%-16s %3d legs, %4d segments", fl, length(unique(d$leg_id[use])), sum(vapply(segs, function(x) x$flight[1] == fl, NA))))
}
S <- do.call(rbind, segs)
S$region <- tag_region(data.frame(lat = S$lat, lon = S$lon), box = URBAN_BOX, basin_lat = BASIN_LAT)$region
S$id <- sprintf("%s_L%02d_s%03d", S$flight, S$leg_id, S$seg)
S <- S[, c("id", setdiff(names(S), "id"))]
num <- vapply(S, is.numeric, NA); S[num] <- lapply(S[num], function(x) round(x, 4))
write.csv(S, file.path(INV_OUT, "obs_segments.csv"), row.names = FALSE)
FL <- do.call(rbind, fl_rows); rownames(FL) <- NULL
write.csv(FL, file.path(INV_OUT, "obs_flights.csv"), row.names = FALSE)
cat(sprintf("\n%d segments of %d s from %d flights (%d urban, %d edge, %d basin); C2H6 on %d, CO on %d\n",
            nrow(S), SEG_S, length(unique(S$flight)), sum(S$region == "urban"), sum(S$region == "edge"),
            sum(S$region == "basin"), sum(is.finite(S$c2h6)), sum(is.finite(S$co))))
cat("wrote obs_segments.csv, obs_flights.csv to ", INV_OUT, "\n", sep = "")
