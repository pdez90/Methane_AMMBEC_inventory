# 04_priors.R — prior emission fields for the inversion ------------------------------------------
# Reads each inventory on its native grid over DOMAIN and splits it into the inversion components
# (config.R): the three named waste facilities, the rest of waste, basin and urban oil and gas
# (GRA2PES OG split by the urban mask), post-meter (RES), agriculture, and other (total minus the
# sectors). CO is the total CO field, used only as a transport tracer.
#
#   GRA2PES v1.1 and v2.0beta: 4-km Lambert grid, hourly (24 UTC hours), vertical levels summed,
#     mol km-2 h-1 -> umol m-2 s-1. Sector files exist for weekdays only, so the CH4 components use
#     the weekday fields every day (methane has almost no weekly cycle); CO uses the weekday,
#     Saturday or Sunday total where that tree exists (v1.1 allspec), else weekday (v2.0beta).
#   EPA gridded GHGI 2020: 0.1-degree, time-invariant, molec cm-2 s-1 -> umol m-2 s-1; no CO.
#
#   Rscript scripts/04_priors.R [v1.1] [v2.0beta] [epa]      (default: all three)
#   Out: <PRIOR_DIR>/prior_<name>.rds, prior_box_totals.csv
# ----------------------------------------------------------------------------------------------
proj <- if (file.exists("config.R")) "." else ".."; source(file.path(proj, "config.R"))
source(file.path(proj, "R", "grids.R"))
suppressPackageStartupMessages(library(ncdf4))

which_priors <- commandArgs(TRUE); if (!length(which_priors)) which_priors <- c("v1.1", "v2.0beta", "epa")
MONTH <- Sys.getenv("METHANE_GRA2PES_MONTH", "202307")

## ---- GRA2PES ----------------------------------------------------------------------------------
gra_files <- function(dir, daytype) {
  f <- list.files(file.path(dir, daytype), pattern = "\\.nc$", full.names = TRUE)
  f[order(grepl("12to23Z", f))]                                   # 00to11Z first, then 12to23Z
}
# Hourly field [nx, ny, 24] (UTC hour 0..23) of one variable over the index window, levels summed.
gra_read <- function(files, var, win) {
  out <- array(0, c(win$nx, win$ny, 24)); got <- 0L
  for (fp in files) {
    nc <- nc_open(fp); on.exit(nc_close(nc), add = TRUE)
    if (!var %in% names(nc$var)) { nc_close(nc); on.exit(); next }
    u <- nc$var[[var]]$units; if (!grepl("mole", u)) stop(var, " in ", basename(fp), " has units ", u)
    tu <- ncatt_get(nc, "time", "units")$value; t0 <- as.POSIXct(sub("hours since ", "", tu), tz = "UTC")
    hrs <- as.integer(format(t0 + 3600 * ncvar_get(nc, "time"), "%H", tz = "UTC"))
    dims <- vapply(nc$var[[var]]$dim, function(d) d$name, "")
    st <- c(x = win$x0, y = win$y0, level = 1, time = 1)[dims]; ct <- c(x = win$nx, y = win$ny, level = -1, time = -1)[dims]
    a <- ncvar_get(nc, var, start = st, count = ct, collapse_degen = FALSE)
    names(dim(a)) <- dims
    a <- aperm(a, match(c("x", "y", "level", "time"), dims))
    a[!is.finite(a)] <- 0
    a <- apply(a, c(1, 2, 4), sum)                                     # sum levels -> [x, y, time]
    for (k in seq_along(hrs)) out[, , hrs[k] + 1] <- a[, , k]
    got <- got + length(hrs); nc_close(nc); on.exit()
  }
  if (got < 24) stop(var, ": only ", got, " hours found in ", paste(basename(files), collapse = ", "))
  out * GRA2PES_TO_UMOL
}
gra_window <- function(any_file) {
  nc <- nc_open(any_file); on.exit(nc_close(nc))
  x <- ncvar_get(nc, "x"); y <- ncvar_get(nc, "y")
  lat <- ncvar_get(nc, "lat"); lon <- ncvar_get(nc, "lon"); if (max(lon, na.rm = TRUE) > 180) lon <- lon - 360
  m <- lat >= DOMAIN["ymn"] & lat <= DOMAIN["ymx"] & lon >= DOMAIN["xmn"] & lon <= DOMAIN["xmx"]
  ix <- range(which(apply(m, 1, any))); iy <- range(which(apply(m, 2, any)))
  lc <- ncatt_get(nc, "lambert_conformal_conic")
  sp <- ncatt_get(nc, "lambert_conformal_conic", "standard_parallel")$value
  list(x0 = ix[1], nx = diff(ix) + 1, y0 = iy[1], ny = diff(iy) + 1,
       x = x[ix[1]:ix[2]], y = y[iy[1]:iy[2]], lat = lat[ix[1]:ix[2], iy[1]:iy[2]], lon = lon[ix[1]:ix[2], iy[1]:iy[2]],
       lcc = list(lat0 = ncatt_get(nc, "lambert_conformal_conic", "latitude_of_projection_origin")$value,
                  lon0 = ncatt_get(nc, "lambert_conformal_conic", "longitude_of_central_meridian")$value,
                  std1 = sp[1], std2 = sp[2], R = ncatt_get(nc, "lambert_conformal_conic", "earth_radius")$value, fe = 0, fn = 0))
}
build_gra2pes <- function(ver) {
  sec_dir <- file.path(GRA2PES_DIR, "sectors", ver)
  tot_ch4 <- if (ver == "v1.1") file.path(GRA2PES_DIR, "ch4only", MONTH) else file.path(GRA2PES_DIR, "v2total", MONTH)
  tot_co  <- if (ver == "v1.1") file.path(GRA2PES_DIR, "allspec", MONTH) else file.path(GRA2PES_DIR, "v2total", MONTH)
  ref <- gra_files(file.path(sec_dir, "OG", MONTH), "weekdy")[1]
  if (is.na(ref)) stop("no GRA2PES ", ver, " OG sector files under ", sec_dir)
  win <- gra_window(ref)
  # check the spherical LCC used for particles against the file's own coordinates
  chk <- lcc_forward(win$lat, win$lon, win$lcc)
  err <- max(abs(chk$x - matrix(win$x, win$nx, win$ny)), abs(chk$y - matrix(win$y, win$nx, win$ny, byrow = TRUE)))
  if (err > 200) stop("LCC check failed: max position error ", round(err), " m")
  message(sprintf("GRA2PES %s: window %d x %d cells, LCC check max error %.1f m", ver, win$nx, win$ny, err))
  sec <- list()
  for (s in c("WASTE", "OG", "RES", "AG")) sec[[s]] <- gra_read(gra_files(file.path(sec_dir, s, MONTH), "weekdy"), "HC01", win)
  sec$TOTAL <- gra_read(gra_files(tot_ch4, "weekdy"), "HC01", win)
  g <- list(type = "lcc", x = win$x, y = win$y, lcc = win$lcc, lat = win$lat, lon = win$lon, hours = 0:23,
            meta = list(name = paste0("gra2pes_", ver), month = MONTH, units = "umol m-2 s-1"))
  masks <- grid_masks(g); g$comps <- make_components(sec, masks)
  g[["co"]] <- list()
  for (dt in c("weekdy", "satdy", "sundy")) {
    f <- gra_files(tot_co, dt)
    if (length(f)) g[["co"]][[dt]] <- tryCatch(gra_read(f, "CO", win), error = function(e) NULL)
  }
  g[["co"]] <- g[["co"]][!vapply(g[["co"]], is.null, NA)]
  if (!length(g[["co"]])) stop("no CO found for GRA2PES ", ver, " under ", tot_co)
  g$box_t_hr <- component_box_t_hr(g, masks); g$region_t_hr <- component_region_t_hr(g, masks); g$masks <- masks
  g
}

## ---- EPA GHGI -----------------------------------------------------------------------------------
build_epa <- function() {
  nc <- nc_open(GHGI_FILE); on.exit(nc_close(nc))
  glat <- ncvar_get(nc, "lat"); glon <- ncvar_get(nc, "lon"); area <- ncvar_get(nc, "grid_cell_area")   # cm2
  ilat <- which(glat >= DOMAIN["ymn"] & glat <= DOMAIN["ymx"]); ilon <- which(glon >= DOMAIN["xmn"] & glon <= DOMAIN["xmx"])
  rd <- function(v) { a <- ncvar_get(nc, v); if (length(dim(a)) == 3) a <- a[, , 1]; a <- a[ilon, ilat]; a[!is.finite(a)] <- 0; a * GHGI_TO_UMOL }
  sumv <- function(vs) Reduce(`+`, lapply(vs, rd))
  allv <- grep("^emi_ch4_", names(nc$var), value = TRUE)
  og_prod <- c("emi_ch4_1B2a_Petroleum_Systems_Exploration", "emi_ch4_1B2a_Petroleum_Systems_Production",
               "emi_ch4_1B2b_Natural_Gas_Exploration", "emi_ch4_1B2b_Natural_Gas_Production",
               "emi_ch4_1B2b_Natural_Gas_Processing", "emi_ch4_1B2b_Natural_Gas_TransmissionStorage")
  waste <- c("emi_ch4_5A1_Landfills_MSW", "emi_ch4_5A1_Landfills_Industrial", "emi_ch4_5B1_Composting",
             "emi_ch4_5D_Wastewater_Treatment_Domestic", "emi_ch4_5D_Wastewater_Treatment_Industrial")
  ag <- c("emi_ch4_3A_Enteric_Fermentation", "emi_ch4_3B_Manure_Management")
  # EPA separates distribution from production, so OG here = production..transmission and the urban
  # gas term = distribution; they are then split by the same urban mask as GRA2PES for comparability.
  sec <- list(WASTE = sumv(waste), OG = sumv(og_prod) + rd("emi_ch4_1B2b_Natural_Gas_Distribution"),
              RES = rd("emi_ch4_Supp_1B2b_PostMeter"), AG = sumv(ag), TOTAL = sumv(allv))
  sec <- lapply(sec, function(a) array(a, c(dim(a), 1)))
  g <- list(type = "latlon", x = glon[ilon], y = glat[ilat], lat = matrix(glat[ilat], length(ilon), length(ilat), byrow = TRUE),
            lon = matrix(glon[ilon], length(ilon), length(ilat)), hours = NA,
            cell_area_m2 = area[ilon, ilat] * 1e-4, meta = list(name = "epa_ghgi_2020", units = "umol m-2 s-1"))
  masks <- grid_masks(g); g$comps <- make_components(sec, masks); g[["co"]] <- NULL
  g$box_t_hr <- component_box_t_hr(g, masks); g$region_t_hr <- component_region_t_hr(g, masks); g$masks <- masks
  g
}

## ---- run ------------------------------------------------------------------------------------------
tot <- list(); regs <- list()
for (w in which_priors) {
  g <- if (w == "epa") build_epa() else build_gra2pes(w)
  saveRDS(g, file.path(PRIOR_DIR, paste0("prior_", g$meta$name, ".rds")), compress = "xz")
  tot[[g$meta$name]] <- g$box_t_hr; regs[[g$meta$name]] <- g$region_t_hr
  cat(sprintf("\n%s  (t CH4/h inside the Paper 1 box, daily mean)\n", g$meta$name))
  print(round(g$box_t_hr, 3)); cat(sprintf("  total %.2f t/h\n", sum(g$box_t_hr)))
}
B <- do.call(rbind, lapply(names(tot), function(n) data.frame(prior = n, component = names(tot[[n]]), box_t_hr = round(tot[[n]], 4))))
write.csv(B, file.path(PRIOR_DIR, "prior_box_totals.csv"), row.names = FALSE)
RT <- do.call(rbind, lapply(names(tot), function(n) {
  m <- regs[[n]]
  data.frame(prior = n, component = rownames(m), round(m, 4), row.names = NULL) }))
write.csv(RT, file.path(PRIOR_DIR, "prior_region_totals.csv"), row.names = FALSE)
cat("\nprior totals by reporting region (t CH4/h):\n")
print(aggregate(cbind(paper1_box, obs_box, djb) ~ prior, data = RT, FUN = sum), row.names = FALSE)
cat("fossil share by region:\n")
print(do.call(rbind, lapply(split(RT, RT$prior), function(s) { f <- s$component %in% c("og_basin", "og_urban", "postmeter")
  data.frame(prior = s$prior[1], paper1_box = round(sum(s$paper1_box[f]) / sum(s$paper1_box), 2),
             obs_box = round(sum(s$obs_box[f]) / sum(s$obs_box), 2), djb = round(sum(s$djb[f]) / sum(s$djb), 2)) })), row.names = FALSE)
cat("\nwrote prior_*.rds and prior_box_totals.csv to", PRIOR_DIR, "\n")
