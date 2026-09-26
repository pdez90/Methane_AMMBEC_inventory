# grids.R — emission grids on their native projections, and mapping particles onto them --------
# A "grid" is list(type, x, y, lcc (for type "lcc"), lon, lat [nx x ny], comps = list(name =
# array [nx, ny, nh] in umol m-2 s-1), co = list(daytype = array [nx, ny, 24]) or NULL, hours
# (length nh: UTC hour of each slice, or NA for a time-invariant field), meta).
# ----------------------------------------------------------------------------------------------

#' Spherical Lambert conformal conic, forward (Snyder 1987, eqs. 15-1 to 15-4). Returns metres.
lcc_forward <- function(lat, lon, p) {
  d <- pi / 180; f1 <- p$std1 * d; f2 <- p$std2 * d; f0 <- p$lat0 * d
  n <- if (abs(f1 - f2) < 1e-10) sin(f1) else log(cos(f1) / cos(f2)) / log(tan(pi / 4 + f2 / 2) / tan(pi / 4 + f1 / 2))
  F <- cos(f1) * tan(pi / 4 + f1 / 2)^n / n
  rho <- p$R * F / tan(pi / 4 + lat * d / 2)^n; rho0 <- p$R * F / tan(pi / 4 + f0 / 2)^n
  th <- n * ((lon - p$lon0 + 180) %% 360 - 180) * d
  list(x = rho * sin(th) + p$fe, y = rho0 - rho * cos(th) + p$fn)
}

#' Map (lon, lat) to grid cell indices; ok = inside the grid.
grid_index <- function(g, lon, lat) {
  if (g$type == "lcc") {
    xy <- lcc_forward(lat, lon, g$lcc); dx <- g$x[2] - g$x[1]; dy <- g$y[2] - g$y[1]
    ix <- round((xy$x - g$x[1]) / dx) + 1L; iy <- round((xy$y - g$y[1]) / dy) + 1L
  } else {
    dx <- g$x[2] - g$x[1]; dy <- g$y[2] - g$y[1]
    ix <- round((lon - g$x[1]) / dx) + 1L; iy <- round((lat - g$y[1]) / dy) + 1L
  }
  ok <- is.finite(ix) & is.finite(iy) & ix >= 1 & ix <= length(g$x) & iy >= 1 & iy <= length(g$y)
  list(ix = ix, iy = iy, ok = ok)
}

hav_km <- function(lat1, lon1, lat2, lon2) {
  p <- pi / 180; a <- sin((lat2 - lat1) * p / 2)^2 + cos(lat1 * p) * cos(lat2 * p) * sin((lon2 - lon1) * p / 2)^2
  2 * 6371.0088 * asin(pmin(1, sqrt(a)))
}

#' Masks shared by every prior: urban (Paper 1 box grown by URBAN_BUFFER_KM) and one per facility.
grid_masks <- function(g) {
  lat <- g$lat; lon <- g$lon
  b <- URBAN_BOX; buf_lat <- URBAN_BUFFER_KM / 111.2; buf_lon <- URBAN_BUFFER_KM / (111.2 * cos(39.7 * pi / 180))
  urban <- lat >= b$lat_s - buf_lat & lat <= b$lat_n + buf_lat & lon >= b$lon_w - buf_lon & lon <= b$lon_e + buf_lon
  box <- lat >= b$lat_s & lat <= b$lat_n & lon >= b$lon_w & lon <= b$lon_e
  # AREA-WEIGHTED region masks for lat/lon grids (EPA GHGI, 0.1 deg): the box is 6.5 x 4.5 cells, so
  # selecting cells by their centres sums a footprint shifted relative to the box (Paper 1 script 14:
  # 2.62 t/h by centres vs 2.80 t/h area-weighted). Each cell is weighted by the fraction of its area
  # inside the region, assuming uniform emissions within a cell -- identical to Paper 1. For the 1 km
  # LCC grids (GRA2PES) the edge cells are negligible and centre selection (0/1 weights) is kept.
  frac_box <- function(lo_lat, hi_lat, lo_lon, hi_lon) {
    if (g$type == "lcc") return(lat >= lo_lat & lat <= hi_lat & lon >= lo_lon & lon <= hi_lon)
    dlon <- abs(g$x[2] - g$x[1]); dlat <- abs(g$y[2] - g$y[1])
    ov <- function(c, lo, hi, d) pmax(0, pmin(c + d / 2, hi) - pmax(c - d / 2, lo)) / d
    ov(lon, lo_lon, hi_lon, dlon) * ov(lat, lo_lat, hi_lat, dlat)
  }
  box_w <- frac_box(b$lat_s, b$lat_n, b$lon_w, b$lon_e)
  fac <- lapply(seq_len(nrow(FACILITIES)), function(i) {     # cells within the radius, and always the nearest cell
    d <- hav_km(lat, lon, FACILITIES$lat[i], FACILITIES$lon[i]); d <= FACILITIES$radius_km[i] | d == min(d) })
  names(fac) <- FACILITIES$name
  # reporting regions (t/h totals, posterior = alpha x prior): the whole observed domain and the DJ Basin
  # part of it (north of BASIN_LAT, as Paper 1's DJB-core definition). og_basin is zero inside the Paper 1
  # box by construction, so without these the basin factor would never appear in any reported total.
  ob <- OBS_BOX
  obs <- frac_box(ob[["ymn"]], ob[["ymx"]], ob[["xmn"]], ob[["xmx"]])
  djb <- frac_box(BASIN_LAT, ob[["ymx"]], ob[["xmn"]], ob[["xmx"]])
  # cell areas (m2) and, for relocated components (config RELOCATE), the single destination cell
  area <- if (g$type == "lcc") matrix((g$x[2] - g$x[1]) * (g$y[2] - g$y[1]), nrow(lat), ncol(lat)) else g$cell_area_m2
  dest <- lapply(RELOCATE, function(p) { d <- hav_km(lat, lon, p[["lat"]], p[["lon"]]); d == min(d) })
  list(urban = urban, box = box_w, fac = fac, area = area, dest = dest,
       regions = list(paper1_box = box_w, obs_box = obs, djb = djb))
}

#' Move a component's mass (per time slice) into one destination cell, conserving mass.
relocate_mass <- function(a, area, dest) {
  out <- array(0, dim(a)); w <- which(dest)[1]; ai <- as.vector(area)
  for (h in seq_len(dim(a)[3])) { m <- sum(a[, , h] * area); out[, , h][w] <- m / ai[w] }
  out
}

#' Split sector fields into the inversion components (see config.R). `sec` is a list with
#' WASTE, OG, RES, AG, TOTAL arrays [nx, ny, nh]. Returns list of component arrays.
make_components <- function(sec, masks) {
  m3 <- function(m, a) array(rep(as.vector(m), dim(a)[3]), dim(a))
  z <- function(a) { a[!is.finite(a)] <- 0; a }
  W <- z(sec$WASTE); O <- z(sec$OG); Rr <- z(sec$RES); A <- z(sec$AG); Tt <- z(sec$TOTAL)
  anyfac <- Reduce(`|`, masks$fac)
  comps <- list()
  for (nm in names(masks$fac)) {
    comps[[nm]] <- W * m3(masks$fac[[nm]], W)
    if (nm %in% names(masks$dest)) comps[[nm]] <- relocate_mass(comps[[nm]], masks$area, masks$dest[[nm]])
  }
  comps$waste    <- W * m3(!anyfac, W)
  comps$og_basin <- O * m3(!masks$urban, O)
  comps$og_urban <- O * m3(masks$urban, O)
  comps$postmeter <- Rr
  comps$ag       <- A
  oth <- Tt - W - O - Rr - A; oth[oth < 0] <- 0; comps$other <- oth   # (pmax would drop dims)
  comps
}

#' Box total (t CH4 per hour, daily mean) of each component inside the Paper 1 box, for reporting
#' posterior emissions as alpha x prior.
component_box_t_hr <- function(g, masks, mask = masks$box) {
  area_m2 <- if (g$type == "lcc") rep((g$x[2] - g$x[1]) * (g$y[2] - g$y[1]), length(g$lon)) else as.vector(g$cell_area_m2)
  sapply(g$comps, function(a) {
    daily <- apply(a, c(1, 2), mean)                     # umol m-2 s-1
    w <- as.numeric(as.vector(mask))                     # 0/1 (LCC) or area fraction inside the region (lat/lon)
    sum(as.vector(daily) * w * area_m2) * 1e-6 * MW[["CH4"]] * 3600 / 1e6
  })
}
#' Component totals (t/h) for every reporting region in masks$regions: matrix [component, region].
component_region_t_hr <- function(g, masks) sapply(masks$regions, function(m) component_box_t_hr(g, masks, m))
