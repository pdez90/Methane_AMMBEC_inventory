# 10_figures.R — Figures 1–6 of the Paper 2 manuscript, from the run outputs ------------------------------
#   Rscript scripts/10_figures.R            (after run_paper.R; env METHANE_MAIN_ZI=0.8 METHANE_RANGE_ZI=0.8,1.0)
#   Out: <INV_OUT>/figures/paper2/fig{1..6}_*.png (300 dpi) and .pdf
# Base graphics only. Every panel reads a CSV/RDS that 04–09 wrote; nothing is recomputed here.
#   1  Map: receptors by region, metro box, observed domain, DJ Basin line, facilities; GRA2PES v2.0b prior
#   2  Prior vs posterior scaling factors by component and prior (CH4+C2H6), coloured by regime
#   3  Transport sweep: totals and fossil shares vs ZISCALE, with the aircraft-profile constraint; aircraft vs HRRR MLH
#   4  Identifiability: information gain and fossil-share sd by tracer set and prior; 06b regime map
#   5  Background sensitivity: box total and fossil share by background case and prior
#   6  Decomposition of the observed urban enhancement by transport setting (signal / background / leg / residual)
proj <- if (file.exists("config.R")) "." else ".."; source(file.path(proj, "config.R")); source(file.path(proj, "R", "grids.R"))
MAIN_ZI  <- as.numeric(Sys.getenv("METHANE_MAIN_ZI", "0.8"))
RANGE_ZI <- as.numeric(strsplit(Sys.getenv("METHANE_RANGE_ZI", "0.8,1.0"), ",")[[1]])
MODEL <- Sys.getenv("METHANE_INV_MODEL", "v2"); CFG_MAIN <- "ch4_c2h6"
FIG <- file.path(INV_OUT, "figures", "paper2"); dir.create(FIG, FALSE, TRUE)
runs_root <- file.path(INV_OUT, "runs")
zi_of <- function(d) { m <- regmatches(d, regexpr("_zi([0-9.]+)$", d)); if (length(m)) as.numeric(sub("_zi", "", m)) else 1 }
rd <- function(...) { f <- file.path(...); if (file.exists(f)) read.csv(f, stringsAsFactors = FALSE) else NULL }
dirs <- list.dirs(runs_root, recursive = FALSE, full.names = FALSE); dirs <- dirs[grepl("^np1000_h10", dirs)]
main_dir <- dirs[sapply(dirs, zi_of) == MAIN_ZI][1]; base_dir <- "np1000_h10"   # 05c/06/06b write to the unscaled dir
PRIORS <- c("epa_ghgi_2020", "gra2pes_v1.1", "gra2pes_v2.0beta")
PRIOR_LAB <- c(epa_ghgi_2020 = "EPA GHGI 2020", gra2pes_v1.1 = "GRA2PES v1.1", gra2pes_v2.0beta = "GRA2PES v2.0b")
PCOL <- c(epa_ghgi_2020 = "#009E73", gra2pes_v1.1 = "#CC79A7", gra2pes_v2.0beta = "#0072B2")
RCOL <- c(urban = "#D55E00", edge = "#E69F00", basin = "#0072B2")
CFG_LAB <- c(ch4 = "CH4", ch4_c2h6 = "CH4 + C2H6", ch4_c2h6_co = "CH4 + C2H6 + CO")
COMP_LAB <- c(dads = "DADS landfill", tower_road = "Tower Road landfill", metro_complex = "Metro / Suncor", waste = "Other waste",
              og_basin = "O&G, DJ Basin", og_urban = "O&G, urban", postmeter = "Post-meter gas", ag = "Agriculture", other = "Other")
REGCOL <- c("well identified" = "#009E73", "partially identified" = "#E69F00", "not identified" = "#999999")
dev_open <- function(name, w, h) {   # PNG at 300 dpi and PDF, same content
  list(png = function() png(file.path(FIG, paste0(name, ".png")), w, h, units = "in", res = 300, type = if (capabilities("cairo")) "cairo" else "quartz"),
       pdf = function() pdf(file.path(FIG, paste0(name, ".pdf")), w, h))
}
draw_both <- function(name, w, h, fun) { for (o in dev_open(name, w, h)) { o(); fun(); dev.off() }; cat("wrote", name, "\n") }
panel_label <- function(l) { u <- par("usr"); if (par("xlog")) u[1:2] <- 10^u[1:2]; mtext(l, side = 3, line = 0.4, at = u[1] - 0.04 * (u[2] - u[1]), adj = 1, font = 2, cex = 1.05) }

## ---- Figure 1: map -------------------------------------------------------------------------------------
obs <- rd(INV_OUT, "obs_segments.csv")
g2 <- { f <- file.path(PRIOR_DIR, "prior_gra2pes_v2.0beta.rds"); if (file.exists(f)) readRDS(f) else NULL }
fig1 <- function() {
  layout(matrix(1:2, 1), widths = c(1, 1.15)); par(mar = c(4, 4, 2, 1), mgp = c(2.2, 0.7, 0))
  b <- URBAN_BOX; ob <- OBS_BOX
  plot(NA, xlim = c(ob[["xmn"]] - 0.05, ob[["xmx"]] + 0.05), ylim = c(ob[["ymn"]] - 0.05, ob[["ymx"]] + 0.05),
       xlab = expression("Longitude (" * degree * "E)"), ylab = expression("Latitude (" * degree * "N)"), asp = 1 / cos(40.2 * pi / 180), las = 1)
  rect(ob[["xmn"]], ob[["ymn"]], ob[["xmx"]], ob[["ymx"]], border = "grey40", lty = 2)
  segments(ob[["xmn"]], BASIN_LAT, ob[["xmx"]], BASIN_LAT, col = "grey40", lty = 3)
  if (!is.null(obs)) points(obs$lon, obs$lat, pch = 16, cex = 0.35, col = adjustcolor(RCOL[obs$region], 0.7))
  rect(b$lon_w, b$lat_s, b$lon_e, b$lat_n, border = "black", lwd = 1.6)
  points(FACILITIES$lon, FACILITIES$lat, pch = 24, bg = "white", cex = 1.1)
  text(FACILITIES$lon, FACILITIES$lat, c("DADS", "Tower Rd", "Metro"), pos = c(1, 3, 2), cex = 0.65, offset = 0.5)
  legend("topright", c("urban", "urban edge", "basin", "Denver metro box", "observed domain", "DJ Basin boundary (40.05 N)"), pch = c(16, 16, 16, NA, NA, NA), lty = c(NA, NA, NA, 1, 2, 3),
         col = c(RCOL, "black", "grey40", "grey40"), pt.cex = 0.9, cex = 0.65, bg = "white", box.col = "grey70")
  panel_label("a")
  if (!is.null(g2)) {
    tot <- Reduce(`+`, lapply(g2$comps, function(a) apply(a, c(1, 2), mean)))          # umol m-2 s-1, daily mean
    xk <- g2$x / 1000; yk <- g2$y / 1000; z <- log10(pmax(tot, 1e-4))
    cols <- hcl.colors(64, "YlOrRd", rev = TRUE)
    oc4 <- lcc_forward(c(ob[["ymn"]], ob[["ymn"]], ob[["ymx"]], ob[["ymx"]]), c(ob[["xmn"]], ob[["xmx"]], ob[["xmn"]], ob[["xmx"]]), g2$lcc)
    image(xk, yk, z, col = cols, xlab = "x (km, LCC)", ylab = "y (km, LCC)", las = 1, asp = 1, zlim = c(-2.5, max(z, na.rm = TRUE)),
          xlim = range(oc4$x / 1000) + c(-40, 40), ylim = range(oc4$y / 1000) + c(-40, 40))
    rb <- function(lat_s, lat_n, lon_w, lon_e, ...) { cc <- lcc_forward(c(lat_s, lat_s, lat_n, lat_n, lat_s), c(lon_w, lon_e, lon_e, lon_w, lon_w), g2$lcc)
      lines(cc$x / 1000, cc$y / 1000, ...) }
    rb(b$lat_s, b$lat_n, b$lon_w, b$lon_e, lwd = 1.6); rb(ob[["ymn"]], ob[["ymx"]], ob[["xmn"]], ob[["xmx"]], lty = 2, col = "grey30")
    bl <- lcc_forward(c(BASIN_LAT, BASIN_LAT), c(ob[["xmn"]], ob[["xmx"]]), g2$lcc); lines(bl$x / 1000, bl$y / 1000, lty = 3, col = "grey30")
    fc <- lcc_forward(FACILITIES$lat, FACILITIES$lon, g2$lcc); points(fc$x / 1000, fc$y / 1000, pch = 24, bg = "white", cex = 1)
    if (!is.null(obs)) { oc <- lcc_forward(obs$lat, obs$lon, g2$lcc); points(oc$x / 1000, oc$y / 1000, pch = ".", col = adjustcolor("black", 0.5)) }
    # colour key
    usr <- par("usr"); zz <- seq(-2.5, max(z, na.rm = TRUE), length.out = 64)
    xl <- usr[2] - 0.06 * diff(usr[1:2]); yb <- usr[3] + 0.08 * diff(usr[3:4]); yt <- yb + 0.35 * diff(usr[3:4])
    rect(xl, yb + (seq_along(zz) - 1) / 64 * (yt - yb), xl + 0.02 * diff(usr[1:2]), yb + seq_along(zz) / 64 * (yt - yb), col = cols, border = NA)
    at <- c(-2, -1, 0, 1); ok <- at >= min(zz) & at <= max(zz)
    text(xl - 0.005 * diff(usr[1:2]), yb + (at[ok] - min(zz)) / diff(range(zz)) * (yt - yb), c("0.01", "0.1", "1", "10")[ok], adj = 1, cex = 0.6)
    text(xl, yt + 0.03 * diff(usr[3:4]), expression(CH[4] ~ (mu * mol ~ m^-2 ~ s^-1)), adj = 1, cex = 0.6)
    panel_label("b"); mtext("GRA2PES v2.0b prior, daily mean", side = 3, line = 0.4, cex = 0.8)
  }
}
draw_both("fig1_map", 9.5, 4.6, fig1)

## ---- Figure 2: scaling factors -------------------------------------------------------------------------
A <- rd(runs_root, main_dir, paste0("posterior_alpha_", MODEL, ".csv"))
IC <- rd(runs_root, main_dir, paste0("posterior_identifiability_components_", MODEL, ".csv"))
fig2 <- function() {
  par(mfrow = c(1, 3), mar = c(4, 9, 3, 1), mgp = c(2.2, 0.7, 0))
  comps <- names(COMP_LAB)
  for (p in PRIORS) {
    a <- A[A$prior == p & A$config == CFG_MAIN, ]; a <- a[match(comps, a$component), ]
    r <- IC[IC$prior == p & IC$config == CFG_MAIN, ]; reg <- r$regime[match(comps, r$component)]
    y <- rev(seq_along(comps)); ok <- is.finite(a$q50)
    plot(NA, xlim = c(0.03, 300), ylim = c(0.5, length(comps) + 0.5), log = "x", xaxt = "n", yaxt = "n", ylab = "", las = 1,
         xlab = expression(alpha ~ "(posterior / prior)"), main = PRIOR_LAB[p], cex.main = 0.95)
    axis(1, at = c(0.1, 1, 10, 100), labels = c("0.1", "1", "10", "100"))
    abline(v = 1, col = "grey50"); abline(v = c(0.1, 10, 100), col = "grey90")
    axis(2, at = y, labels = COMP_LAB, las = 1, cex.axis = 0.8)
    segments(a$q05[ok], y[ok], a$q95[ok], y[ok], col = REGCOL[reg[ok]], lwd = 2)
    points(a$q50[ok], y[ok], pch = 21, bg = REGCOL[reg[ok]], cex = 1.2)
    if (any(!ok)) text(0.03, y[!ok], "absent from prior", adj = 0, cex = 0.65, col = "grey50")
    big <- ok & a$q50 > 100; if (any(big)) text(a$q50[big], y[big] + 0.35, sprintf("alpha = %.0f", a$q50[big]), cex = 0.65)
    if (p == PRIORS[1]) legend("topright", names(REGCOL), pch = 21, pt.bg = REGCOL, cex = 0.7, bg = "white", box.col = "grey70")
    panel_label(letters[match(p, PRIORS)])
  }
}
draw_both("fig2_scaling_factors", 11, 4.2, fig2)

## ---- Figure 3: transport sweep + aircraft MLH ----------------------------------------------------------
RT <- do.call(rbind, lapply(dirs, function(d) { r <- rd(runs_root, d, paste0("posterior_region_totals_", MODEL, ".csv")); if (is.null(r)) return(NULL); r$ziscale <- zi_of(d); r }))
RT <- RT[RT$config == CFG_MAIN, ]
BLH <- rd(runs_root, base_dir, "blh_aircraft_profiles.csv")
blh_txt <- file.path(runs_root, base_dir, "blh_aircraft_summary.txt")
blh_zi <- if (file.exists(blh_txt)) { l <- grep("implied ZISCALE", readLines(blh_txt), value = TRUE); if (length(l)) as.numeric(sub(".*uncensored\\): *([0-9.]+).*", "\\1", l[1])) else NA } else NA
sweep_panel <- function(reg, col, ylab, main, share = FALSE, lab) {
  x <- RT[RT$region == reg, ]; zs <- sort(unique(x$ziscale))
  yr <- if (share) c(0, 1) else c(0, max(x[[if (share) "fossil_share_q95" else "E_q95"]]) * 1.05)
  plot(NA, xlim = c(0.3, 1.25), ylim = yr, xlab = "Mixed-layer height scaling (ZISCALE)", ylab = ylab, las = 1, main = main, cex.main = 0.95)
  rect(min(RANGE_ZI), yr[1] - 1, max(RANGE_ZI), yr[2] * 2, col = adjustcolor("grey80", 0.5), border = NA)
  if (is.finite(blh_zi)) { abline(v = blh_zi, lty = 2, col = "grey30"); text(blh_zi, yr[2] * 0.98, sprintf("aircraft >= %.2f", blh_zi), adj = c(-0.05, 1), cex = 0.65, col = "grey30") }
  for (p in PRIORS) { s <- x[x$prior == p, ]; s <- s[order(s$ziscale), ]
    if (share) { lines(s$ziscale, s$fossil_share_q50, col = PCOL[p], lwd = 2); arrows(s$ziscale, s$fossil_share_q05, s$ziscale, s$fossil_share_q95, angle = 90, code = 3, length = 0.02, col = PCOL[p]) }
    else { lines(s$ziscale, s$E_q50, col = PCOL[p], lwd = 2); arrows(s$ziscale, s$E_q05, s$ziscale, s$E_q95, angle = 90, code = 3, length = 0.02, col = PCOL[p]) }
    points(s$ziscale, if (share) s$fossil_share_q50 else s$E_q50, pch = 21, bg = PCOL[p]) }
  panel_label(lab)
}
fig3 <- function() {
  par(mfrow = c(2, 3), mar = c(4, 4.2, 2.5, 1), mgp = c(2.2, 0.7, 0))
  sweep_panel("paper1_box", PCOL, expression("Emission (t CH"[4] ~ "h"^-1 * ")"), "Metro box (urban)", FALSE, "a")
  legend("topleft", PRIOR_LAB, col = PCOL, lwd = 2, pch = 21, pt.bg = PCOL, cex = 0.7, bg = "white", box.col = "grey70")
  sweep_panel("obs_box", PCOL, expression("Emission (t CH"[4] ~ "h"^-1 * ")"), "Observed domain", FALSE, "b")
  sweep_panel("djb", PCOL, expression("Emission (t CH"[4] ~ "h"^-1 * ")"), "DJ Basin", FALSE, "c")
  sweep_panel("obs_box", PCOL, "Fossil share", "Observed domain", TRUE, "d")
  sweep_panel("djb", PCOL, "Fossil share", "DJ Basin", TRUE, "e")
  # (f) aircraft vs HRRR mixed-layer height
  if (!is.null(BLH)) {
    unc <- BLH[!is.na(BLH$zi_theta) & !is.na(BLH$hrrr_mlht), ]; cen <- BLH[is.na(BLH$zi_theta) & !is.na(BLH$zi_lower_bound) & !is.na(BLH$hrrr_mlht), ]
    lim <- c(0, max(c(unc$zi_theta, unc$hrrr_mlht, cen$zi_lower_bound, cen$hrrr_mlht), na.rm = TRUE) * 1.05)
    plot(unc$hrrr_mlht, unc$zi_theta, xlim = lim, ylim = lim, pch = 21, bg = "#0072B2", xlab = "HRRR mixed-layer height (m AGL)",
         ylab = "Aircraft mixed-layer height (m AGL)", las = 1, main = "Aircraft profiles vs HRRR", cex.main = 0.95)
    points(cen$hrrr_mlht, cen$zi_lower_bound, pch = 24, bg = "white", cex = 0.9)
    arrows(cen$hrrr_mlht, cen$zi_lower_bound, cen$hrrr_mlht, cen$zi_lower_bound + 0.06 * lim[2], length = 0.04, col = "grey40")
    abline(0, 1, col = "grey50"); if (is.finite(blh_zi)) abline(0, blh_zi, lty = 2, col = "grey30")
    if (any(!is.na(unc$lidar_blh))) points(unc$hrrr_mlht, unc$lidar_blh, pch = 4, col = "#D55E00")
    legend("topleft", c(sprintf("resolved top (n = %d)", nrow(unc)), sprintf("lower bound, top above profile (n = %d)", nrow(cen)), "lidar at same time", "1:1", if (is.finite(blh_zi)) sprintf("slope %.2f", blh_zi)),
           pch = c(21, 24, 4, NA, NA), pt.bg = c("#0072B2", "white", NA, NA, NA), col = c("black", "black", "#D55E00", "grey50", "grey30"), lty = c(NA, NA, NA, 1, 2), cex = 0.62, bg = "white", box.col = "grey70")
    panel_label("f")
  }
}
draw_both("fig3_transport_sweep", 11, 7, fig3)

## ---- Figure 4: identifiability --------------------------------------------------------------------------
IO <- rd(runs_root, main_dir, paste0("posterior_identifiability_overall_", MODEL, ".csv"))
RG <- rd(runs_root, base_dir, paste0("regimes_overall_", MODEL, ".csv"))
grouped <- function(val, ylab, main, lab, ylim = NULL) {
  cfgs <- names(CFG_LAB); m <- sapply(PRIORS, function(p) sapply(cfgs, function(cf) { v <- IO[[val]][IO$prior == p & IO$config == cf]; if (length(v)) v[1] else NA }))
  bp <- barplot(m, beside = TRUE, names.arg = PRIOR_LAB[PRIORS], col = c("#BBBBBB", "#0072B2", "#56B4E9"), border = NA, las = 1, ylab = ylab, main = main, cex.main = 0.95,
                ylim = if (is.null(ylim)) c(0, max(m, na.rm = TRUE) * 1.35) else ylim, cex.names = 0.8)
  text(bp, m, ifelse(is.na(m), "", formatC(m, format = "f", digits = if (val == "sd_fossil_share") 2 else 1)), pos = 3, cex = 0.65, xpd = TRUE)
  legend("topleft", CFG_LAB, fill = c("#BBBBBB", "#0072B2", "#56B4E9"), border = NA, cex = 0.7, bg = "white", box.col = "grey70")
  panel_label(lab)
}
fig4 <- function() {
  layout(matrix(c(1, 2, 3, 3), 2, byrow = TRUE), heights = c(1, 0.85)); par(mar = c(4, 4.2, 2.5, 1), mgp = c(2.2, 0.7, 0))
  grouped("info_gain_bits", "Information gain (bits)", "Information gained from the data", "a")
  grouped("sd_fossil_share", "Posterior sd of box fossil share", "Uncertainty of the urban fossil share", "b")
  if (!is.null(RG)) {   # 06b: gain from ethane by endmember contrast x noise (sampling = all), mean over replicates
    r <- RG[RG$config == CFG_MAIN & RG$sampling == "all", ]
    t <- aggregate(cbind(dIG_vs_ch4_bits, sd_fossil_share_ratio_vs_ch4) ~ contrast + noise, data = r, FUN = mean)
    cl <- unique(t$contrast); nl <- unique(t$noise)
    x <- match(t$contrast, cl); y <- match(t$noise, nl)
    plot(NA, xlim = c(0.5, length(cl) + 0.5), ylim = c(0.5, length(nl) + 0.5), xaxt = "n", yaxt = "n", xlab = "Ethane endmember contrast (basin vs urban gas)", ylab = "Ethane noise",
         main = "Synthetic experiment: what ethane adds (06b)", cex.main = 0.95)
    CLAB <- c(none = "none (same ratio)", paper1 = "Paper 1 endmembers", wide = "wide"); NLAB <- c(low = "low", high = "high")
    axis(1, seq_along(cl), ifelse(cl %in% names(CLAB), CLAB[cl], cl), cex.axis = 0.8); axis(2, seq_along(nl), ifelse(nl %in% names(NLAB), NLAB[nl], nl), las = 1, cex.axis = 0.8)
    sz <- 0.8 + 2.2 * (t$dIG_vs_ch4_bits - min(t$dIG_vs_ch4_bits)) / max(1e-9, diff(range(t$dIG_vs_ch4_bits)))
    cc <- hcl.colors(20, "Blues 3", rev = TRUE)[pmax(1, pmin(20, round(20 * (1 - t$sd_fossil_share_ratio_vs_ch4) / 0.6)))]
    points(x, y, pch = 21, cex = sz * 2.2, bg = cc)
    text(x, y, sprintf("+%.1f bits\nsd x%.2f", t$dIG_vs_ch4_bits, t$sd_fossil_share_ratio_vs_ch4), cex = 0.62, pos = 4, offset = 1.4)
    mtext("circle size: information gain from ethane; colour: reduction of fossil-share sd relative to CH4 alone", side = 1, line = 3, cex = 0.6)
    panel_label("c")
  }
}
draw_both("fig4_identifiability", 9, 7.5, fig4)

## ---- Figure 5: background sensitivity --------------------------------------------------------------------
BG <- rd(runs_root, main_dir, paste0("background_sensitivity_", MODEL, ".csv"))
fig5 <- function() {
  if (is.null(BG)) return(invisible())
  cases <- c("base", "q01", "q10", "clean", "sd2", "sd10", "leg1.5", "leg6", "drift"); cases <- cases[cases %in% BG$case]
  CL <- c(base = "base (p05, sd 5 ppb, leg 3 ppb)", q01 = "1st percentile", q10 = "10th percentile", clean = "clean-air centre", sd2 = "prior sd 2 ppb", sd10 = "prior sd 10 ppb", leg1.5 = "leg scale 1.5 ppb", leg6 = "leg scale 6 ppb", drift = "within-flight drift")
  par(mfrow = c(1, 2), mar = c(4, 10, 2.5, 1), mgp = c(2.2, 0.7, 0))
  y0 <- rev(seq_along(cases)); off <- c(-0.22, 0, 0.22); names(off) <- PRIORS
  plot(NA, xlim = c(0, max(BG$E_q95) * 1.05), ylim = c(0.5, length(cases) + 0.5), yaxt = "n", ylab = "", xlab = expression("Metro box total (t CH"[4] ~ "h"^-1 * ")"), las = 1, main = "Urban total", cex.main = 0.95)
  axis(2, y0, CL[cases], las = 1, cex.axis = 0.75); abline(h = y0 + 0.5, col = "grey92")
  for (p in PRIORS) { s <- BG[BG$prior == p, ]; s <- s[match(cases, s$case), ]; yy <- y0 + off[p]
    segments(s$E_q05, yy, s$E_q95, yy, col = PCOL[p]); points(s$E_q50, yy, pch = 21, bg = PCOL[p], cex = 0.9) }
  legend("topright", PRIOR_LAB, pch = 21, pt.bg = PCOL, col = PCOL, cex = 0.7, bg = "white", box.col = "grey70"); panel_label("a")
  plot(NA, xlim = c(0, 1), ylim = c(0.5, length(cases) + 0.5), yaxt = "n", ylab = "", xlab = "Metro box fossil share", las = 1, main = "Urban fossil share", cex.main = 0.95)
  axis(2, y0, CL[cases], las = 1, cex.axis = 0.75); abline(h = y0 + 0.5, col = "grey92")
  for (p in PRIORS) { s <- BG[BG$prior == p, ]; s <- s[match(cases, s$case), ]; points(s$fossil_q50, y0 + off[p], pch = 21, bg = PCOL[p], cex = 0.9) }
  panel_label("b")
}
draw_both("fig5_background", 10, 4.5, fig5)

## ---- Figure 6: decomposition of the urban enhancement by transport ----------------------------------------
DC <- do.call(rbind, lapply(dirs, function(d) { x <- rd(runs_root, d, paste0("posterior_decomposition_", MODEL, ".csv")); if (is.null(x)) return(NULL); x$ziscale <- zi_of(d); x }))
fig6 <- function() {
  if (is.null(DC)) return(invisible())
  d <- DC[DC$config == CFG_MAIN & DC$region == "urban", ]; zs <- sort(unique(d$ziscale))
  par(mfrow = c(1, 3), mar = c(4, 4.2, 2.5, 1), mgp = c(2.2, 0.7, 0))
  SC <- c(fitted_signal = "#0072B2", bg_shift = "#E69F00", abs_leg = "#CC79A7")
  for (p in PRIORS) {
    s <- d[d$prior == p, ]; s <- s[match(zs, s$ziscale), ]
    m <- rbind(s$fitted_signal, s$bg_shift, s$abs_leg); colnames(m) <- zs
    obs <- mean(s$obs_above_p05, na.rm = TRUE)
    bp <- barplot(m, beside = TRUE, col = SC, border = NA, las = 1, ylim = c(0, max(m, obs, na.rm = TRUE) * 1.25),
                  xlab = "ZISCALE", ylab = "Urban CH4 above flight 5th percentile (ppb, median over segments)", main = PRIOR_LAB[p], cex.main = 0.95)
    abline(h = obs, lty = 2, col = "black"); text(bp[1, 1], obs, sprintf("observed %.1f ppb", obs), adj = c(0, -0.4), cex = 0.7)
    text(colMeans(bp), apply(m, 2, max), sprintf("%.0f%%", 100 * s$signal_share), pos = 3, cex = 0.7)
    xr <- range(bp[, zs %in% RANGE_ZI]); rect(xr[1] - 0.7, -1, xr[2] + 0.7, par("usr")[4], border = "grey40", lty = 3)
    if (p == PRIORS[1]) legend("topleft", c("fitted emission signal", "background shift above p05", "|leg offset|", "observed enhancement", "% = sum(signal) / sum(observed) over urban segments"),
                               fill = c(SC, NA, NA), border = NA, lty = c(NA, NA, NA, 2, NA), cex = 0.62, bg = "white", box.col = "grey70")
    panel_label(letters[match(p, PRIORS)])
  }
}
draw_both("fig6_decomposition", 11, 4.2, fig6)
cat("figures written to ", FIG, "\n", sep = "")
