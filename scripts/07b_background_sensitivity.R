# 07b_background_sensitivity.R — how much do the posterior totals depend on the background? ----------
# The OSSE (06) found the background prior to be the main systematic: totals are biased and under-
# covered even when transport is perfect. On the real data, every v2 fit moves the flight backgrounds
# +3 to +5 ppb above their 5th percentile, and at ZISCALE 1.2 the backgrounds and leg offsets absorb
# half the urban enhancement. This script refits the main configuration with the background treated
# differently, one choice at a time, for every prior, and tabulates the Paper 1 box total, the
# regional totals and the fossil share.
#
# Cases (base = model v2: background centre = flight 5th percentile, sd 5 ppb, leg-offset scale 3 ppb):
#   base        as 07
#   q01, q10    background centre at the flight 1st / 10th percentile
#   clean       centre = median of the flight's segments with the lowest 10% prior-predicted signal
#               (the air the inventory says is least influenced: a model-selected background)
#   sd2, sd10   background prior sd 2 / 10 ppb (tight / loose)
#   leg1.5, leg6  leg-offset scale 1.5 / 6 ppb
#   drift       one background per half flight (allows the background to change during the flight)
# Spread of the box total across cases (log sd, and the min-max ratio) is the background
# uncertainty to carry into the paper next to the transport (ZISCALE) and prior spreads.
#
#   Rscript scripts/07b_background_sensitivity.R                 (all cases, config ch4_c2h6)
#   METHANE_BG_CASES=base,drift,clean METHANE_BG_CONFIG=ch4 Rscript scripts/07b_background_sensitivity.R
#   Runs on whatever transport the environment selects (METHANE_ZISCALE etc.), like 07.
#   Out: <RUN_DIR>/background_sensitivity_<model>.csv, background_sensitivity_summary_<model>.csv,
#        figures/background_sensitivity_<model>.png
# ----------------------------------------------------------------------------------------------
proj <- if (file.exists("config.R")) "." else ".."; source(file.path(proj, "config.R"))
source(file.path(proj, "R", "inversion.R")); source(file.path(proj, "R", "identifiability.R")); source(file.path(proj, "R", "diagnostics.R"))

ALL_CASES <- c("base", "q01", "q10", "clean", "sd2", "sd10", "leg1.5", "leg6", "drift")
CASES <- strsplit(Sys.getenv("METHANE_BG_CASES", paste(ALL_CASES, collapse = ",")), ",")[[1]]
if (length(bad <- setdiff(CASES, ALL_CASES))) stop("unknown case(s): ", paste(bad, collapse = ", "))
CONFIG <- Sys.getenv("METHANE_BG_CONFIG", "ch4_c2h6")
USE <- list(ch4 = character(0), ch4_c2h6 = "c2h6", ch4_c2h6_co = c("c2h6", "co"))[[CONFIG]]
if (is.null(USE) && CONFIG != "ch4") stop("METHANE_BG_CONFIG must be ch4, ch4_c2h6 or ch4_c2h6_co")
cat(sprintf("background sensitivity: model %s, transport %s, config %s, cases %s\n",
            INV_MODEL, TRANSPORT, CONFIG, paste(CASES, collapse = " ")))

priors <- commandArgs(TRUE)
if (!length(priors)) priors <- sub("^jacobian_(.*)\\.rds$", "\\1", list.files(RUN_DIR, pattern = "^jacobian_.*\\.rds$"))
if (!length(priors)) stop("no jacobian_*.rds in ", RUN_DIR, " (run 05 for this transport first)")
S0 <- read.csv(file.path(INV_OUT, "obs_segments.csv"), stringsAsFactors = FALSE)
mod <- inv_model(); out <- list()

for (pr in priors) {
  J <- readRDS(file.path(RUN_DIR, paste0("jacobian_", pr, ".rds")))
  S <- S0[match(J$ids, S0$id), ]; S$leg_key <- paste(S$flight, S$leg_id)
  seen <- colSums(J$H) > 1e-6 * sum(J$H)
  J$H <- J$H[, seen, drop = FALSE]; K <- J$components <- J$components[seen]; J$box_t_hr <- J$box_t_hr[K]
  if (CONFIG == "ch4_c2h6_co" && !isTRUE(J$has_co)) { message("skip ", pr, ": no CO"); next }
  for (cs in CASES) {
    st <- inv_settings(K); Sx <- S; bq <- 0.05; bmu <- NULL
    if (cs == "q01") bq <- 0.01
    if (cs == "q10") bq <- 0.10
    if (cs == "sd2")  st$b_ch4_sd <- 2
    if (cs == "sd10") st$b_ch4_sd <- 10
    if (cs == "leg1.5") st$leg_ch4_scale <- 1.5
    if (cs == "leg6")   st$leg_ch4_scale <- 6
    if (cs == "drift") {                        # background per half flight, split at the median time
      t <- as.POSIXct(Sx$time_utc, tz = "UTC")
      half <- ave(as.numeric(t), Sx$flight, FUN = function(v) as.integer(v > median(v)))
      Sx$flight <- paste0(Sx$flight, "_h", half)
    }
    if (cs == "clean") {                        # model-selected clean air per flight
      sig <- rowSums(J$H); fl <- factor(Sx$flight)
      bmu <- tapply(seq_len(nrow(Sx)), fl, function(i) {
        k <- i[sig[i] <= quantile(sig[i], 0.10)]; median(Sx$ch4[k], na.rm = TRUE) })
    }
    sd <- make_stan_data(Sx, J, use = USE, settings = st, bg_q = bq, b_ch4_mu = bmu)
    t0 <- Sys.time(); fit <- fit_inversion(mod, sd)
    dr <- fit$draws(c("E_total", "fossil_share"), format = "draws_matrix")
    dc <- fit_decomposition(fit, Sx, sd); ub <- dc$by_region[dc$by_region$region == "urban", ]
    rg <- region_posterior(fit, K, J$region_t_hr)
    rgv <- function(r, col) if (!is.null(rg) && r %in% rg$region) rg[[col]][rg$region == r] else NA_real_
    dg <- fit$diagnostic_summary(quiet = TRUE)
    row <- data.frame(prior = pr, case = cs, config = CONFIG, transport = TRANSPORT,
      bg_centre_median = round(median(sd$b_ch4_mu), 1), n_backgrounds = sd$F,
      E_prior = round(sum(J$box_t_hr), 2),
      E_q05 = round(quantile(dr[, "E_total"], 0.05), 2), E_q50 = round(median(dr[, "E_total"]), 2),
      E_q95 = round(quantile(dr[, "E_total"], 0.95), 2),
      fossil_q50 = round(median(dr[, "fossil_share"]), 3),
      obs_box_q50 = rgv("obs_box", "E_q50"), djb_q50 = rgv("djb", "E_q50"),
      djb_fossil_share = rgv("djb", "fossil_share_q50"),
      bg_shift_median = round(median(dc$bg_shift), 2), urban_signal_share = round(ub$signal_share, 2),
      urban_abs_leg = round(ub$abs_leg, 2), leg_scale = round(dc$nuisance[["sl_ch4"]], 2),
      phi = round(dc$nuisance[["phi"]], 2), divergent = sum(dg$num_divergent),
      max_rhat = round(max(fit$summary("alpha")$rhat), 3), row.names = NULL)
    out[[length(out) + 1]] <- row
    with(row, cat(sprintf("%-17s %-7s box %.1f t/h [%.1f-%.1f]  fossil %.2f  DJB %.1f  bg shift %+.1f  urban signal %.0f%%  |leg| %.1f  div %d rhat %.3f  (%.1f min)\n",
      prior, case, E_q50, E_q05, E_q95, fossil_q50, djb_q50, bg_shift_median, 100 * urban_signal_share, urban_abs_leg,
      divergent, max_rhat, as.numeric(difftime(Sys.time(), t0, units = "mins")))))
  }
}
B <- do.call(rbind, out)
write.csv(B, inv_file("background_sensitivity.csv"), row.names = FALSE)

## summary: spread across background cases, per prior, compared with the spread across priors -----
lsd <- function(v) sd(log(pmax(v, 1e-9)))
SUM <- do.call(rbind, lapply(split(B, B$prior), function(s) data.frame(prior = s$prior[1],
  n_cases = nrow(s), box_min = min(s$E_q50), box_base = s$E_q50[s$case == "base"][1], box_max = max(s$E_q50),
  box_max_over_min = round(max(s$E_q50) / min(s$E_q50), 2), box_log_sd = round(lsd(s$E_q50), 3),
  fossil_min = min(s$fossil_q50), fossil_max = max(s$fossil_q50),
  djb_min = min(s$djb_q50, na.rm = TRUE), djb_max = max(s$djb_q50, na.rm = TRUE))))
PR <- do.call(rbind, lapply(split(B, B$case), function(s) data.frame(case = s$case[1],
  box_prior_spread_log_sd = round(lsd(s$E_q50), 3),
  prior_dependence = round(lsd(s$E_q50) / lsd(s$E_prior), 2))))
write.csv(SUM, inv_file("background_sensitivity_summary.csv"), row.names = FALSE)
cat("\nPaper 1 box total (posterior median, t/h) by background case and prior:\n")
W <- reshape(B[, c("prior", "case", "E_q50")], idvar = "case", timevar = "prior", direction = "wide")
names(W) <- sub("^E_q50\\.", "", names(W)); print(W, row.names = FALSE)
cat("\nspread across background cases, per prior:\n"); print(SUM, row.names = FALSE)
cat("\nspread across priors, per background case (prior dependence 0 = data-determined):\n"); print(PR, row.names = FALSE)
cat(sprintf("\nbackground uncertainty in the box total: median over priors of the across-case log sd = %.2f (~x%.2f)\n",
            median(SUM$box_log_sd), exp(median(SUM$box_log_sd))))

png(file.path(RUN_DIR, "figures", paste0("background_sensitivity_", INV_MODEL, ".png")), 1500, 800, res = 150)
par(mar = c(6, 4.5, 3, 1)); cols <- c("#1b9e77", "#d95f02", "#7570b3", "#e7298a")
pp <- unique(B$prior); cc <- intersect(ALL_CASES, unique(B$case)); x <- seq_along(cc)
yl <- range(c(B$E_q05, B$E_q95), na.rm = TRUE)
plot(NA, xlim = c(0.5, length(cc) + 0.5), ylim = yl, xaxt = "n", xlab = "", ylab = "Paper 1 box total (t CH4/h)",
     main = sprintf("Background sensitivity (%s, %s, %s)", CONFIG, TRANSPORT, INV_MODEL))
axis(1, x, cc, las = 2)
for (j in seq_along(pp)) { s <- B[B$prior == pp[j], ]; i <- match(s$case, cc); off <- (j - (length(pp) + 1) / 2) * 0.18
  segments(i + off, s$E_q05, i + off, s$E_q95, col = cols[j], lwd = 2); points(i + off, s$E_q50, pch = 19, col = cols[j]) }
legend("topright", bty = "n", pch = 19, col = cols[seq_along(pp)], legend = pp, cex = 0.8)
dev.off()
cat("wrote background_sensitivity_*.csv and figures/background_sensitivity_", INV_MODEL, ".png to ", RUN_DIR, "\n", sep = "")
