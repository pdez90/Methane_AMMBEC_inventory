# 06_osse.R — synthetic-data tests: what can these flights identify? -----------------------------
# Observing-system simulation experiments with the real flight tracks and footprints. A "true"
# set of scale factors, ethane ratios, transport factors, backgrounds and errors is drawn, synthetic
# CH4 / C2H6 / CO segment values are generated with the forward model of stan/inversion.stan, and
# the model is fitted three ways:
#
#   ch4           CH4 only
#   ch4_c2h6      CH4 + ethane
#   ch4_c2h6_co   CH4 + ethane + CO
#
# For each component: does the 90% interval cover the truth, and how much narrower is the
# posterior than the prior (error reduction = 1 - posterior sd / prior sd of log alpha)? Also the
# posterior correlation between basin and urban oil and gas (can the data tell them apart?) and
# the recovered fossil share in the Paper 1 box. Replicate 1 uses a Paper-1-like truth relative to
# GRA2PES v2.0beta (landfills far too high, the Metro complex too low); the others are random.
#
#   Rscript scripts/06_osse.R [prior=gra2pes_v2.0beta] [n_rep=4]
#   Out: <RUN_DIR>/osse_results.csv, osse_summary.csv, osse_identifiability_components.csv,
#        osse_identifiability_overall.csv (DOFS, information gain, fossil/biogenic correlation, and
#        the gain from each tracer vs CH4 alone), osse_posterior_correlations.csv, figures/osse_recovery.png
# ----------------------------------------------------------------------------------------------
proj <- if (file.exists("config.R")) "." else ".."; source(file.path(proj, "config.R"))
source(file.path(proj, "R", "inversion.R")); source(file.path(proj, "R", "identifiability.R"))
if (nzchar(Sys.getenv("CMDSTAN"))) set_cmdstan_path(Sys.getenv("CMDSTAN"))

a <- commandArgs(TRUE)
PRIOR <- if (length(a) >= 1) a[1] else "gra2pes_v2.0beta"
N_REP <- if (length(a) >= 2) as.integer(a[2]) else 4L
set.seed(20240713)
REALISTIC <- Sys.getenv("METHANE_OSSE_REALISTIC", "1") == "1"   # heavy-tailed errors, larger leg offsets
cat("OSSE: model", INV_MODEL, if (REALISTIC) "(realistic errors: t3, phi 0.5, leg sd 5)" else "(Gaussian errors)", "\n")

J <- readRDS(file.path(RUN_DIR, paste0("jacobian_", PRIOR, ".rds")))
S <- read.csv(file.path(INV_OUT, "obs_segments.csv"), stringsAsFactors = FALSE)
S <- S[match(J$ids, S$id), ]; S$leg_key <- paste(S$flight, S$leg_id)
K <- J$components; st <- inv_settings(K)
# components the flights cannot see at all (no footprint overlap) are dropped from the test
seen <- colSums(J$H) > 1e-6 * sum(J$H)
if (any(!seen)) message("not seen by any footprint, dropped: ", paste(K[!seen], collapse = ", "))
J$H <- J$H[, seen, drop = FALSE]; J$components <- K <- K[seen]; J$box_t_hr <- J$box_t_hr[K]; st <- inv_settings(K)
fl <- factor(S$flight); lg <- factor(S$leg_key)

truth_for <- function(rep) {
  al <- if (rep == 1) {
    t <- c(dads = 0.1, tower_road = 0.15, metro_complex = 3, metro_wwtp = 3, waste = 1, og_basin = 1.5, og_urban = 1.3,
           postmeter = 1, ag = 1, other = 1)
    t[K]
  } else setNames(exp(rnorm(length(K), 0, 0.7)), K)
  al[is.na(al)] <- 1
  rp <- sapply(K, function(k) if (k %in% names(R_PRIOR)) R_PRIOR[[k]] else c(mean = 0, sd = 0))
  r <- pmax(0, rp["mean", ] + rp["sd", ] * rnorm(length(K)))
  list(alpha = al, r = r, tau = exp(rnorm(nlevels(fl), 0, 0.3)), gamma = exp(rnorm(1, 0, 0.1)),
       b = list(ch4 = 1960 + rnorm(nlevels(fl), 0, 10), c2h6 = 1.5 + rnorm(nlevels(fl), 0, 0.3), co = 100 + rnorm(nlevels(fl), 0, 8)),
       u = list(ch4 = rnorm(nlevels(lg), 0, if (REALISTIC) 5 else 3), c2h6 = rnorm(nlevels(lg), 0, 0.2), co = rnorm(nlevels(lg), 0, 3)),
       s = c(ch4 = 5, c2h6 = 0.3, co = 5), phi = if (REALISTIC) 0.5 else 0.3)
}
simulate <- function(tr) {
  f <- as.integer(fl); l <- as.integer(lg)
  e_ch4 <- tr$tau[f] * as.vector(J$H %*% tr$alpha)
  e_c2 <- tr$tau[f] * as.vector(J$H %*% (tr$alpha * tr$r))
  e_co <- tr$gamma * tr$tau[f] * J$Hco
  # REALISTIC: Student-t(3) errors, i.e. occasional unresolved plumes several sd off, as the real
  # residuals show (26 Sep: v1 residual sd 8-27 ppb against s ~1 ppb)
  nz <- function(e, s) (if (REALISTIC) rt(length(e), 3) else rnorm(length(e))) * sqrt(s^2 + (tr$phi * e)^2)
  data.frame(flight = S$flight, leg_key = S$leg_key,
             ch4 = tr$b$ch4[f] + e_ch4 + tr$u$ch4[l] + nz(e_ch4, tr$s["ch4"]),
             c2h6 = tr$b$c2h6[f] + e_c2 + tr$u$c2h6[l] + nz(e_c2, tr$s["c2h6"]),
             co = tr$b$co[f] + e_co + tr$u$co[l] + nz(e_co, tr$s["co"]))
}
CONFIGS <- list(ch4 = character(0), ch4_c2h6 = "c2h6", ch4_c2h6_co = c("c2h6", "co"))
mod <- inv_model()
res <- list(); summ <- list(); idc <- list(); ido <- list(); idr <- list(); rgs <- list()
for (rep in seq_len(N_REP)) {
  tr <- truth_for(rep); obs <- simulate(tr)
  fos_true <- sum((tr$alpha * J$box_t_hr)[K %in% FOSSIL_COMPONENTS]) / sum(tr$alpha * J$box_t_hr)
  for (cf in names(CONFIGS)) {
    sd <- make_stan_data(obs, J, use = CONFIGS[[cf]], settings = st)
    fit <- fit_inversion(mod, sd, seed = rep)
    dg <- fit$diagnostic_summary(quiet = TRUE)
    A <- summarise_alpha(fit, K)
    A$truth <- tr$alpha[K]; A$covered <- A$truth >= A$q05 & A$truth <= A$q95
    A$error_reduction <- 1 - A$sd_log_alpha / st$la_sd[K]
    A$config <- cf; A$rep <- rep
    res[[length(res) + 1]] <- A
    dr <- fit$draws(c("alpha", "fossil_share"), format = "draws_matrix")
    cr <- if (all(c("og_basin", "og_urban") %in% K))
      cor(log(dr[, sprintf("alpha[%d]", which(K == "og_basin"))]), log(dr[, sprintf("alpha[%d]", which(K == "og_urban"))]))[1] else NA
    fs <- as.numeric(dr[, "fossil_share"])
    idf <- identifiability(fit, K, st$la_sd)
    idc[[length(idc) + 1]] <- cbind(rep = rep, config = cf, idf$component, truth = tr$alpha[K], covered = A$covered)
    ido[[length(ido) + 1]] <- cbind(rep = rep, config = cf, idf$overall)
    idr[[length(idr) + 1]] <- cbind(rep = rep, config = cf, corr_long(idf$corr))
    rg <- region_posterior(fit, K, J$region_t_hr, truth = tr$alpha)
    if (!is.null(rg)) rgs[[length(rgs) + 1]] <- cbind(rep = rep, config = cf, rg)
    summ[[length(summ) + 1]] <- data.frame(rep = rep, config = cf, coverage_90 = mean(A$covered),
      median_error_reduction = median(A$error_reduction), corr_og_basin_urban = round(cr, 2),
      fossil_true = round(fos_true, 3), fossil_q50 = round(median(fs), 3),
      fossil_q05 = round(quantile(fs, 0.05), 3), fossil_q95 = round(quantile(fs, 0.95), 3),
      divergent = sum(dg$num_divergent), max_rhat = round(max(fit$summary("alpha")$rhat), 3))
    cat(sprintf("rep %d %-12s coverage %.2f  error reduction %.2f  corr(basin,urban OG) %5.2f  fossil %.2f [%.2f-%.2f] true %.2f  div %d  rhat %.3f\n",
                rep, cf, mean(A$covered), median(A$error_reduction), cr, median(fs), quantile(fs, 0.05), quantile(fs, 0.95),
                fos_true, sum(dg$num_divergent), max(fit$summary("alpha")$rhat)))
  }
}
RES <- do.call(rbind, res); SUM <- do.call(rbind, summ); rownames(SUM) <- NULL
num <- vapply(RES, is.numeric, NA); RES[num] <- lapply(RES[num], function(x) signif(x, 4))
write.csv(RES, inv_file("osse_results.csv"), row.names = FALSE)
write.csv(SUM, inv_file("osse_summary.csv"), row.names = FALSE)

# identifiability: per component, overall, and what each tracer added relative to CH4 alone
IDC <- do.call(rbind, idc); IDO <- tracer_gain(do.call(rbind, ido)); IDR <- do.call(rbind, idr)
write.csv(IDC, inv_file("osse_identifiability_components.csv"), row.names = FALSE)
write.csv(IDO, inv_file("osse_identifiability_overall.csv"), row.names = FALSE)
write.csv(IDR, inv_file("osse_posterior_correlations.csv"), row.names = FALSE)
if (length(rgs)) {
  RGd <- do.call(rbind, rgs); write.csv(RGd, inv_file("osse_region_totals.csv"), row.names = FALSE)
  cat("\nregion totals: 90% coverage of the true total and fossil share, and median |error| in fossil share:\n")
  print(do.call(rbind, lapply(split(RGd, list(RGd$region, RGd$config), drop = TRUE), function(s)
    data.frame(region = s$region[1], config = s$config[1], cover_E = mean(s$covered_E), cover_share = mean(s$covered_share),
               share_abs_err = round(median(abs(s$fossil_share_q50 - s$fossil_share_true)), 3),
               share_90w = round(median(s$fossil_share_q95 - s$fossil_share_q05), 3)))), row.names = FALSE)
}
agg <- function(v) round(mean(v, na.rm = TRUE), 2)
cat("\nidentifiability (mean over replicates):\n")
print(do.call(rbind, lapply(names(CONFIGS), function(cf) { s <- IDO[IDO$config == cf, ]
  data.frame(config = cf, DOFS = agg(s$dofs), info_gain_bits = agg(s$info_gain_bits), dIG_vs_ch4 = agg(s$dIG_vs_ch4_bits),
             corr_fossil_biogenic = agg(s$corr_fossil_biogenic), sd_fossil_share = agg(s$sd_fossil_share),
             n_well = agg(s$n_well), n_partial = agg(s$n_partial), n_not = agg(s$n_not)) })), row.names = FALSE)
cat("\nregime by component (most common over replicates):\n")
RG <- tapply(IDC$regime, list(IDC$component, IDC$config), function(v) names(sort(table(v), decreasing = TRUE))[1])
print(RG[K, names(CONFIGS), drop = FALSE])

# error reduction by component and configuration, averaged over replicates
ER <- tapply(RES$error_reduction, list(RES$component, RES$config), mean)[K, names(CONFIGS), drop = FALSE]
cat("\nmean error reduction (1 - posterior/prior sd of log alpha):\n"); print(round(ER, 2))

png(file.path(RUN_DIR, "figures", paste0("osse_recovery_", INV_MODEL, ".png")), 1800, 700, res = 150)
par(mfrow = c(1, 2), mar = c(7, 4.5, 2, 1))
cols <- c(ch4 = "#999999", ch4_c2h6 = "#4575b4", ch4_c2h6_co = "#d73027")
barplot(t(ER), beside = TRUE, col = cols, las = 2, ylim = c(0, 1), ylab = "error reduction in log alpha",
        legend.text = c("CH4", "CH4 + C2H6", "CH4 + C2H6 + CO"), args.legend = list(x = "topright", bty = "n", cex = 0.8))
title(sprintf("(a) information gain by component (%s, %d replicates)", PRIOR, N_REP), cex.main = 0.9)
r1 <- RES[RES$rep == 1, ]; x <- match(r1$component, K) + c(ch4 = -0.2, ch4_c2h6 = 0, ch4_c2h6_co = 0.2)[r1$config]
plot(x, r1$q50, log = "y", ylim = range(c(r1$q05, r1$q95, r1$truth)), pch = 19, col = cols[r1$config], xaxt = "n",
     xlab = "", ylab = "alpha (posterior median, 90% interval)")
segments(x, r1$q05, x, r1$q95, col = cols[r1$config]); axis(1, seq_along(K), K, las = 2)
points(match(K, K), r1$truth[match(K, r1$component)], pch = 4, cex = 1.4, lwd = 2)
abline(h = 1, lty = 3, col = "gray50")
title("(b) replicate 1 (Paper-1-like truth, x = truth)", cex.main = 0.9)
invisible(dev.off())
cat("\nwrote osse_results.csv, osse_summary.csv, osse_identifiability_*.csv, osse_posterior_correlations.csv,",
    "figures/osse_recovery.png to", RUN_DIR, "\n")
