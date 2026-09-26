# 07_invert.R — the inversion on the measured segments ------------------------------------------
# Fits stan/inversion.stan to the observed CH4 / C2H6 / CO segment means for each prior (GRA2PES
# v1.1, v2.0beta, EPA) and each tracer set (CH4; +C2H6; +C2H6+CO). Priors that agree after the
# fit mean the data, not the prior, set the answer. Read the OSSE (06) first: it says which
# components these flights can constrain at all.
#
#   Rscript scripts/07_invert.R [prior ...]           (default: every jacobian_*.rds)
#   Out: <RUN_DIR>/posterior_alpha.csv, posterior_summary.csv, fit_<prior>_<config>.rds,
#        posterior_identifiability_{components,overall}.csv, posterior_correlations.csv,
#        posterior_prior_dependence.csv (how much the answer still depends on which inventory)
# ----------------------------------------------------------------------------------------------
proj <- if (file.exists("config.R")) "." else ".."; source(file.path(proj, "config.R"))
source(file.path(proj, "R", "inversion.R")); source(file.path(proj, "R", "identifiability.R")); source(file.path(proj, "R", "diagnostics.R"))
cat("model", INV_MODEL, "(set METHANE_INV_MODEL=v1 for the 24 Sep Gaussian model)\n")
DC_all <- list()
if (nzchar(Sys.getenv("CMDSTAN"))) set_cmdstan_path(Sys.getenv("CMDSTAN"))

priors <- commandArgs(TRUE)
if (!length(priors)) priors <- sub("^jacobian_(.*)\\.rds$", "\\1", list.files(RUN_DIR, pattern = "^jacobian_.*\\.rds$"))
S0 <- read.csv(file.path(INV_OUT, "obs_segments.csv"), stringsAsFactors = FALSE)
CONFIGS <- list(ch4 = character(0), ch4_c2h6 = "c2h6", ch4_c2h6_co = c("c2h6", "co"))
mod <- inv_model(); A_all <- list(); S_all <- list(); IC_all <- list(); IO_all <- list(); IR_all <- list(); RG_all <- list()
for (pr in priors) {
  J <- readRDS(file.path(RUN_DIR, paste0("jacobian_", pr, ".rds")))
  S <- S0[match(J$ids, S0$id), ]; S$leg_key <- paste(S$flight, S$leg_id)
  seen <- colSums(J$H) > 1e-6 * sum(J$H)
  J$H <- J$H[, seen, drop = FALSE]; K <- J$components <- J$components[seen]; J$box_t_hr <- J$box_t_hr[K]
  for (cf in names(CONFIGS)) {
    if ("co" %in% CONFIGS[[cf]] && !isTRUE(J$has_co)) next
    sd <- make_stan_data(S, J, use = CONFIGS[[cf]])
    fit <- fit_inversion(mod, sd)
    fit$save_object(inv_file(sprintf("fit_%s_%s.rds", pr, cf)))
    A <- summarise_alpha(fit, K); A$prior <- pr; A$config <- cf; A$E_prior <- J$box_t_hr[K]
    A_all[[length(A_all) + 1]] <- A
    idf <- identifiability(fit, K, inv_settings(K)$la_sd)
    IC_all[[length(IC_all) + 1]] <- cbind(prior = pr, config = cf, idf$component)
    IO_all[[length(IO_all) + 1]] <- cbind(prior = pr, config = cf, idf$overall,
      E_fossil_prior = round(sum(J$box_t_hr[K %in% FOSSIL_COMPONENTS]), 3), E_total_prior = round(sum(J$box_t_hr), 3))
    IR_all[[length(IR_all) + 1]] <- cbind(prior = pr, config = cf, corr_long(idf$corr))
    rg <- region_posterior(fit, K, J$region_t_hr)
    if (!is.null(rg)) RG_all[[length(RG_all) + 1]] <- cbind(prior = pr, config = cf, rg)
    dc <- fit_decomposition(fit, S, sd); print_decomposition(dc, paste(pr, cf))
    DC_all[[length(DC_all) + 1]] <- cbind(prior = pr, config = cf, dc$by_region)
    dr <- fit$draws(c("E_total", "fossil_share"), format = "draws_matrix")
    dg <- fit$diagnostic_summary(quiet = TRUE)
    S_all[[length(S_all) + 1]] <- data.frame(prior = pr, config = cf, n_obs = nrow(S),
      E_prior = round(sum(J$box_t_hr), 2), E_q05 = quantile(dr[, "E_total"], 0.05), E_q50 = median(dr[, "E_total"]),
      E_q95 = quantile(dr[, "E_total"], 0.95), fossil_q05 = quantile(dr[, "fossil_share"], 0.05),
      fossil_q50 = median(dr[, "fossil_share"]), fossil_q95 = quantile(dr[, "fossil_share"], 0.95),
      divergent = sum(dg$num_divergent), max_rhat = max(fit$summary("alpha")$rhat))
    with(S_all[[length(S_all)]], cat(sprintf("%-18s %-12s box total %.1f t/h [%.1f-%.1f] (prior %.1f); fossil %.2f [%.2f-%.2f]; div %d, rhat %.3f\n",
      pr, cf, E_q50, E_q05, E_q95, E_prior, fossil_q50, fossil_q05, fossil_q95, divergent, max_rhat)))
  }
}
A <- do.call(rbind, A_all); SM <- do.call(rbind, S_all); rownames(SM) <- NULL
num <- vapply(A, is.numeric, NA); A[num] <- lapply(A[num], function(x) signif(x, 4))
num <- vapply(SM, is.numeric, NA); SM[num] <- lapply(SM[num], function(x) signif(x, 4))
write.csv(A, inv_file("posterior_alpha.csv"), row.names = FALSE)
write.csv(SM, inv_file("posterior_summary.csv"), row.names = FALSE)
# ---- identifiability and prior dependence ----------------------------------------------------
IC <- do.call(rbind, IC_all); IR <- do.call(rbind, IR_all)
IO <- tracer_gain(do.call(rbind, IO_all), group = "prior")
write.csv(IC, inv_file("posterior_identifiability_components.csv"), row.names = FALSE)
write.csv(IO, inv_file("posterior_identifiability_overall.csv"), row.names = FALSE)
write.csv(IR, inv_file("posterior_correlations.csv"), row.names = FALSE)
# Prior dependence: the three inventories disagree ~10-fold on the box total. If the posteriors
# agree, the data set the answer. Ratio = spread (sd of log) of posterior medians across priors /
# spread of the prior values; 1 = the data changed nothing, 0 = fully data-determined.
if (length(unique(SM$prior)) >= 2) {
  PD <- do.call(rbind, lapply(split(SM, SM$config), function(s) {
    lp <- function(v) sd(log(pmax(v, 1e-9)))
    pf <- IO$E_fossil_prior[match(s$prior, IO$prior)] / IO$E_total_prior[match(s$prior, IO$prior)]
    data.frame(config = s$config[1], n_priors = nrow(s),
      total_prior_spread = round(lp(s$E_prior), 3), total_post_spread = round(lp(s$E_q50), 3),
      total_prior_dependence = round(lp(s$E_q50) / lp(s$E_prior), 2),
      fossil_share_prior_range = paste(round(range(pf), 2), collapse = "-"),
      fossil_share_post_range = paste(round(range(s$fossil_q50), 2), collapse = "-"),
      fossil_share_prior_dependence = round(diff(range(s$fossil_q50)) / max(diff(range(pf)), 1e-9), 2))
  }))
  write.csv(PD, inv_file("posterior_prior_dependence.csv"), row.names = FALSE)
  cat("\nprior dependence across inventories (0 = data-determined, 1 = prior-determined):\n"); print(PD, row.names = FALSE)
}
if (length(RG_all)) {
  write.csv(do.call(rbind, DC_all), inv_file("posterior_decomposition.csv"), row.names = FALSE)
  RG <- do.call(rbind, RG_all); write.csv(RG, inv_file("posterior_region_totals.csv"), row.names = FALSE)
  cat("\nposterior totals by region (t CH4/h, median [90%]) and fossil share:\n")
  print(within(RG, { E = sprintf("%.2f [%.2f-%.2f]", E_q50, E_q05, E_q95); fossil = sprintf("%.2f [%.2f-%.2f]", fossil_share_q50, fossil_share_q05, fossil_share_q95) })[
        , c("prior", "config", "region", "E_prior", "E", "fossil")], row.names = FALSE)
}
cat("\nidentifiability by prior and tracer set:\n")
print(IO[, c("prior", "config", "dofs", "info_gain_bits", "dIG_vs_ch4_bits", "corr_fossil_biogenic",
             "sd_fossil_share", "sd_fossil_share_ratio_vs_ch4", "n_well", "n_partial", "n_not")], row.names = FALSE)
cat("wrote posterior_*_", INV_MODEL, ".csv (alpha, summary, region_totals, decomposition, identifiability, correlations,",
    " prior_dependence) to ", RUN_DIR, "\n", sep = "")
