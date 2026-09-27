# 13_predictive_checks.R — does each prior structure actually predict the tracers? -----------------------------
# Reads the saved fits of 07 (fit_<prior>_<config>_<model>.rds in RUN_DIR) and, from the posterior draws:
#   (A) posterior predictive means for CH4 and C2H6 at every segment; residual sd, correlation between predicted
#       and observed enhancement, and the fraction of the observed C2H6 enhancement explained, by region;
#   (B) the pointwise log-likelihood matrix and out-of-sample predictive accuracy: PSIS-LOO (package loo, if
#       installed) and WAIC (always). Priors are compared at the same tracer set only — the data are then
#       identical, so elpd differences are meaningful; across tracer sets they are not (different data);
#   (C) coverage of the OSSE (06) 90% intervals for region totals and fossil shares, by tracer set.
#   Rscript scripts/13_predictive_checks.R          (env METHANE_ZISCALE selects the transport, as for 07)
#   Out: <RUN_DIR>/predictive_checks_<model>.csv, model_comparison_<model>.csv, osse_coverage_<model>.csv
proj <- if (file.exists("config.R")) "." else ".."; source(file.path(proj, "config.R"))
source(file.path(proj, "R", "inversion.R"))
N_DRAWS <- as.integer(Sys.getenv("METHANE_PPC_DRAWS", "400"))
S0 <- read.csv(file.path(INV_OUT, "obs_segments.csv"), stringsAsFactors = FALSE)
# main run plus the prior variants that live in sub-directories (v11waste_epa, v11waste_v2, metro_wwtp): each is a
# separate source representation and is compared on the same data
DIRS <- c(main = RUN_DIR, setNames(list.dirs(RUN_DIR, recursive = FALSE), basename(list.dirs(RUN_DIR, recursive = FALSE))))
DIRS <- DIRS[names(DIRS) %in% c("main", "v11waste_epa", "v11waste_v2", "metro_wwtp")]
fits <- unlist(lapply(names(DIRS), function(v) { f <- list.files(DIRS[[v]], pattern = sprintf("^fit_.*_%s\\.rds$", INV_MODEL)); if (length(f)) paste(v, f, sep = "|") else NULL }))
if (!length(fits)) stop("no fit_*_", INV_MODEL, ".rds in ", RUN_DIR)
has_loo <- requireNamespace("loo", quietly = TRUE)
if (!has_loo) cat("package 'loo' not installed: reporting WAIC only (install.packages('loo') for PSIS-LOO)\n")
lse <- function(v) { m <- max(v); m + log(mean(exp(v - m))) }
dt_ll <- function(y, mu, sc, nu) dt((y - mu) / sc, df = nu, log = TRUE) - log(sc)
PC <- list(); MC <- list()
LL <- list()
for (f0 in fits) {
  v <- sub("\\|.*$", "", f0); f <- sub("^.*\\|", "", f0); RD <- DIRS[[v]]
  m <- regmatches(f, regexec(sprintf("^fit_(.*)_(ch4(?:_c2h6)?(?:_co)?)_%s\\.rds$", INV_MODEL), f, perl = TRUE))[[1]]
  if (!length(m)) next
  pr <- if (v == "main") m[2] else paste0(m[2], "+", sub("^v11waste_", "waste_", v)); cf <- m[3]; use <- list(ch4 = character(0), ch4_c2h6 = "c2h6", ch4_c2h6_co = c("c2h6", "co"))[[cf]]
  J <- readRDS(file.path(RD, paste0("jacobian_", m[2], ".rds")))
  S <- S0[match(J$ids, S0$id), ]; S$leg_key <- paste(S$flight, S$leg_id)
  seen <- colSums(J$H) > 1e-6 * sum(J$H)
  J$H <- J$H[, seen, drop = FALSE]; K <- J$components <- J$components[seen]; J$box_t_hr <- J$box_t_hr[K]
  sd <- make_stan_data(S, J, use = use)
  fit <- readRDS(file.path(RD, f))
  D <- fit$draws(c("alpha", "r", "tau", "b_ch4_z", "b_c2h6_z", "sl_ch4", "sl_c2h6", "u_ch4_z", "u_c2h6_z", "s_ch4", "s_c2h6", "phi"), format = "draws_matrix")
  if (nrow(D) > N_DRAWS) D <- D[round(seq(1, nrow(D), length.out = N_DRAWS)), , drop = FALSE]
  D <- unclass(D); attr(D, "nchains") <- NULL                       # plain matrix: column subsets drop to vectors
  g <- function(v, n) D[, sprintf("%s[%d]", v, seq_len(n)), drop = FALSE]
  al <- g("alpha", sd$K); rr <- g("r", sd$K); ta <- g("tau", sd$F)
  Hal <- al %*% t(sd$H); Hr <- (al * rr) %*% t(sd$H)                                   # draws x N
  tf <- ta[, sd$flight, drop = FALSE]
  enh4 <- tf * Hal; mu4 <- matrix(sd$b_ch4_mu[sd$flight], nrow(D), sd$N, byrow = TRUE) + sd$b_ch4_sd * g("b_ch4_z", sd$F)[, sd$flight] + enh4 + D[, "sl_ch4"] * g("u_ch4_z", sd$L)[, sd$leg]
  sc4 <- sqrt(D[, "s_ch4"]^2 + (D[, "phi"] * enh4)^2)
  ll4 <- dt_ll(matrix(sd$y_ch4, nrow(D), sd$N, byrow = TRUE), mu4, sc4, sd$nu)
  ll <- ll4; obs_all <- list(ch4 = sd$y_ch4); pred_all <- list(ch4 = colMeans(mu4)); bg_all <- list(ch4 = sd$b_ch4_mu[sd$flight])
  if (sd$N_c2h6 > 0) {
    i <- sd$idx_c2h6; enh2 <- tf[, i] * Hr[, i]
    mu2 <- matrix(sd$b_c2h6_mu[sd$flight[i]], nrow(D), length(i), byrow = TRUE) + sd$b_c2h6_sd * g("b_c2h6_z", sd$F)[, sd$flight[i]] + enh2 + D[, "sl_c2h6"] * g("u_c2h6_z", sd$L)[, sd$leg[i]]
    sc2 <- sqrt(D[, "s_c2h6"]^2 + (D[, "phi"] * enh2)^2)
    ll2 <- dt_ll(matrix(sd$y_c2h6, nrow(D), length(i), byrow = TRUE), mu2, sc2, sd$nu)
    ll <- cbind(ll, ll2); obs_all$c2h6 <- sd$y_c2h6; pred_all$c2h6 <- colMeans(mu2); bg_all$c2h6 <- sd$b_c2h6_mu[sd$flight[i]]
    reg2 <- S$region[i]
  }
  # (A) per region, per tracer
  for (tr in names(obs_all)) {
    reg <- if (tr == "ch4") S$region else reg2
    for (r in c("urban", "edge", "basin")) { k <- reg == r; if (!any(k)) next
      eo <- obs_all[[tr]][k] - bg_all[[tr]][k]; ep <- pred_all[[tr]][k] - bg_all[[tr]][k]; res <- obs_all[[tr]][k] - pred_all[[tr]][k]
      PC[[length(PC) + 1]] <- data.frame(prior = pr, config = cf, tracer = tr, region = r, n = sum(k),
        obs_enh_median = round(median(eo), 3), pred_enh_median = round(median(ep), 3), resid_sd = round(sd(res), 3),
        resid_median_abs = round(median(abs(res)), 3), cor_pred_obs_enh = round(cor(ep, eo), 3),
        share_explained = round(sum(pmax(ep, 0)) / sum(pmax(eo, 0)), 3)) }
  }
  # (B) elpd
  lppd <- sum(apply(ll, 2, lse)); p_waic <- sum(apply(ll, 2, var)); elpd_waic <- lppd - p_waic
  LL[[paste(pr, cf)]] <- ll
  row <- data.frame(prior = pr, config = cf, n_tracer_obs = ncol(ll), n_ch4 = sd$N, n_c2h6 = sd$N_c2h6, lppd = round(lppd, 1), p_waic = round(p_waic, 1), elpd_waic = round(elpd_waic, 1),
                    elpd_loo = NA, se_elpd_loo = NA, p_loo = NA, pareto_k_gt_0.7 = NA)
  if (has_loo) { lo <- tryCatch(loo::loo(ll, r_eff = NA), error = function(e) NULL)
    if (!is.null(lo)) { row$elpd_loo <- round(lo$estimates["elpd_loo", "Estimate"], 1); row$se_elpd_loo <- round(lo$estimates["elpd_loo", "SE"], 1)
      row$p_loo <- round(lo$estimates["p_loo", "Estimate"], 1); row$pareto_k_gt_0.7 <- sum(lo$diagnostics$pareto_k > 0.7) } }
  MC[[length(MC) + 1]] <- row
  cat(sprintf("%-17s %-12s elpd_waic %8.1f  p_waic %5.1f  elpd_loo %8s  (k>0.7: %s)\n", pr, cf, elpd_waic, p_waic, ifelse(is.na(row$elpd_loo), "-", row$elpd_loo), ifelse(is.na(row$pareto_k_gt_0.7), "-", row$pareto_k_gt_0.7)))
}
PC <- do.call(rbind, PC); MC <- do.call(rbind, MC)
write.csv(PC, inv_file("predictive_checks.csv"), row.names = FALSE); write.csv(MC, inv_file("model_comparison.csv"), row.names = FALSE)
cat("\nposterior predictive checks (urban segments):\n"); print(PC[PC$region == "urban", ], row.names = FALSE)
cat("\nmodel comparison within each tracer set (same tracer observations within a set; units = CH4 + C2H6 [+ CO excluded] segment observations):\n")
CMP <- list()
for (cf in unique(MC$config)) { m <- MC[MC$config == cf, ]; e <- if (all(!is.na(m$elpd_loo))) m$elpd_loo else m$elpd_waic
  cat(sprintf("  %-12s %s\n", cf, paste(sprintf("%s %+.1f", m$prior, e - max(e)), collapse = " | ")))
  if (has_loo && nrow(m) >= 2) {   # paired elpd differences with their standard error (loo_compare)
    lo <- lapply(m$prior, function(p) tryCatch(loo::loo(LL[[paste(p, cf)]], r_eff = NA), error = function(e) NULL)); names(lo) <- m$prior
    lo <- lo[!vapply(lo, is.null, NA)]
    if (length(lo) >= 2) { cc <- loo::loo_compare(lo); cc <- as.data.frame(cc[, c("elpd_diff", "se_diff")]); cc$prior <- rownames(cc); cc$config <- cf
      CMP[[cf]] <- cc; cat("    paired LOO differences (elpd_diff +- se_diff, vs the best):\n"); for (i in seq_len(nrow(cc))) cat(sprintf("      %-22s %+7.1f +- %.1f\n", cc$prior[i], cc$elpd_diff[i], cc$se_diff[i])) } }
}
if (length(CMP)) write.csv(do.call(rbind, CMP), inv_file("model_comparison_paired.csv"), row.names = FALSE)

## (C) OSSE coverage --------------------------------------------------------------------------------------------
of <- file.path(TRANSPORT_DIR, sprintf("osse_region_totals_%s.csv", INV_MODEL))
if (!file.exists(of)) of <- file.path(INV_OUT, "runs", "np1000_h10", sprintf("osse_region_totals_%s.csv", INV_MODEL))
if (file.exists(of)) {
  O <- read.csv(of, stringsAsFactors = FALSE)
  CV <- aggregate(cbind(covered_E, covered_share) ~ config + region, data = O, FUN = function(v) round(mean(v, na.rm = TRUE), 2))
  nrep <- aggregate(rep ~ config + region, data = O, FUN = function(v) length(unique(v)))
  CV$n_rep <- nrep$rep[match(paste(CV$config, CV$region), paste(nrep$config, nrep$region))]
  err <- aggregate(cbind(abs_err_share = abs(fossil_share_q50 - fossil_share_true), rel_err_E = abs(E_q50 - E_true) / E_true) ~ config + region, data = O, FUN = function(v) round(median(v), 3))
  CV <- merge(CV, err); write.csv(CV, inv_file("osse_coverage.csv"), row.names = FALSE)
  cat(sprintf("\nOSSE coverage of 90%% posterior intervals (%s, %d replicates):\n", basename(dirname(of)), max(CV$n_rep))); print(CV, row.names = FALSE)
}
cat("wrote predictive_checks, model_comparison and osse_coverage csv to", RUN_DIR, "\n")
