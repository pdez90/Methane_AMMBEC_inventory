# 12_endmember_sensitivity.R — how much do the ethane:methane endmembers control the answer? ------------------
# The fossil / non-fossil split rests on the ratio priors (basin gas 0.080 +- 0.015, delivered urban gas
# 0.110 +- 0.012, "other" 0.02 +- 0.02). Each case below refits the main configuration (CH4 + C2H6, the
# transport the environment selects) with the priors shifted, widened or collapsed, for every prior inventory:
#   base          as in config.R
#   basin_lo/hi   basin 0.065 / 0.095           urban_lo/hi  urban gas 0.090 / 0.130
#   no_contrast   basin = urban gas = 0.095 +- 0.015 (ethane cannot separate basin from urban gas)
#   sd_half/sd_x2 all ratio sds halved / doubled
#   other_zero    "other" fixed at 0 (treated as biogenic)     other_fossil  "other" 0.080 +- 0.015 (treated as gas)
#   Rscript scripts/12_endmember_sensitivity.R [prior names]     (env METHANE_EM_CASES=base,no_contrast,... to subset)
#   Out: <RUN_DIR>/endmember_sensitivity_<model>.csv, endmember_sensitivity_summary_<model>.csv
proj <- if (file.exists("config.R")) "." else ".."; source(file.path(proj, "config.R"))
source(file.path(proj, "R", "inversion.R")); source(file.path(proj, "R", "identifiability.R"))
ALL <- c("base", "basin_lo", "basin_hi", "urban_lo", "urban_hi", "no_contrast", "sd_half", "sd_x2", "other_zero", "other_fossil")
CASES <- strsplit(Sys.getenv("METHANE_EM_CASES", paste(ALL, collapse = ",")), ",")[[1]]
if (length(bad <- setdiff(CASES, ALL))) stop("unknown case(s): ", paste(bad, collapse = ", "))
priors <- commandArgs(TRUE)
if (!length(priors)) priors <- sub("^jacobian_(.*)\\.rds$", "\\1", list.files(RUN_DIR, pattern = "^jacobian_.*\\.rds$"))
S0 <- read.csv(file.path(INV_OUT, "obs_segments.csv"), stringsAsFactors = FALSE)
R0 <- R_PRIOR; mod <- inv_model(); out <- list()
set_case <- function(cs) {
  R <- R0; g <- names(R)
  if (cs == "basin_lo") R$og_basin["mean"] <- 0.065
  if (cs == "basin_hi") R$og_basin["mean"] <- 0.095
  if (cs %in% c("urban_lo", "urban_hi")) for (k in c("og_urban", "postmeter")) R[[k]]["mean"] <- if (cs == "urban_lo") 0.090 else 0.130
  if (cs == "no_contrast") for (k in c("og_basin", "og_urban", "postmeter")) R[[k]] <- c(mean = 0.095, sd = 0.015)
  if (cs == "sd_half") for (k in g) R[[k]]["sd"] <- R[[k]]["sd"] / 2
  if (cs == "sd_x2")   for (k in g) R[[k]]["sd"] <- R[[k]]["sd"] * 2
  if (cs == "other_zero")   R$other <- c(mean = 0, sd = 0)
  if (cs == "other_fossil") R$other <- c(mean = 0.080, sd = 0.015)
  R
}
cat(sprintf("endmember sensitivity: model %s, transport %s, cases %s\n", INV_MODEL, TRANSPORT, paste(CASES, collapse = " ")))
for (pr in priors) {
  J <- readRDS(file.path(RUN_DIR, paste0("jacobian_", pr, ".rds")))
  S <- S0[match(J$ids, S0$id), ]; S$leg_key <- paste(S$flight, S$leg_id)
  seen <- colSums(J$H) > 1e-6 * sum(J$H)
  J$H <- J$H[, seen, drop = FALSE]; K <- J$components <- J$components[seen]; J$box_t_hr <- J$box_t_hr[K]
  la_sd <- inv_settings(K)$la_sd
  for (cs in CASES) {
    R_PRIOR <<- set_case(cs)
    sd <- make_stan_data(S, J, use = "c2h6")
    t0 <- Sys.time(); fit <- fit_inversion(mod, sd)
    dr <- fit$draws(c("E_total", "fossil_share"), format = "draws_matrix")
    rg <- region_posterior(fit, K, J$region_t_hr); rv <- function(r, col) if (!is.null(rg) && r %in% rg$region) rg[[col]][rg$region == r] else NA_real_
    A <- fit$draws("alpha", format = "draws_matrix")[, sprintf("alpha[%d]", seq_along(K)), drop = FALSE]
    ob <- if ("og_basin" %in% K && !is.null(J$region_t_hr)) median(A[, which(K == "og_basin")]) * J$region_t_hr["og_basin", "djb"] else NA
    io <- identifiability(fit, K, la_sd)$overall
    rr <- fit$summary("r")$median; names(rr) <- K
    dg <- fit$diagnostic_summary(quiet = TRUE)
    row <- data.frame(prior = pr, case = cs, transport = TRANSPORT,
      r_basin = round(unname(rr["og_basin"]), 3), r_urban = round(unname(rr["og_urban"]), 3), r_other = round(unname(rr["other"]), 3),
      box_q50 = round(median(dr[, "E_total"]), 2), box_fossil = round(median(dr[, "fossil_share"]), 3),
      obs_box_q50 = rv("obs_box", "E_q50"), obs_box_fossil = rv("obs_box", "fossil_share_q50"),
      djb_q50 = rv("djb", "E_q50"), djb_fossil = rv("djb", "fossil_share_q50"), djb_og_t_hr = round(ob, 2),
      info_gain_bits = io$info_gain_bits, sd_fossil_share = io$sd_fossil_share,
      divergent = sum(dg$num_divergent), max_rhat = round(max(fit$summary("alpha")$rhat), 3), row.names = NULL)
    out[[length(out) + 1]] <- row
    with(row, cat(sprintf("%-17s %-12s box %.1f fossil %.2f | domain %.1f fossil %.2f | DJB %.1f fossil %.2f og %.1f | r %.3f/%.3f/%.3f | IG %.1f div %d (%.1f min)\n",
      prior, case, box_q50, box_fossil, obs_box_q50, obs_box_fossil, djb_q50, djb_fossil, djb_og_t_hr, r_basin, r_urban, r_other, info_gain_bits, divergent,
      as.numeric(difftime(Sys.time(), t0, units = "mins")))))
  }
}
R_PRIOR <<- R0
B <- do.call(rbind, out); write.csv(B, inv_file("endmember_sensitivity.csv"), row.names = FALSE)
rng <- function(v) sprintf("%.2f-%.2f", min(v, na.rm = TRUE), max(v, na.rm = TRUE))
SUM <- do.call(rbind, lapply(split(B, B$prior), function(s) data.frame(prior = s$prior[1], n_cases = nrow(s),
  box_fossil_range = rng(s$box_fossil), obs_box_fossil_range = rng(s$obs_box_fossil), djb_fossil_range = rng(s$djb_fossil),
  djb_og_range_t_hr = rng(s$djb_og_t_hr), box_total_range = rng(s$box_q50), obs_box_total_range = rng(s$obs_box_q50),
  fossil_no_contrast = s$obs_box_fossil[s$case == "no_contrast"][1], fossil_other_fossil = s$obs_box_fossil[s$case == "other_fossil"][1])))
write.csv(SUM, inv_file("endmember_sensitivity_summary.csv"), row.names = FALSE)
cat("\nrange across endmember cases, per prior:\n"); print(SUM, row.names = FALSE)
cat("wrote endmember_sensitivity_*.csv to", RUN_DIR, "\n")
