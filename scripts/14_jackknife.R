# 14_jackknife.R — leave-one-flight-out sampling uncertainty of the totals and fossil shares ---------------------
# The sensitivity ranges of 09 measure how the answer moves across model configurations; they say nothing about
# how much it depends on which flights happened to be flown. Refitting the main configuration (CH4 + C2H6, the
# transport the environment selects) with each flight left out gives a jackknife standard error for every
# reported quantity: se = sqrt((n-1)/n * sum((theta_i - mean)^2)). ~1 min per fit; 22 flights x priors.
#   Rscript scripts/14_jackknife.R [prior names]         (env METHANE_JK_MAXFLIGHTS to test on a few)
#   Out: <RUN_DIR>/jackknife_<model>.csv (every refit), jackknife_summary_<model>.csv
proj <- if (file.exists("config.R")) "." else ".."; source(file.path(proj, "config.R"))
source(file.path(proj, "R", "inversion.R")); source(file.path(proj, "R", "identifiability.R"))
priors <- commandArgs(TRUE)
if (!length(priors)) priors <- sub("^jacobian_(.*)\\.rds$", "\\1", list.files(RUN_DIR, pattern = "^jacobian_.*\\.rds$"))
MAXF <- as.integer(Sys.getenv("METHANE_JK_MAXFLIGHTS", "999"))
S0 <- read.csv(file.path(INV_OUT, "obs_segments.csv"), stringsAsFactors = FALSE)
mod <- inv_model(); out <- list()
for (pr in priors) {
  J0 <- readRDS(file.path(RUN_DIR, paste0("jacobian_", pr, ".rds")))
  S <- S0[match(J0$ids, S0$id), ]; S$leg_key <- paste(S$flight, S$leg_id)
  seen <- colSums(J0$H) > 1e-6 * sum(J0$H)
  J0$H <- J0$H[, seen, drop = FALSE]; K <- J0$components <- J0$components[seen]; J0$box_t_hr <- J0$box_t_hr[K]
  flights <- head(sort(unique(S$flight)), MAXF)
  for (fl in c("none", flights)) {
    keep <- if (fl == "none") rep(TRUE, nrow(S)) else S$flight != fl
    J <- J0; J$H <- J0$H[keep, , drop = FALSE]; J$ids <- J0$ids[keep]; if (!is.null(J0$Hco)) J$Hco <- J0$Hco[keep]
    Sx <- S[keep, ]
    t0 <- Sys.time(); fit <- fit_inversion(mod, make_stan_data(Sx, J, use = "c2h6"))
    dr <- fit$draws(c("E_total", "fossil_share"), format = "draws_matrix")
    rg <- region_posterior(fit, K, J$region_t_hr); rv <- function(r, col) if (!is.null(rg) && r %in% rg$region) rg[[col]][rg$region == r] else NA_real_
    A <- fit$draws("alpha", format = "draws_matrix")[, sprintf("alpha[%d]", seq_along(K)), drop = FALSE]
    ob <- if ("og_basin" %in% K) median(A[, which(K == "og_basin")]) * J$region_t_hr["og_basin", "djb"] else NA
    dg <- fit$diagnostic_summary(quiet = TRUE)
    row <- data.frame(prior = pr, left_out = fl, n_obs = nrow(Sx), transport = TRANSPORT,
      box = round(median(dr[, "E_total"]), 3), box_fossil = round(median(dr[, "fossil_share"]), 3),
      obs_box = rv("obs_box", "E_q50"), obs_box_fossil = rv("obs_box", "fossil_share_q50"),
      djb = rv("djb", "E_q50"), djb_fossil = rv("djb", "fossil_share_q50"), djb_og = round(ob, 3),
      divergent = sum(dg$num_divergent), max_rhat = round(max(fit$summary("alpha")$rhat), 3), row.names = NULL)
    out[[length(out) + 1]] <- row
    with(row, cat(sprintf("%-17s -%-14s box %.2f f %.2f | domain %.1f f %.2f | DJB %.1f f %.2f og %.1f | div %d (%.1f min)\n",
      prior, fl, box, box_fossil, obs_box, obs_box_fossil, djb, djb_fossil, djb_og, divergent, as.numeric(difftime(Sys.time(), t0, units = "mins")))))
  }
}
B <- do.call(rbind, out); write.csv(B, inv_file("jackknife.csv"), row.names = FALSE)
qs <- c("box", "box_fossil", "obs_box", "obs_box_fossil", "djb", "djb_fossil", "djb_og")
SUM <- do.call(rbind, lapply(split(B, B$prior), function(s) { full <- s[s$left_out == "none", ]; jk <- s[s$left_out != "none", ]; n <- nrow(jk)
  do.call(rbind, lapply(qs, function(q) { v <- jk[[q]]; se <- sqrt((n - 1) / n * sum((v - mean(v))^2))
    data.frame(prior = s$prior[1], quantity = q, full = full[[q]], n_flights = n, jackknife_se = round(se, 3),
               jackknife_rel_se = round(se / full[[q]], 3), min_leave_one_out = round(min(v), 3), max_leave_one_out = round(max(v), 3),
               most_influential_flight = jk$left_out[which.max(abs(v - full[[q]]))]) })) }))
write.csv(SUM, inv_file("jackknife_summary.csv"), row.names = FALSE)
cat("\nleave-one-flight-out jackknife (main configuration):\n"); print(SUM, row.names = FALSE)
cat("wrote jackknife_*.csv to", RUN_DIR, "\n")
