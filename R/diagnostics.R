# diagnostics.R — where does the observed enhancement go in a fit? ---------------------------------
# Splits every segment's observed CH4 into flight background, source signal (tau * H alpha), leg
# offset and residual, using posterior medians, and summarises by region. If the backgrounds or
# leg offsets carry much of what the flights saw above their 5th percentile, the alphas are biased
# low: that is what the v1 model did on the real data (26 Sep).
fit_decomposition <- function(fit, S, sd) {
  med <- function(v) apply(fit$draws(v, format = "draws_matrix"), 2, median)
  dr <- fit$draws(c("alpha", "tau"), format = "draws_matrix")
  al <- dr[, grep("^alpha\\[", colnames(dr)), drop = FALSE]; ta <- dr[, grep("^tau\\[", colnames(dr)), drop = FALSE]
  sig <- vapply(seq_len(sd$N), function(i) median(ta[, sd$flight[i]] * (al %*% sd$H[i, ])), 0)
  u <- fit$draws(c("sl_ch4", "u_ch4_z"), format = "draws_matrix")
  leg <- vapply(seq_len(sd$N), function(i) median(u[, "sl_ch4"] * u[, sprintf("u_ch4_z[%d]", sd$leg[i])]), 0)
  bsh <- sd$b_ch4_sd * med("b_ch4_z")                      # background shift from the prior centre, ppb
  enh_obs <- sd$y_ch4 - sd$b_ch4_mu[sd$flight]
  res <- enh_obs - bsh[sd$flight] - sig - leg
  nuis <- med(c("s_ch4", "phi", "sl_ch4"))
  by <- do.call(rbind, lapply(c("urban", "edge", "basin"), function(r) { k <- S$region == r
    data.frame(region = r, n = sum(k), obs_above_p05 = median(enh_obs[k]), prior_signal = median(rowSums(sd$H[k, , drop = FALSE])),
               fitted_signal = median(sig[k]), bg_shift = median(bsh[sd$flight][k]), leg_offset = median(leg[k]),
               abs_leg = median(abs(leg[k])), resid_sd = sd(res[k]),
               signal_share = sum(sig[k]) / sum(pmax(enh_obs[k], 0))) }))
  list(by_region = by, nuisance = nuis, bg_shift = bsh)
}
print_decomposition <- function(dc, label) {
  b <- dc$by_region
  cat(sprintf("  %s: s %.1f ppb, phi %.2f, leg scale %.1f ppb, background shift median %+.1f ppb (range %+.1f to %+.1f)\n", label,
              dc$nuisance["s_ch4"], dc$nuisance["phi"], dc$nuisance["sl_ch4"], median(dc$bg_shift), min(dc$bg_shift), max(dc$bg_shift)))
  for (i in seq_len(nrow(b))) cat(sprintf("    %-5s obs-p05 %5.1f | fitted signal %5.1f (%.0f%% of enhancement) | bg %+5.1f | leg |%4.1f| | resid sd %4.1f\n",
    b$region[i], b$obs_above_p05[i], b$fitted_signal[i], 100 * b$signal_share[i], b$bg_shift[i], b$abs_leg[i], b$resid_sd[i]))
}
