# identifiability.R — how much do the observations (and each co-tracer) tell us? --------------------
# Paper 2's framing: which source components can CH4 alone identify, and what does adding C2H6
# (and CO) change? Computed from the posterior draws of one fit, in log(alpha) space, where the
# prior is Gaussian (log alpha_k ~ N(0, la_sd_k)) so the Gaussian formulas below are exact for the
# prior and a second-moment approximation for the posterior.
#
#   error reduction_k     1 - sd_post / sd_prior of log alpha_k
#   DOFS                  degrees of freedom for signal = K - tr(S_prior^-1 S_post)   (Rodgers 2000)
#   information gain      H(prior) - H(posterior), H = 0.5 log det(2 pi e S), in bits
#                         -> IG(C2H6) = IG(ch4_c2h6) - IG(ch4) is "what measuring ethane added"
#   fossil/biogenic corr  posterior correlation of log E_fossil and log E_biogenic (box totals):
#                         strongly negative = the data fix the total but cannot split it
#   regime_k              well identified / partially identified / not identified
# ----------------------------------------------------------------------------------------------
REGIME_ER_WELL <- 0.5      # error reduction needed to call a component well identified
REGIME_ER_MIN  <- 0.2      # below this it is not identified
REGIME_CORR    <- 0.5      # |posterior corr| with another component that makes it only partial

.bits <- 1 / log(2)
.gauss_H <- function(S) 0.5 * as.numeric(determinant(2 * pi * exp(1) * S, logarithm = TRUE)$modulus) * .bits

#' Identifiability of one fit. fit: CmdStanMCMC; K: component names; la_sd: prior sd of log alpha.
identifiability <- function(fit, K, la_sd, fossil = FOSSIL_COMPONENTS) {
  dr <- fit$draws(c("alpha", "E_post", "E_total", "fossil_share"), format = "draws_matrix")
  la <- log(dr[, sprintf("alpha[%d]", seq_along(K)), drop = FALSE]); colnames(la) <- K
  Sp <- stats::cov(la); S0 <- diag(la_sd[K]^2, length(K))
  C <- stats::cor(la); diag(C) <- NA
  er <- 1 - sqrt(diag(Sp)) / la_sd[K]
  maxc <- apply(abs(C), 1, max, na.rm = TRUE); maxc[!is.finite(maxc)] <- 0
  partner <- K[apply(abs(C), 1, function(v) if (all(is.na(v))) NA_integer_ else which.max(v))]
  regime <- ifelse(er < REGIME_ER_MIN, "not identified",
            ifelse(er >= REGIME_ER_WELL & maxc < REGIME_CORR, "well identified", "partially identified"))
  E <- dr[, sprintf("E_post[%d]", seq_along(K)), drop = FALSE]
  Ef <- rowSums(E[, K %in% fossil, drop = FALSE]); Eb <- rowSums(E[, !(K %in% fossil), drop = FALSE])
  lc <- function(a, b) if (all(a > 0) && all(b > 0) && sd(log(a)) > 0 && sd(log(b)) > 0) cor(log(a), log(b)) else NA_real_
  list(
    component = data.frame(component = K, error_reduction = round(er, 3),
      averaging_kernel_diag = round(1 - diag(Sp) / la_sd[K]^2, 3),
      max_abs_corr = round(maxc, 2), most_correlated_with = partner, regime = regime, row.names = NULL),
    overall = data.frame(K = length(K),
      dofs = round(length(K) - sum(diag(solve(S0, Sp))), 2),
      info_gain_bits = round(.gauss_H(S0) - .gauss_H(Sp), 2),
      corr_fossil_biogenic = round(lc(Ef, Eb), 2),
      sd_log_E_fossil = round(sd(log(pmax(Ef, 1e-9))), 3), sd_log_E_biogenic = round(sd(log(pmax(Eb, 1e-9))), 3),
      sd_log_E_total = round(sd(log(dr[, "E_total"])), 3), sd_fossil_share = round(sd(dr[, "fossil_share"]), 3),
      n_well = sum(regime == "well identified"), n_partial = sum(regime == "partially identified"),
      n_not = sum(regime == "not identified")),
    corr = C)
}

#' Long table of a posterior correlation matrix (upper triangle), for CSV output.
corr_long <- function(C) {
  ij <- which(upper.tri(C), arr.ind = TRUE)
  data.frame(a = rownames(C)[ij[, 1]], b = colnames(C)[ij[, 2]], corr = round(C[ij], 3))
}

#' Add "what each tracer added" columns: differences to the CH4-only row within the same group.
tracer_gain <- function(OV, group = c("rep"), base = "ch4") {
  do.call(rbind, lapply(split(OV, OV[group], drop = TRUE), function(s) {
    b <- s[s$config == base, ][1, ]
    s$dIG_vs_ch4_bits <- round(s$info_gain_bits - b$info_gain_bits, 2)
    s$dDOFS_vs_ch4 <- round(s$dofs - b$dofs, 2)
    s$sd_fossil_share_ratio_vs_ch4 <- round(s$sd_fossil_share / b$sd_fossil_share, 2)
    s
  }))
}

#' Posterior totals by reporting region (J$region_t_hr: [component, region] prior t/h).
#' Returns per region: total and fossil t/h (median, 90% interval), fossil share, and the posterior
#' fossil/biogenic correlation there. `truth` (named alpha) adds the true values for OSSEs.
region_posterior <- function(fit, K, RT, truth = NULL, fossil = FOSSIL_COMPONENTS) {
  if (is.null(RT)) return(NULL)
  A <- fit$draws("alpha", format = "draws_matrix")[, sprintf("alpha[%d]", seq_along(K)), drop = FALSE]
  RT <- RT[K, , drop = FALSE]; f <- K %in% fossil
  q <- function(v, p) if (all(is.na(v))) NA_real_ else unname(quantile(v, p, na.rm = TRUE))
  med <- function(v) q(v, 0.5)
  do.call(rbind, lapply(colnames(RT), function(r) {
    Ef <- as.vector(A[, f, drop = FALSE] %*% RT[f, r]); Eb <- as.vector(A[, !f, drop = FALSE] %*% RT[!f, r])
    Et <- Ef + Eb; sh <- ifelse(Et > 0, Ef / Et, NA)
    out <- data.frame(region = r, E_prior = round(sum(RT[, r]), 3),
      E_q05 = round(q(Et, .05), 3), E_q50 = round(med(Et), 3), E_q95 = round(q(Et, .95), 3),
      fossil_prior = round(sum(RT[f, r]), 3), Efos_q50 = round(med(Ef), 3),
      fossil_share_q05 = round(q(sh, .05), 3), fossil_share_q50 = round(med(sh), 3), fossil_share_q95 = round(q(sh, .95), 3),
      corr_fossil_biogenic = round(if (all(Ef > 0) && all(Eb > 0)) cor(log(Ef), log(Eb)) else NA, 2))
    if (!is.null(truth)) { tf <- sum(truth[K][f] * RT[f, r]); tt <- sum(truth[K] * RT[, r])
      out$E_true <- round(tt, 3); out$fossil_share_true <- if (tt > 0) round(tf / tt, 3) else NA
      out$covered_E <- isTRUE(tt >= q(Et, .05) & tt <= q(Et, .95))
      out$covered_share <- if (tt > 0) isTRUE((tf / tt) >= q(sh, .05) & (tf / tt) <= q(sh, .95)) else NA }
    out
  }))
}
