# 06b_identifiability_regimes.R — what controls whether ethane resolves the source split? ----------
# The OSSE (06) asks what THESE flights can identify. This asks WHY, by changing one thing at a
# time on the real flight tracks and footprints and refitting CH4-only and CH4+C2H6:
#
#   contrast   ethane:methane contrast between basin and urban gas (truth AND prior):
#                paper1  basin 0.080 / urban 0.110 (as config.R)
#                none    both 0.095 (fossil sources indistinguishable by ethane)
#                wide    basin 0.060 / urban 0.130
#   noise      ethane model-data noise s_C2H6 in the synthetic data: 0.3 ppb (as 06) or 1.0 ppb
#   sampling   all flights, or a random half of them (the same half in every cell)
#
# Every cell uses the same true alphas, transport factors, backgrounds and noise draws for a given
# replicate (common random numbers), so differences between cells come from the factor alone.
# Metrics (R/identifiability.R): DOFS, information gain, what ethane added (dIG), fossil/biogenic
# posterior correlation, sd of the fossil share, and the regime of each component.
#
#   Rscript scripts/06b_identifiability_regimes.R [prior=gra2pes_v2.0beta] [n_rep=2]
#   (3 x 2 x 2 cells x 2 tracer sets x n_rep fits; about 1-2 h at n_rep = 2)
#   Out: <RUN_DIR>/regimes_overall.csv, regimes_components.csv, figures/identifiability_regimes.png
# ----------------------------------------------------------------------------------------------
proj <- if (file.exists("config.R")) "." else ".."; source(file.path(proj, "config.R"))
source(file.path(proj, "R", "inversion.R")); source(file.path(proj, "R", "identifiability.R"))

a <- commandArgs(TRUE)
PRIOR <- if (length(a) >= 1) a[1] else "gra2pes_v2.0beta"
N_REP <- if (length(a) >= 2) as.integer(a[2]) else 2L

J0 <- readRDS(file.path(RUN_DIR, paste0("jacobian_", PRIOR, ".rds")))
S0 <- read.csv(file.path(INV_OUT, "obs_segments.csv"), stringsAsFactors = FALSE)
S0 <- S0[match(J0$ids, S0$id), ]; S0$leg_key <- paste(S0$flight, S0$leg_id)
seen <- colSums(J0$H) > 1e-6 * sum(J0$H)
J0$H <- J0$H[, seen, drop = FALSE]; K <- J0$components <- J0$components[seen]; J0$box_t_hr <- J0$box_t_hr[K]
st <- inv_settings(K); R_PRIOR0 <- R_PRIOR
if (is.null(J0$Hco)) J0$Hco <- rep(0, nrow(S0))

CONTRAST <- list(paper1 = c(og_basin = 0.080, og_urban = 0.110, postmeter = 0.110),
                 none   = c(og_basin = 0.095, og_urban = 0.095, postmeter = 0.095),
                 wide   = c(og_basin = 0.060, og_urban = 0.130, postmeter = 0.130))
NOISE <- c(low = 0.3, high = 1.0)
set.seed(99); fl_all <- unique(S0$flight); HALF <- sample(fl_all, ceiling(length(fl_all) / 2))
SAMPLING <- list(all = fl_all, half = HALF)
CONFIGS <- list(ch4 = character(0), ch4_c2h6 = "c2h6")
mod <- inv_model()

# replicate-level random draws, shared by all cells (common random numbers)
draw_rep <- function(rep) {
  set.seed(1000 + rep); nf <- length(fl_all); nl <- length(unique(S0$leg_key))
  list(alpha = setNames(exp(rnorm(length(K), 0, 0.7)), K), z_r = setNames(rnorm(length(K)), K),
       tau = setNames(exp(rnorm(nf, 0, 0.3)), fl_all), gamma = exp(rnorm(1, 0, 0.1)),
       b_ch4 = setNames(1960 + rnorm(nf, 0, 10), fl_all), b_c2 = setNames(1.5 + rnorm(nf, 0, 0.3), fl_all),
       b_co = setNames(100 + rnorm(nf, 0, 8), fl_all),
       u = matrix(rnorm(3 * nl), nl, 3, dimnames = list(unique(S0$leg_key), NULL)),
       e = matrix(rnorm(3 * nrow(S0)), nrow(S0), 3))
}
rprior_for <- function(cn) { rp <- R_PRIOR0; for (k in names(CONTRAST[[cn]])) rp[[k]]["mean"] <- CONTRAST[[cn]][[k]]; rp }

OV <- list(); CP <- list(); t0 <- Sys.time()
for (rep in seq_len(N_REP)) {
  D <- draw_rep(rep)
  for (cn in names(CONTRAST)) for (nz in names(NOISE)) for (sm in names(SAMPLING)) {
    R_PRIOR <- rprior_for(cn)                               # make_stan_data reads the global
    keep <- S0$flight %in% SAMPLING[[sm]]; S <- S0[keep, ]
    J <- J0; J$H <- J0$H[keep, , drop = FALSE]; J$Hco <- J0$Hco[keep]; J$ids <- J0$ids[keep]
    rm_ <- sapply(K, function(k) if (k %in% names(R_PRIOR)) R_PRIOR[[k]] else c(mean = 0, sd = 0))
    r <- pmax(0, rm_["mean", ] + rm_["sd", ] * D$z_r[K])
    tau <- D$tau[S$flight]; ee <- D$e[keep, , drop = FALSE]; uu <- D$u[S$leg_key, , drop = FALSE]
    e_ch4 <- tau * as.vector(J$H %*% D$alpha); e_c2 <- tau * as.vector(J$H %*% (D$alpha * r))
    e_co <- D$gamma * tau * ifelse(is.finite(J$Hco), J$Hco, 0)
    sdv <- function(e, s) sqrt(s^2 + (0.3 * e)^2)
    obs <- data.frame(flight = S$flight, leg_key = S$leg_key,
      ch4 = D$b_ch4[S$flight] + e_ch4 + 3 * uu[, 1] + sdv(e_ch4, 5) * ee[, 1],
      c2h6 = D$b_c2[S$flight] + e_c2 + 0.2 * uu[, 2] + sdv(e_c2, NOISE[[nz]]) * ee[, 2],
      co = D$b_co[S$flight] + e_co + 3 * uu[, 3] + sdv(e_co, 5) * ee[, 3])
    fos_true <- sum((D$alpha * J$box_t_hr)[K %in% FOSSIL_COMPONENTS]) / sum(D$alpha * J$box_t_hr)
    for (cf in names(CONFIGS)) {
      fit <- fit_inversion(mod, make_stan_data(obs, J, use = CONFIGS[[cf]], settings = st), seed = rep)
      idf <- identifiability(fit, K, st$la_sd)
      fs <- as.numeric(fit$draws("fossil_share", format = "draws_matrix"))
      cell <- data.frame(rep = rep, contrast = cn, noise = nz, sampling = sm, config = cf, n_obs = nrow(obs))
      OV[[length(OV) + 1]] <- cbind(cell, idf$overall, fossil_true = round(fos_true, 3),
        fossil_q50 = round(median(fs), 3), fossil_err = round(median(fs) - fos_true, 3),
        divergent = sum(fit$diagnostic_summary(quiet = TRUE)$num_divergent))
      CP[[length(CP) + 1]] <- cbind(cell, idf$component)
      cat(sprintf("[%s] rep %d %-6s noise %-4s %-4s %-9s DOFS %4.1f  IG %5.1f bits  corr(F,B) %5.2f  sd(fossil) %.3f  (%.0f min)\n",
        format(Sys.time(), "%H:%M"), rep, cn, nz, sm, cf, idf$overall$dofs, idf$overall$info_gain_bits,
        idf$overall$corr_fossil_biogenic, idf$overall$sd_fossil_share, as.numeric(difftime(Sys.time(), t0, units = "mins"))))
    }
  }
}
R_PRIOR <- R_PRIOR0
OVd <- do.call(rbind, OV); CPd <- do.call(rbind, CP)
OVd$cell_id <- paste(OVd$rep, OVd$contrast, OVd$noise, OVd$sampling)
OVd <- tracer_gain(OVd, group = "cell_id"); OVd$cell_id <- NULL
write.csv(OVd, inv_file("regimes_overall.csv"), row.names = FALSE)
write.csv(CPd, inv_file("regimes_components.csv"), row.names = FALSE)

# ---- summary: what ethane adds, by factor --------------------------------------------------------
E2 <- OVd[OVd$config == "ch4_c2h6", ]
m <- function(v) round(mean(v, na.rm = TRUE), 2)
tab <- aggregate(cbind(dIG_vs_ch4_bits, dDOFS_vs_ch4, sd_fossil_share_ratio_vs_ch4, corr_fossil_biogenic) ~ contrast + noise + sampling,
                 data = E2, FUN = m)
cat("\nwhat ethane adds (CH4+C2H6 vs CH4 alone), mean over replicates:\n"); print(tab, row.names = FALSE)
for (f in c("contrast", "noise", "sampling")) {
  cat(sprintf("\n  by %s:\n", f))
  print(setNames(aggregate(E2[c("dIG_vs_ch4_bits", "sd_fossil_share_ratio_vs_ch4")], by = list(E2[[f]]), FUN = m),
                 c(f, "dIG_vs_ch4_bits", "sd_fossil_share_ratio_vs_ch4")), row.names = FALSE)
}
RG <- with(CPd[CPd$contrast == "paper1" & CPd$noise == "low" & CPd$sampling == "all", ],
           tapply(regime, list(component, config), function(v) names(sort(table(v), decreasing = TRUE))[1]))
cat("\nregimes at the Paper-1 contrast, low noise, all flights:\n"); print(RG[K, , drop = FALSE])

# ---- figure ------------------------------------------------------------------------------------
png(file.path(RUN_DIR, "figures", paste0("identifiability_regimes_", INV_MODEL, ".png")), 1800, 700, res = 150)
par(mfrow = c(1, 2), mar = c(6, 4.5, 2.5, 1))
cells <- with(tab, paste(contrast, noise, sampling, sep = "/"))
barplot(tab$dIG_vs_ch4_bits, names.arg = cells, las = 2, col = "#4575b4", cex.names = 0.7,
        ylab = "information added by ethane (bits)")
title("(a) information gain from ethane, by cell", cex.main = 0.9)
barplot(tab$sd_fossil_share_ratio_vs_ch4, names.arg = cells, las = 2, col = "#d73027", cex.names = 0.7,
        ylab = "sd(fossil share): CH4+C2H6 / CH4", ylim = c(0, max(1.1, tab$sd_fossil_share_ratio_vs_ch4, na.rm = TRUE)))
abline(h = 1, lty = 3); title("(b) fossil-share uncertainty with ethane, relative to CH4 alone", cex.main = 0.9)
invisible(dev.off())
cat("\nwrote regimes_overall.csv, regimes_components.csv, figures/identifiability_regimes.png to", RUN_DIR, "\n")
