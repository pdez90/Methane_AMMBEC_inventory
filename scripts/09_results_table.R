# 09_results_table.R — assemble the Paper 2 main results from the runs already on disk ------------------
# Reads every completed transport run under <INV_OUT>/runs/np1000_h10* (07 outputs), the background
# sensitivity (07b, ZISCALE 0.8) and the aircraft BLH summary (05c), and writes the tables the paper
# needs. Nothing is refit here.
#
#   T1  main result: ZISCALE 0.8, CH4+C2H6, by prior x region (prior, posterior [90%], fossil share)
#   T2  transport sweep: ZISCALE x prior (CH4+C2H6): box / obs_box / DJB totals, fossil share, fit rho,
#       urban fitted-signal share (the diagnostic that shows 0.37 and 1.2 failing)
#   T3  tracer sets at ZISCALE 0.8: config x prior: box total, fossil share, DOFS, IG, dIG vs CH4-only
#   T4  uncertainty budget for each region total and fossil share: spread across priors, across transport
#       (0.8 vs 1.0 and the full sweep), across background cases (07b), and the posterior 90% width
#   T5  headline ensemble: priors x transport {0.8, 1.0} x background {base, drift}: median and range
#
#   Rscript scripts/09_results_table.R            (main transport 0.8, range 0.8-1.0; override with
#   METHANE_MAIN_ZI=0.8 METHANE_RANGE_ZI=0.8,1.0)
#   Out: <INV_OUT>/results/table*.csv, results_summary.md
# ----------------------------------------------------------------------------------------------------
proj <- if (file.exists("config.R")) "." else ".."; source(file.path(proj, "config.R"))
MAIN_ZI  <- as.numeric(Sys.getenv("METHANE_MAIN_ZI", "0.8"))
RANGE_ZI <- as.numeric(strsplit(Sys.getenv("METHANE_RANGE_ZI", "0.8,1.0"), ",")[[1]])
MODEL <- Sys.getenv("METHANE_INV_MODEL", "v2"); CFG_MAIN <- "ch4_c2h6"
OUT <- file.path(INV_OUT, "results"); dir.create(OUT, FALSE, TRUE)
runs_root <- file.path(INV_OUT, "runs")
zi_of <- function(d) { m <- regmatches(d, regexpr("_zi([0-9.]+)$", d)); if (length(m)) as.numeric(sub("_zi", "", m)) else 1 }
rd <- function(...) { f <- file.path(...); if (file.exists(f)) read.csv(f, stringsAsFactors = FALSE) else NULL }
f1 <- function(x, d = 1) formatC(x, format = "f", digits = d)
ci <- function(q50, q05, q95, d = 1) sprintf("%s [%s-%s]", f1(q50, d), f1(q05, d), f1(q95, d))
lsd <- function(v) { v <- v[is.finite(v) & v > 0]; if (length(v) < 2) NA_real_ else sd(log(v)) }
PRIOR_LAB <- c(epa_ghgi_2020 = "EPA GHGI 2020", gra2pes_v1.1 = "GRA2PES v1.1", gra2pes_v2.0beta = "GRA2PES v2.0b")
REG_LAB <- c(paper1_box = "Paper 1 box (urban)", obs_box = "Observed domain", djb = "DJ Basin")

## ---- gather every completed transport run ----------------------------------------------------------
dirs <- list.dirs(runs_root, recursive = FALSE, full.names = FALSE)
dirs <- dirs[grepl("^np1000_h10", dirs)]
RT <- list(); ID <- list(); FS <- list(); DC <- list()
for (d in dirs) {
  r <- rd(runs_root, d, paste0("posterior_region_totals_", MODEL, ".csv")); if (is.null(r)) next
  r$ziscale <- zi_of(d); r$transport <- d; RT[[d]] <- r
  i <- rd(runs_root, d, paste0("posterior_identifiability_overall_", MODEL, ".csv")); if (!is.null(i)) { i$ziscale <- zi_of(d); ID[[d]] <- i }
  s <- rd(runs_root, d, "forward_summary.csv"); if (!is.null(s)) { s$ziscale <- zi_of(d); FS[[d]] <- s }
  dc <- rd(runs_root, d, paste0("posterior_decomposition_", MODEL, ".csv")); if (!is.null(dc)) { dc$ziscale <- zi_of(d); DC[[d]] <- dc }
}
if (!length(RT)) stop("no posterior_region_totals_", MODEL, ".csv under ", runs_root)
RT <- do.call(rbind, RT); ID <- do.call(rbind, ID); FS <- do.call(rbind, FS); DC <- do.call(rbind, DC)
cat(sprintf("transport runs found: %s\n", paste(sort(unique(RT$ziscale)), collapse = " ")))
if (!MAIN_ZI %in% RT$ziscale) stop("main ZISCALE ", MAIN_ZI, " has no run")
main_dir <- dirs[sapply(dirs, zi_of) == MAIN_ZI][1]
BG <- rd(runs_root, main_dir, paste0("background_sensitivity_", MODEL, ".csv"))
blh_txt <- file.path(runs_root, "np1000_h10", "blh_aircraft_summary.txt")
blh_zi <- if (file.exists(blh_txt)) { l <- grep("implied ZISCALE", readLines(blh_txt), value = TRUE)
  if (length(l)) as.numeric(sub(".*uncensored\\): *([0-9.]+).*", "\\1", l[1])) else NA } else NA

## ---- T1: main result ---------------------------------------------------------------------------------
m <- RT[RT$ziscale == MAIN_ZI & RT$config == CFG_MAIN, ]
T1 <- data.frame(prior = PRIOR_LAB[m$prior], region = REG_LAB[m$region],
  prior_t_hr = f1(m$E_prior), posterior_t_hr = ci(m$E_q50, m$E_q05, m$E_q95),
  ratio_post_prior = f1(m$E_q50 / m$E_prior, 2),
  prior_fossil_share = f1(m$fossil_prior / m$E_prior, 2),
  posterior_fossil_share = ci(m$fossil_share_q50, m$fossil_share_q05, m$fossil_share_q95, 2),
  fossil_t_hr = f1(m$Efos_q50), row.names = NULL)
T1 <- T1[order(match(m$region, names(REG_LAB)), match(m$prior, names(PRIOR_LAB))), ]
write.csv(T1, file.path(OUT, "table1_main_results.csv"), row.names = FALSE)

## ---- T2: transport sweep -----------------------------------------------------------------------------
sw <- RT[RT$config == CFG_MAIN, ]
w <- function(reg, col) { x <- sw[sw$region == reg, ]; tapply(x[[col]], list(x$ziscale, x$prior), function(v) v[1]) }
box <- w("paper1_box", "E_q50"); obs <- w("obs_box", "E_q50"); djb <- w("djb", "E_q50")
fbox <- w("paper1_box", "fossil_share_q50"); fobs <- w("obs_box", "fossil_share_q50"); fdjb <- w("djb", "fossil_share_q50")
rho <- if (!is.null(FS)) tapply(FS$rho_ch4, list(FS$ziscale, FS$prior), function(v) v[1]) else NULL
usig <- if (!is.null(DC)) { u <- DC[DC$config == CFG_MAIN & DC$region == "urban", ]; tapply(u$signal_share, list(u$ziscale, u$prior), function(v) v[1]) } else NULL
T2 <- do.call(rbind, lapply(rownames(box), function(z) data.frame(ziscale = as.numeric(z),
  box_EPA_v11_v2 = paste(f1(box[z, names(PRIOR_LAB)]), collapse = " / "),
  obs_box = paste(f1(obs[z, names(PRIOR_LAB)]), collapse = " / "),
  djb = paste(f1(djb[z, names(PRIOR_LAB)]), collapse = " / "),
  fossil_box = paste(f1(fbox[z, names(PRIOR_LAB)], 2), collapse = " / "),
  fossil_obs_box = paste(f1(fobs[z, names(PRIOR_LAB)], 2), collapse = " / "),
  fossil_djb = paste(f1(fdjb[z, names(PRIOR_LAB)], 2), collapse = " / "),
  rho_fit = if (!is.null(rho) && z %in% rownames(rho)) paste(f1(rho[z, names(PRIOR_LAB)], 2), collapse = " / ") else NA,
  urban_signal_share = if (!is.null(usig) && z %in% rownames(usig)) paste(f1(100 * usig[z, names(PRIOR_LAB)], 0), collapse = " / ") else NA)))
T2 <- T2[order(T2$ziscale), ]; write.csv(T2, file.path(OUT, "table2_transport_sweep.csv"), row.names = FALSE)

## ---- T3: tracer sets at the main transport ---------------------------------------------------------
im <- ID[ID$ziscale == MAIN_ZI, ]; rm3 <- RT[RT$ziscale == MAIN_ZI & RT$region == "paper1_box", ]
T3 <- merge(rm3[, c("prior", "config", "E_q50", "E_q05", "E_q95", "fossil_share_q50", "fossil_share_q05", "fossil_share_q95")],
            im[, c("prior", "config", "dofs", "info_gain_bits", "dIG_vs_ch4_bits", "sd_fossil_share", "n_well", "n_partial", "n_not")], by = c("prior", "config"))
T3 <- T3[order(match(T3$prior, names(PRIOR_LAB)), match(T3$config, c("ch4", "ch4_c2h6", "ch4_c2h6_co"))), ]
T3o <- data.frame(prior = PRIOR_LAB[T3$prior], tracers = T3$config, box_t_hr = ci(T3$E_q50, T3$E_q05, T3$E_q95),
  fossil_share = ci(T3$fossil_share_q50, T3$fossil_share_q05, T3$fossil_share_q95, 2), DOFS = f1(T3$dofs, 1),
  info_gain_bits = f1(T3$info_gain_bits, 1), dIG_vs_ch4 = f1(T3$dIG_vs_ch4_bits, 1), sd_fossil_share = f1(T3$sd_fossil_share, 3),
  regimes_well_partial_not = paste(T3$n_well, T3$n_partial, T3$n_not, sep = "/"), row.names = NULL)
write.csv(T3o, file.path(OUT, "table3_tracer_sets.csv"), row.names = FALSE)

## ---- T4: uncertainty budget ---------------------------------------------------------------------------
budget <- function(reg, col, lab) {
  x <- RT[RT$config == CFG_MAIN & RT$region == reg, ]
  at_main <- x[x$ziscale == MAIN_ZI, ]; in_range <- x[x$ziscale %in% RANGE_ZI, ]
  post_w <- median(log(at_main[[sub("q50", "q95", col)]] / at_main[[sub("q50", "q05", col)]])) / (2 * qnorm(0.95))  # 90% -> 1 sigma log
  prior_sd <- lsd(at_main[[col]])
  tr_range <- sapply(split(in_range, in_range$prior), function(s) lsd(s[[col]]))           # per prior, across range ZI
  tr_full  <- sapply(split(x, x$prior), function(s) lsd(s[[col]]))
  bgcol <- c(paper1_box = "E_q50", obs_box = "obs_box_q50", djb = "djb_q50")[reg]
  bg_sd <- if (!is.null(BG) && col == "E_q50" && bgcol %in% names(BG)) median(sapply(split(BG, BG$prior), function(s) lsd(s[[bgcol]]))) else NA
  bg_drift <- if (!is.null(BG) && col == "E_q50" && bgcol %in% names(BG)) median(sapply(split(BG, BG$prior), function(s) s[[bgcol]][s$case == "drift"] / s[[bgcol]][s$case == "base"])) else NA
  data.frame(quantity = lab, posterior_1sigma_log = round(post_w, 3), across_priors_log_sd = round(prior_sd, 3),
    across_transport_range_log_sd = round(median(tr_range, na.rm = TRUE), 3), across_transport_full_sweep_log_sd = round(median(tr_full, na.rm = TRUE), 3),
    across_background_log_sd = round(bg_sd, 3), drift_over_base = round(bg_drift, 2),
    dominant = c("posterior", "prior", "transport", "background")[which.max(c(post_w, prior_sd, median(tr_range, na.rm = TRUE), ifelse(is.na(bg_sd), 0, bg_sd)))])
}
T4 <- rbind(budget("paper1_box", "E_q50", "Paper 1 box total"), budget("obs_box", "E_q50", "Observed-domain total"), budget("djb", "E_q50", "DJ Basin total"))
# fossil shares: absolute spreads, not log
fb <- function(reg, lab) { x <- RT[RT$config == CFG_MAIN & RT$region == reg, ]; a <- x[x$ziscale == MAIN_ZI, ]; r <- x[x$ziscale %in% RANGE_ZI, ]
  bgc <- c(paper1_box = "fossil_q50", obs_box = NA, djb = "djb_fossil_share")[reg]
  data.frame(quantity = lab, posterior_1sigma = round(median((a$fossil_share_q95 - a$fossil_share_q05) / (2 * qnorm(.95))), 3),
    across_priors_range = sprintf("%.2f-%.2f", min(a$fossil_share_q50), max(a$fossil_share_q50)),
    across_transport_range_sd = round(median(sapply(split(r, r$prior), function(s) sd(s$fossil_share_q50))), 3),
    across_background_sd = if (!is.null(BG) && !is.na(bgc)) round(median(sapply(split(BG, BG$prior), function(s) sd(s[[bgc]]))), 3) else NA) }
T4f <- rbind(fb("paper1_box", "Paper 1 box fossil share"), fb("obs_box", "Observed-domain fossil share"), fb("djb", "DJ Basin fossil share"))
write.csv(T4, file.path(OUT, "table4_budget_totals.csv"), row.names = FALSE); write.csv(T4f, file.path(OUT, "table4_budget_fossil.csv"), row.names = FALSE)

## ---- T5: headline ensemble ----------------------------------------------------------------------------
# members: prior x ZISCALE in RANGE x background {base, drift}. Base = 07 run at that ZISCALE; drift only exists
# at MAIN_ZI (07b), applied as a ratio to the other transport(s).
ens <- RT[RT$config == CFG_MAIN & RT$ziscale %in% RANGE_ZI, ]
drift_ratio <- if (!is.null(BG)) sapply(split(BG, BG$prior), function(s) c(paper1_box = s$E_q50[s$case == "drift"] / s$E_q50[s$case == "base"],
  obs_box = s$obs_box_q50[s$case == "drift"] / s$obs_box_q50[s$case == "base"], djb = s$djb_q50[s$case == "drift"] / s$djb_q50[s$case == "base"])) else NULL
T5 <- do.call(rbind, lapply(names(REG_LAB), function(reg) { e <- ens[ens$region == reg, ]
  tot <- e$E_q50; if (!is.null(drift_ratio)) tot <- c(tot, e$E_q50 * drift_ratio[reg, e$prior])
  fs <- e$fossil_share_q50
  data.frame(region = REG_LAB[reg], n_members = length(tot), total_median = f1(median(tot)), total_range = sprintf("%s-%s", f1(min(tot)), f1(max(tot))),
    total_90pct_members = sprintf("%s-%s", f1(min(e$E_q05)), f1(max(e$E_q95))),
    fossil_share_median = f1(median(fs), 2), fossil_share_range = sprintf("%.2f-%.2f", min(fs), max(fs)),
    fossil_share_excl_v11 = { g <- e$fossil_share_q50[e$prior != "gra2pes_v1.1"]; sprintf("%.2f (%.2f-%.2f)", median(g), min(g), max(g)) },
    biogenic_share_median = f1(1 - median(fs), 2), row.names = NULL) }))
write.csv(T5, file.path(OUT, "table5_headline.csv"), row.names = FALSE)

## ---- T6: prior dependence by region (main transport, CH4+C2H6) ------------------------------------------
# Totals: ratio of the across-inventory spread of posterior medians (sd of log) to the same spread of
# the priors; 0 = data-determined, 1 = prior-determined. The second column compares the across-prior
# spread with the posterior 1-sigma width instead ("do the inventories still disagree by more than the
# posterior width?"): < 1 means the three posteriors overlap. Fossil shares use absolute ranges; when
# the priors happen to agree (domain: 0.38-0.42) the ratio is not informative and only the ranges are
# reported.
T6 <- do.call(rbind, lapply(names(REG_LAB), function(reg) {
  x <- RT[RT$config == CFG_MAIN & RT$ziscale == MAIN_ZI & RT$region == reg, ]
  lsd <- function(v) sd(log(pmax(v, 1e-9))); pw <- mean(log(x$E_q95 / x$E_q05) / (2 * qnorm(0.95)))
  pf <- x$fossil_prior / x$E_prior
  data.frame(region = REG_LAB[reg], prior_spread_log = round(lsd(x$E_prior), 3), post_spread_log = round(lsd(x$E_q50), 3),
    total_prior_dependence = round(lsd(x$E_q50) / lsd(x$E_prior), 2), post_spread_over_post_sigma = round(lsd(x$E_q50) / pw, 2),
    fossil_prior_range = sprintf("%.2f-%.2f", min(pf), max(pf)), fossil_post_range = sprintf("%.2f-%.2f", min(x$fossil_share_q50), max(x$fossil_share_q50)),
    fossil_prior_dependence = if (diff(range(pf)) < 0.1) NA else round(diff(range(x$fossil_share_q50)) / diff(range(pf)), 2))
}))
write.csv(T6, file.path(OUT, "table6_prior_dependence.csv"), row.names = FALSE)

## ---- print + markdown ---------------------------------------------------------------------------------
knitr_free <- function(d) { d[] <- lapply(d, as.character); c(paste("|", paste(names(d), collapse = " | "), "|"),
  paste("|", paste(rep("---", ncol(d)), collapse = " | "), "|"), apply(d, 1, function(r) paste("|", paste(r, collapse = " | "), "|"))) }
md <- c(sprintf("# Paper 2 results (%s, %s; main transport ZISCALE %.2f, range %s)", MODEL, CFG_MAIN, MAIN_ZI, paste(RANGE_ZI, collapse = "-")),
  if (is.finite(blh_zi)) sprintf("Aircraft-implied ZISCALE (05c, lower bound): %.2f", blh_zi) else "", "",
  sprintf("## T1 main result (ZISCALE %.2f, CH4+C2H6)", MAIN_ZI), knitr_free(T1), "",
  "## T2 transport sweep (values EPA / v1.1 / v2)", knitr_free(T2), "",
  "## T3 tracer sets at the main transport", knitr_free(T3o), "",
  "## T4 uncertainty budget: totals (log sd; 1 sigma)", knitr_free(T4), "", "## T4 uncertainty budget: fossil shares", knitr_free(T4f), "",
  "## T5 headline ensemble (priors x transport range x background base/drift)", knitr_free(T5), "",
  "## T6 prior dependence by region (main transport, CH4+C2H6)", knitr_free(T6))
writeLines(md, file.path(OUT, "results_summary.md"))
for (t in list(T1, T2, T3o, T4, T4f, T5, T6)) { print(t, row.names = FALSE); cat("\n") }
cat("wrote table1-6 csv and results_summary.md to ", OUT, "\n", sep = "")
