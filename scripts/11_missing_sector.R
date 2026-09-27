# 11_missing_sector.R — two checks the reviewers asked for --------------------------------------------------
#   (A) three-way split of every posterior total: fossil (og_basin + og_urban + postmeter) / biogenic (waste,
#       landfills, wastewater, agriculture) / "other" (the unclassified remainder, ethane ratio prior 0.02 +- 0.02),
#       so that "biogenic share" is not silently 1 - fossil share.       -> results/table8_three_way_split.csv
#   (B) the missing-sector test: GRA2PES v1.1 with a waste sector added (METHANE_V11_WASTE=epa|v2 runs of
#       04/05/07, in runs/<transport>/v11waste_*/) against v1.1 as published. -> results/table7_missing_sector.csv
#   Rscript scripts/11_missing_sector.R      (env METHANE_MAIN_ZI=0.8)
proj <- if (file.exists("config.R")) "." else ".."; source(file.path(proj, "config.R"))
suppressPackageStartupMessages(library(cmdstanr)); invisible(use_local_cmdstan())   # to read the saved fits
MAIN_ZI <- as.numeric(Sys.getenv("METHANE_MAIN_ZI", "0.8")); MODEL <- Sys.getenv("METHANE_INV_MODEL", "v2"); CFG <- "ch4_c2h6"
OUT <- file.path(INV_OUT, "results"); dir.create(OUT, FALSE, TRUE)
runs_root <- file.path(INV_OUT, "runs")
main_dir <- file.path(runs_root, if (MAIN_ZI == 1) "np1000_h10" else sprintf("np1000_h10_zi%.2f", MAIN_ZI))
rd <- function(...) { f <- file.path(...); if (file.exists(f)) read.csv(f, stringsAsFactors = FALSE) else NULL }
PRIOR_LAB <- c(epa_ghgi_2020 = "EPA GHGI 2020", gra2pes_v1.1 = "GRA2PES v1.1", gra2pes_v2.0beta = "GRA2PES v2.0b")
FOSSIL_COMPONENTS <- c("og_basin", "og_urban", "postmeter")
REG <- c(paper1_box = "Metro box (urban)", obs_box = "Observed domain", djb = "DJ Basin")
BIO <- c("dads", "tower_road", "metro_complex", "metro_wwtp", "waste", "ag")

## (A) three-way split, computed WITHIN each posterior draw (alpha draws x prior region totals) --------------------
# so that totals, shares and their intervals are summaries of the same draws as Table 1, not sums of medians.
three_way <- function(fit_rds, prior_region_csv, prior_name, tag = "") {
  RT <- rd(prior_region_csv); if (is.null(RT) || !file.exists(fit_rds)) return(NULL)
  fit <- readRDS(fit_rds); A <- fit$draws("alpha", format = "draws_matrix"); Rr <- fit$draws("r", format = "draws_matrix")
  K <- ncol(A)
  J <- readRDS(sub("fit_.*$", paste0("jacobian_", prior_name, ".rds"), fit_rds))
  seen <- J$components[colSums(J$H) > 1e-6 * sum(J$H)]     # 07 drops unseen columns before fitting; same rule, same order
  if (length(seen) != K) stop("component mismatch for ", fit_rds)
  r <- RT[RT$prior == prior_name, ]
  do.call(rbind, lapply(names(REG), function(rg) {
    w <- r[[rg]][match(seen, r$component)]; w[is.na(w)] <- 0
    E <- A * matrix(w, nrow(A), K, byrow = TRUE)
    Et <- rowSums(E); Ef <- rowSums(E[, seen %in% FOSSIL_COMPONENTS, drop = FALSE]); Eb <- rowSums(E[, seen %in% BIO, drop = FALSE]); Eo <- if ("other" %in% seen) E[, which(seen == "other")] else 0
    q <- function(v, p) unname(quantile(v, p)); f1 <- function(v) sprintf("%.2f [%.2f-%.2f]", median(v), q(v, .05), q(v, .95))
    ro <- if ("other" %in% seen) median(Rr[, which(seen == "other")]) else NA
    data.frame(variant = tag, prior = PRIOR_LAB[prior_name], region = REG[rg],
      total_t_hr = sprintf("%.1f [%.1f-%.1f]", median(Et), q(Et, .05), q(Et, .95)), fossil_t_hr = round(median(Ef), 1), biogenic_t_hr = round(median(Eb), 1), other_t_hr = round(median(Eo), 1),
      fossil_share = f1(Ef / Et), biogenic_share = f1(Eb / Et), other_share = f1(Eo / Et),
      fossil_share_if_other_fossil = f1((Ef + Eo) / Et), other_ethane_ratio_post = round(ro, 3))
  }))
}
fit_of <- function(dir, p) file.path(dir, sprintf("fit_%s_%s_%s.rds", p, CFG, MODEL))
T8 <- do.call(rbind, lapply(names(PRIOR_LAB), function(p) three_way(fit_of(main_dir, p), file.path(INV_OUT, "prior_region_totals.csv"), p, "main")))
if (!is.null(T8)) { write.csv(T8, file.path(OUT, "table8_three_way_split.csv"), row.names = FALSE)
  cat("\nThree-way split (fossil / biogenic / other), within-draw, at ZISCALE", MAIN_ZI, ":\n"); print(T8[, -1], row.names = FALSE) }

## (B) missing-sector test --------------------------------------------------------------------------------------
rows <- list()
base <- rd(main_dir, paste0("posterior_region_totals_", MODEL, ".csv")); baseA <- rd(main_dir, paste0("posterior_alpha_", MODEL, ".csv"))
if (!is.null(base)) { b <- base[base$prior == "gra2pes_v1.1" & base$config == CFG & base$region == "paper1_box", ]
  ba <- baseA[baseA$prior == "gra2pes_v1.1" & baseA$config == CFG, ]
  rows[["v1.1 as published"]] <- data.frame(variant = "v1.1 as published", prior_box_t_hr = round(b$E_prior, 2), prior_fossil_share = round(b$fossil_prior / b$E_prior, 2),
    post_box_t_hr = sprintf("%.1f [%.1f-%.1f]", b$E_q50, b$E_q05, b$E_q95), post_fossil_share = sprintf("%.2f [%.2f-%.2f]", b$fossil_share_q50, b$fossil_share_q05, b$fossil_share_q95),
    alpha_postmeter = round(ba$q50[ba$component == "postmeter"], 1), alpha_waste_components = NA, domain_fossil_share = round(base$fossil_share_q50[base$prior == "gra2pes_v1.1" & base$config == CFG & base$region == "obs_box"], 2)) }
for (v in c("epa", "v2")) {
  d <- file.path(main_dir, paste0("v11waste_", v)); x <- rd(d, paste0("posterior_region_totals_", MODEL, ".csv")); xa <- rd(d, paste0("posterior_alpha_", MODEL, ".csv"))
  if (is.null(x)) next
  b <- x[x$prior == "gra2pes_v1.1" & x$config == CFG & x$region == "paper1_box", ]; ba <- xa[xa$prior == "gra2pes_v1.1" & xa$config == CFG, ]
  wa <- ba[ba$component %in% c("dads", "tower_road", "metro_complex", "waste"), ]
  rows[[v]] <- data.frame(variant = paste0("v1.1 + waste (", if (v == "epa") "EPA magnitude, v2.0b pattern" else "v2.0b magnitude and pattern", ")"),
    prior_box_t_hr = round(b$E_prior, 2), prior_fossil_share = round(b$fossil_prior / b$E_prior, 2),
    post_box_t_hr = sprintf("%.1f [%.1f-%.1f]", b$E_q50, b$E_q05, b$E_q95), post_fossil_share = sprintf("%.2f [%.2f-%.2f]", b$fossil_share_q50, b$fossil_share_q05, b$fossil_share_q95),
    alpha_postmeter = round(ba$q50[ba$component == "postmeter"], 1), alpha_waste_components = paste(sprintf("%s %.2f", wa$component, wa$q50), collapse = "; "),
    domain_fossil_share = round(x$fossil_share_q50[x$prior == "gra2pes_v1.1" & x$config == CFG & x$region == "obs_box"], 2))
  T8v <- three_way(fit_of(d, "gra2pes_v1.1"), file.path(INV_OUT, "priors", paste0("v11waste_", v), "prior_region_totals.csv"), "gra2pes_v1.1", paste0("v11waste_", v))
  if (!is.null(T8v)) write.csv(T8v, file.path(OUT, sprintf("table8_three_way_split_v11waste_%s.csv", v)), row.names = FALSE)
}
if (length(rows)) { T7 <- do.call(rbind, rows); write.csv(T7, file.path(OUT, "table7_missing_sector.csv"), row.names = FALSE)
  cat("\nMissing-sector test (metro box, CH4+C2H6, ZISCALE", MAIN_ZI, "):\n"); print(T7, row.names = FALSE)
  if (nrow(T7) == 1) cat("  (no v11waste_* runs found: run 04/05/07 with METHANE_V11_WASTE=epa and =v2, see README)\n") }
cat("wrote table7/table8 to", OUT, "\n")
