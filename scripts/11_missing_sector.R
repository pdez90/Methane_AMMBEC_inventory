# 11_missing_sector.R — two checks the reviewers asked for --------------------------------------------------
#   (A) three-way split of every posterior total: fossil (og_basin + og_urban + postmeter) / biogenic (waste,
#       landfills, wastewater, agriculture) / "other" (the unclassified remainder, ethane ratio prior 0.02 +- 0.02),
#       so that "biogenic share" is not silently 1 - fossil share.       -> results/table8_three_way_split.csv
#   (B) the missing-sector test: GRA2PES v1.1 with a waste sector added (METHANE_V11_WASTE=epa|v2 runs of
#       04/05/07, in runs/<transport>/v11waste_*/) against v1.1 as published. -> results/table7_missing_sector.csv
#   Rscript scripts/11_missing_sector.R      (env METHANE_MAIN_ZI=0.8)
proj <- if (file.exists("config.R")) "." else ".."; source(file.path(proj, "config.R"))
MAIN_ZI <- as.numeric(Sys.getenv("METHANE_MAIN_ZI", "0.8")); MODEL <- Sys.getenv("METHANE_INV_MODEL", "v2"); CFG <- "ch4_c2h6"
OUT <- file.path(INV_OUT, "results"); dir.create(OUT, FALSE, TRUE)
runs_root <- file.path(INV_OUT, "runs")
main_dir <- file.path(runs_root, if (MAIN_ZI == 1) "np1000_h10" else sprintf("np1000_h10_zi%.2f", MAIN_ZI))
rd <- function(...) { f <- file.path(...); if (file.exists(f)) read.csv(f, stringsAsFactors = FALSE) else NULL }
PRIOR_LAB <- c(epa_ghgi_2020 = "EPA GHGI 2020", gra2pes_v1.1 = "GRA2PES v1.1", gra2pes_v2.0beta = "GRA2PES v2.0b")
FOSSIL_COMPONENTS <- c("og_basin", "og_urban", "postmeter")
REG <- c(paper1_box = "Metro box (urban)", obs_box = "Observed domain", djb = "DJ Basin")
BIO <- c("dads", "tower_road", "metro_complex", "metro_wwtp", "waste", "ag")

## (A) three-way split: alpha (posterior median) x prior region total per component ----------------------------
three_way <- function(alpha_csv, prior_region_csv, tag = "") {
  A <- rd(alpha_csv); RT <- rd(prior_region_csv); if (is.null(A) || is.null(RT)) return(NULL)
  A <- A[A$config == CFG, ]
  do.call(rbind, lapply(unique(A$prior), function(p) do.call(rbind, lapply(names(REG), function(rg) {
    a <- A[A$prior == p, ]; r <- RT[RT$prior == p, ]; e <- a$q50 * r[[rg]][match(a$component, r$component)]
    e[is.na(e)] <- 0; tot <- sum(e); fos <- sum(e[a$component %in% c("og_basin", "og_urban", "postmeter")]); bio <- sum(e[a$component %in% BIO]); oth <- sum(e[a$component == "other"])
    ro <- a$r_q50[a$component == "other"]
    data.frame(variant = tag, prior = PRIOR_LAB[p], region = REG[rg], total_t_hr = round(tot, 1), fossil_t_hr = round(fos, 1), biogenic_t_hr = round(bio, 1), other_t_hr = round(oth, 1),
               fossil_share = round(fos / tot, 2), biogenic_share = round(bio / tot, 2), other_share = round(oth / tot, 2),
               fossil_share_if_other_fossil = round((fos + oth) / tot, 2), other_ethane_ratio_post = if (length(ro)) round(ro, 3) else NA)
  }))))
}
T8 <- three_way(file.path(main_dir, paste0("posterior_alpha_", MODEL, ".csv")), file.path(INV_OUT, "prior_region_totals.csv"), "main")
if (!is.null(T8)) { write.csv(T8, file.path(OUT, "table8_three_way_split.csv"), row.names = FALSE)
  cat("\nThree-way split (fossil / biogenic / other) at ZISCALE", MAIN_ZI, "-- posterior medians x prior region totals:\n"); print(T8[, -1], row.names = FALSE) }

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
  T8v <- three_way(file.path(d, paste0("posterior_alpha_", MODEL, ".csv")), file.path(INV_OUT, "priors", paste0("v11waste_", v), "prior_region_totals.csv"), paste0("v11waste_", v))
  if (!is.null(T8v)) write.csv(T8v, file.path(OUT, sprintf("table8_three_way_split_v11waste_%s.csv", v)), row.names = FALSE)
}
if (length(rows)) { T7 <- do.call(rbind, rows); write.csv(T7, file.path(OUT, "table7_missing_sector.csv"), row.names = FALSE)
  cat("\nMissing-sector test (metro box, CH4+C2H6, ZISCALE", MAIN_ZI, "):\n"); print(T7, row.names = FALSE)
  if (nrow(T7) == 1) cat("  (no v11waste_* runs found: run 04/05/07 with METHANE_V11_WASTE=epa and =v2, see README)\n") }
cat("wrote table7/table8 to", OUT, "\n")
