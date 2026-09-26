# inversion.R — assemble Stan data, fit with cmdstanr, summarise ---------------------------------
suppressPackageStartupMessages(library(cmdstanr))
invisible(use_local_cmdstan())

FOSSIL_COMPONENTS <- c("og_basin", "og_urban", "postmeter")

#' Prior settings (log-alpha sd, background and error scales). Edit here, not in the model.
#' Two versions (env METHANE_INV_MODEL, default v2):
#'   v1  the 24 Sep model: Gaussian errors, background sd 20 ppb around each flight's 5th percentile,
#'       leg offsets up to ~10 ppb, phi prior scale 0.5. On the real data (26 Sep) this let the
#'       backgrounds (+8 ppb median, +16 ppb on urban legs), leg offsets (|u| ~12 ppb urban) and phi
#'       (~0.95) absorb about half of the observed enhancement, so the alphas came out low.
#'   v2  Student-t errors (nu = 4) so unresolved plumes are outliers rather than evidence for a
#'       large phi; background sd 5 ppb; leg-offset scale 3 ppb (CH4), 0.2 (C2H6), 3 (CO); phi scale 0.3.
INV_MODEL <- Sys.getenv("METHANE_INV_MODEL", "v2")
inv_settings <- function(K_names, model = INV_MODEL) {
  base <- list(
    la_sd = setNames(rep(1.0, length(K_names)), K_names),     # alpha within ~x/2.7 at 1 sd
    b_ch4_sd = 20, b_c2h6_sd = 1.0, b_co_sd = 15,              # ppb, around the flight 5th percentile
    tau_sd = 0.5, gamma_sd = 0.15,
    s_ch4_scale = 20, s_c2h6_scale = 1.5, s_co_scale = 15, phi_scale = 0.5,
    leg_ch4_scale = 10, leg_c2h6_scale = 0.5, leg_co_scale = 8, nu = 1e6)
  if (model == "v2") base[c("b_ch4_sd", "b_c2h6_sd", "b_co_sd", "leg_ch4_scale", "leg_c2h6_scale",
                            "leg_co_scale", "phi_scale", "nu")] <- list(5, 0.3, 5, 3, 0.2, 3, 0.3, 4)
  else if (model != "v1") stop("METHANE_INV_MODEL must be v1 or v2")
  base
}
#' Output-file tag so v1 and v2 results sit side by side: posterior_summary_v2.csv etc.
inv_file <- function(name) { ext <- sub(".*(\\.[a-z]+)$", "\\1", name)
  file.path(RUN_DIR, paste0(sub("\\.[a-z]+$", "", name), "_", INV_MODEL, ext)) }

#' Build the Stan data list.
#' obs: data.frame with flight, leg_key, ch4, c2h6, co (NA allowed for c2h6/co);
#' J: jacobian list (H, Hco, components, box_t_hr); use: c("c2h6", "co") subsets.
#' bg_q: quantile of each flight's observations used as the background prior centre (default the 5th
#' percentile). b_ch4_mu: optional per-flight CH4 background centres (ppb, in levels(factor(obs$flight))
#' order) that override it; used by 07b_background_sensitivity.R.
make_stan_data <- function(obs, J, use = c("c2h6", "co"), settings = inv_settings(J$components),
                           bg_q = 0.05, b_ch4_mu = NULL) {
  K <- J$components; fl <- factor(obs$flight); lg <- factor(obs$leg_key)
  rp <- do.call(rbind, lapply(K, function(k) if (k %in% names(R_PRIOR)) R_PRIOR[[k]] else c(mean = 0, sd = 0)))
  i_c2 <- if ("c2h6" %in% use) which(is.finite(obs$c2h6)) else integer(0)
  i_co <- if ("co" %in% use && isTRUE(J$has_co)) which(is.finite(obs$co) & is.finite(J$Hco)) else integer(0)
  p05 <- function(v) tapply(v, fl, function(x) as.numeric(quantile(x, bg_q, na.rm = TRUE)))
  b_c2 <- p05(obs$c2h6); b_c2[!is.finite(b_c2)] <- 1.5
  b_co <- p05(obs$co);   b_co[!is.finite(b_co)] <- 100
  Hco <- J$Hco; Hco[!is.finite(Hco)] <- 0
  list(N = nrow(obs), K = length(K), F = nlevels(fl), L = nlevels(lg),
       flight = as.integer(fl), leg = as.integer(lg), H = unname(J$H), Hco = Hco, y_ch4 = obs$ch4,
       N_c2h6 = length(i_c2), idx_c2h6 = i_c2, y_c2h6 = obs$c2h6[i_c2],
       N_co = length(i_co), idx_co = i_co, y_co = obs$co[i_co],
       la_mu = rep(0, length(K)), la_sd = unname(settings$la_sd[K]),
       r_mu = unname(rp[, "mean"]), r_sd = unname(rp[, "sd"]),
       K_rfree = sum(rp[, "sd"] > 0), rfree = which(rp[, "sd"] > 0),
       b_ch4_mu = if (is.null(b_ch4_mu)) as.numeric(p05(obs$ch4)) else as.numeric(b_ch4_mu), b_c2h6_mu = as.numeric(b_c2), b_co_mu = as.numeric(b_co),
       b_ch4_sd = settings$b_ch4_sd, b_c2h6_sd = settings$b_c2h6_sd, b_co_sd = settings$b_co_sd,
       tau_sd = settings$tau_sd, gamma_sd = settings$gamma_sd,
       s_ch4_scale = settings$s_ch4_scale, s_c2h6_scale = settings$s_c2h6_scale, s_co_scale = settings$s_co_scale,
       leg_ch4_scale = settings$leg_ch4_scale, leg_c2h6_scale = settings$leg_c2h6_scale, leg_co_scale = settings$leg_co_scale,
       phi_scale = settings$phi_scale, nu = settings$nu,
       E_box = pmax(0, unname(J$box_t_hr[K])), is_fossil = as.integer(K %in% FOSSIL_COMPONENTS))
}

inv_model <- function(stan_file = file.path(INV_ROOT, "stan", "inversion.stan")) cmdstan_model(stan_file)

# show_exceptions = FALSE hides Stan's early-warmup 'Location parameter is inf' rejections (exp(log alpha)
# overflowing on wild initial proposals; harmless, sampler diagnostics are still checked). init = 0.5
# starts the unconstrained parameters closer to 0 so there are fewer of them.
fit_inversion <- function(mod, sdata, chains = 4, iter = as.integer(Sys.getenv("METHANE_INV_ITER", "1000")), seed = 1, ...)
  mod$sample(data = sdata, chains = chains, parallel_chains = min(chains, N_CORES), iter_warmup = iter,
             iter_sampling = iter, seed = seed, refresh = 0, show_messages = FALSE, show_exceptions = FALSE,
             init = 0.5, adapt_delta = 0.95, max_treedepth = 12, ...)

#' Per-component posterior summary of alpha and E_post.
summarise_alpha <- function(fit, K) {
  dr <- fit$draws(c("alpha", "E_post", "r"), format = "draws_matrix")
  q <- function(v) c(q05 = unname(quantile(v, 0.05)), q50 = median(v), q95 = unname(quantile(v, 0.95)))
  la <- log(dr[, paste0("alpha[", seq_along(K), "]"), drop = FALSE])
  data.frame(component = K,
    t(apply(dr[, paste0("alpha[", seq_along(K), "]"), drop = FALSE], 2, q)),
    sd_log_alpha = apply(la, 2, sd),
    E_q50 = apply(dr[, paste0("E_post[", seq_along(K), "]"), drop = FALSE], 2, median),
    r_q50 = apply(dr[, paste0("r[", seq_along(K), "]"), drop = FALSE], 2, median), row.names = NULL)
}
