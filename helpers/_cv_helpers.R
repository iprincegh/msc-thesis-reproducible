# _cv_helpers.R
# Cross-validation helper utilities (10a/10b).
# Sourced by: 10a_temporal_cv.Rmd, 10b_spatial_cv.Rmd

# ── Constants (mirroring 09_fit_model) ──
RHO_TRIM_LOW       <- 0.001
RHO_TRIM_HIGH      <- 50
SNAP_QUEEN          <- 10
SNAP_TIERS          <- c(300, 500, 750, 1000, 1500)
SNAP_MAX            <- 2000
OD_CORR_THRESHOLD   <- 0.6
OD_MIN_ACTIVITY     <- 5.0
OD_MIN_DAYS         <- 10
OD_MIN_ROWS         <- 1000
OD_MAX_FRAG_SIZE    <- 10
OD_MAX_FRAG_DEGREE  <- 1.5
MAX_WORMHOLES       <- 10

MIN_CASES          <- 1L  # match Slater filter in 09_fit_model (MIN_WAVE_CASES = 1)
OUTLIER_QUANTILE   <- 0.999
EXCLUDE_COUNTIES   <- list(UT_wave2 = c(27, 24, 20))

UT_EXCLUDE_TEST <- c("Weber", "Sanpete", "Iron", "Carbon", "Cache",
                     "Duchesne", "Daggett", "Piute", "Rich", "Wayne")

# Per-wave best model (selected by `09_fit_model` via WAIC).
# Hard-coded canonical models used by CV routines.
BEST_MODEL_BY_WAVE <- list(
  NY_wave1 = "NLMob_ST_IV",
  NY_wave2 = "LinMob_ST_IV",
  NY_wave3 = "LinMob_ST_IV",
  UT_wave1 = "NLMob_ST_IV",
  UT_wave2 = "NLMob_ST_IV",
  UT_wave3 = "NLMob_ST_IV"
)

# Scoring functions ---------------------------------------------------------
# Compute point and distributional scoring metrics for counts.
# Primary function: `cv_scores()` returns RMSE/MAE/R2, CRPS, nominal coverages,
# log-score, DSS, calibration slope/intercept, and PIT W1.
cv_scores <- function(y, mu, sd_fit, nb_size = NULL,
                      family = c("auto", "nbinomial", "poisson", "gaussian")) {
  family <- match.arg(family)
  ok <- !is.na(y) & !is.na(mu) & is.finite(mu)
  y <- y[ok]; mu <- mu[ok]; sd_fit <- sd_fit[ok]
  n <- length(y)
  empty <- data.table(RMSE = NA_real_, MAE = NA_real_, R2 = NA_real_,
                       CRPS = NA_real_, Cov50 = NA_real_, Cov80 = NA_real_,
                       Cov95 = NA_real_, Cov95_gauss = NA_real_,
                       logS = NA_real_, DSS = NA_real_,
                       cal_intercept = NA_real_, cal_slope = NA_real_,
                       pit_W1 = NA_real_, N = 0L)
  if (n < 10) return(empty)

  cap <- max(2 * max(y), quantile(c(y, mu), 0.999))
  mu  <- pmin(mu, cap)
  mu_safe <- pmax(mu, 1e-10)

  rmse <- sqrt(mean((y - mu)^2))
  mae  <- mean(abs(y - mu))
  ss_res <- sum((y - mu)^2); ss_tot <- sum((y - mean(y))^2)
  r2 <- if (ss_tot > 0) 1 - ss_res / ss_tot else NA_real_

  has_nb <- !is.null(nb_size) && is.finite(nb_size) && nb_size > 0
  # Resolve family: explicit override wins; otherwise auto-detect from
  # nb_size presence + integer-count test on y. NY-w1/w2 (var/mean ~1.06)
  # fit as Poisson by fit_model auto-switch, so nb_size is NULL there and
  # the Poisson branch below is what we want.
  fam <- if (family != "auto") family
         else if (has_nb) "nbinomial"
         else if (length(y) > 0 && all(y >= 0) && all(y == round(y))) "poisson"
         else "gaussian"

  if (fam == "nbinomial" && has_nb) {
    # Closed-form NB CRPS (scoringRules) when available; vectorised, exact.
    # scoringRules::crps_nbinom requires the `hypergeo` package; if it's not
    # installed or the call errors for any reason, fall back to MC.
    mc_crps <- function() {
      set.seed(42)
      idx <- if (n > 50000) sample.int(n, 50000) else seq_len(n)
      mean(vapply(idx, function(i) {
        m <- mu_safe[i]
        s1 <- rnbinom(200, size = nb_size, mu = m)
        s2 <- rnbinom(200, size = nb_size, mu = m)
        mean(abs(s1 - y[i])) - 0.5 * mean(abs(s1 - s2))
      }, numeric(1)))
    }
    crps <- if (requireNamespace("scoringRules", quietly = TRUE) &&
                requireNamespace("hypergeo", quietly = TRUE)) {
      tryCatch(
        mean(scoringRules::crps_nbinom(y, mu = mu_safe, size = nb_size)),
        error = function(e) mc_crps())
    } else {
      mc_crps()
    }
    # Nominal-coverage at 50%, 80%, 95% under the NB predictive.
    cov_pair <- function(p) {
      a <- (1 - p) / 2
      lo <- qnbinom(a,     size = nb_size, mu = mu_safe)
      hi <- qnbinom(1 - a, size = nb_size, mu = mu_safe)
      mean(y >= lo & y <= hi)
    }
    cov50 <- cov_pair(0.50); cov80 <- cov_pair(0.80); cov95 <- cov_pair(0.95)
    logS  <- -mean(dnbinom(y, size = nb_size, mu = mu_safe, log = TRUE))
    pred_sd <- sqrt(mu_safe + mu_safe^2 / nb_size)
    DSS  <- mean(((y - mu_safe) / pred_sd)^2 + 2 * log(pred_sd))
    # Randomized PIT (Czado et al. 2009) -> Wasserstein-1 to Uniform(0,1).
    rpit <- randomized_pit_nb(y, mu_safe, nb_size)
    rpit <- rpit[is.finite(rpit) & rpit >= 0 & rpit <= 1]
    pit_W1 <- if (length(rpit) >= 50) {
      ord <- sort(rpit); m <- length(ord)
      mean(abs(ord - (seq_len(m) - 0.5) / m))
    } else NA_real_
  } else if (fam == "poisson") {
    # Poisson predictive: variance == mean. Used by NY-w1/w2 (var/mean < 1.1)
    # where fit_model's auto-switch picks Poisson and nb_size is NULL.
    crps <- if (requireNamespace("scoringRules", quietly = TRUE)) {
      mean(scoringRules::crps_pois(y, lambda = mu_safe))
    } else {
      # Fallback: MC estimator on subsample.
      set.seed(42)
      idx <- if (n > 50000) sample.int(n, 50000) else seq_len(n)
      mean(vapply(idx, function(i) {
        m <- mu_safe[i]
        s1 <- rpois(200, lambda = m); s2 <- rpois(200, lambda = m)
        mean(abs(s1 - y[i])) - 0.5 * mean(abs(s1 - s2))
      }, numeric(1)))
    }
    cov_pair_p <- function(p) {
      a <- (1 - p) / 2
      lo <- qpois(a,     lambda = mu_safe)
      hi <- qpois(1 - a, lambda = mu_safe)
      mean(y >= lo & y <= hi)
    }
    cov50 <- cov_pair_p(0.50); cov80 <- cov_pair_p(0.80); cov95 <- cov_pair_p(0.95)
    logS  <- -mean(dpois(y, lambda = mu_safe, log = TRUE))
    pred_sd <- sqrt(mu_safe)
    DSS  <- mean(((y - mu_safe) / pred_sd)^2 + 2 * log(pred_sd))
    rpit <- randomized_pit_pois(y, mu_safe)
    rpit <- rpit[is.finite(rpit) & rpit >= 0 & rpit <= 1]
    pit_W1 <- if (length(rpit) >= 50) {
      ord <- sort(rpit); m <- length(ord)
      mean(abs(ord - (seq_len(m) - 0.5) / m))
    } else NA_real_
  } else {
    sd_safe <- pmax(sd_fit, 1e-6)
    z <- (y - mu) / sd_safe
    crps <- mean(sd_safe * (z * (2 * pnorm(z) - 1) + 2 * dnorm(z) - 1 / sqrt(pi)))
    cov_pair_g <- function(p) {
      zq <- qnorm((1 + p) / 2)
      mean(y >= mu - zq * sd_safe & y <= mu + zq * sd_safe)
    }
    cov50 <- cov_pair_g(0.50); cov80 <- cov_pair_g(0.80); cov95 <- cov_pair_g(0.95)
    logS <- -mean(dnorm(y, mean = mu, sd = sd_safe, log = TRUE))
    DSS  <- mean(((y - mu) / sd_safe)^2 + 2 * log(sd_safe))
    pit_W1 <- NA_real_
  }

  sd_safe_g <- pmax(sd_fit, 1e-6)
  cov95_gauss <- mean(y >= mu - 1.96 * sd_safe_g & y <= mu + 1.96 * sd_safe_g)

  # Calibration intercept/slope: log E[Y] = a + b log(mu_hat).
  # Under perfect calibration a=0, b=1. Quasi-Poisson absorbs over-dispersion
  # without spuriously inflating SEs that we never report here.
  cal_int <- NA_real_; cal_slope <- NA_real_
  if (n >= 30 && sd(log(mu_safe)) > 1e-6) {
    fit_cal <- tryCatch(suppressWarnings(
      stats::glm(y ~ log(mu_safe), family = stats::quasipoisson(link = "log"))
    ), error = function(e) NULL)
    if (!is.null(fit_cal) && length(coef(fit_cal)) == 2L &&
        all(is.finite(coef(fit_cal)))) {
      cal_int   <- unname(coef(fit_cal)[1])
      cal_slope <- unname(coef(fit_cal)[2])
    }
  }

  data.table(RMSE = rmse, MAE = mae, R2 = r2, CRPS = crps,
             Cov50 = cov50, Cov80 = cov80, Cov95 = cov95,
             Cov95_gauss = cov95_gauss,
             logS = logS, DSS = DSS,
             cal_intercept = cal_int, cal_slope = cal_slope,
             pit_W1 = pit_W1, N = n)
}

#' Stratified scoring: apply cv_scores() within each level of `group`.
#' Used for borough / urban-class / wave-phase tables.
cv_scores_stratified <- function(y, mu, sd_fit, nb_size = NULL, group,
                                 group_name = "group", min_n = 20L,
                                 family = c("auto", "nbinomial", "poisson", "gaussian")) {
  family <- match.arg(family)
  stopifnot(length(group) == length(y))
  g <- as.character(group)
  out <- lapply(unique(g[!is.na(g)]), function(lv) {
    idx <- which(g == lv)
    if (length(idx) < min_n) return(NULL)
    s <- cv_scores(y[idx], mu[idx], sd_fit[idx], nb_size = nb_size, family = family)
    s[, (group_name) := lv]
    s
  })
  out <- out[!vapply(out, is.null, logical(1))]
  if (length(out) == 0) return(data.table())
  setcolorder(rbindlist(out, fill = TRUE), group_name)
}

#' LOO diagnostic on INLA CPO output. INLA's CPO_i = p(y_i | y_-i) so
#' sum(log CPO_i) = ELPD_loo. Pareto-k diagnostics (Vehtari, Gelman & Gabry
#' 2017) are NOT applicable here: they require S x N posterior log-likelihood
#' draws for PSIS importance reweighting, whereas INLA's CPO is a
#' deterministic Laplace-style approximation. We report ELPD_loo and SE only.
cpo_loo_diagnostic <- function(cpo_vec, failure_vec = NULL) {
  ok <- !is.na(cpo_vec) & cpo_vec > 0
  if (!is.null(failure_vec)) ok <- ok & (failure_vec == 0 | is.na(failure_vec))
  log_cpo <- log(cpo_vec[ok])

  list(elpd_loo    = sum(log_cpo),
       se_elpd_loo = sqrt(length(log_cpo) * var(log_cpo)),
       n_used      = length(log_cpo),
       k_max        = NA_real_,
       k_med        = NA_real_,
       pct_k_gt_0_7 = NA_real_,
       pct_k_gt_1_0 = NA_real_)
}

#' Area-level scoring: aggregated counts (+ rates per 100K when offset given).
#' Tolerant to two column conventions:
#'   - {obs_count, pred_count, total_E}     (legacy)
#'   - {obs, pred} with optional total_E    (10b spatial CV)
#' Emits both legacy *_count / *_rate columns and short {R2, RMSE, MAE}
#' aliases so downstream callers that average across folds can pull either.
cv_scores_area <- function(area_dt) {
  y_c <- if ("obs_count"  %in% names(area_dt)) area_dt$obs_count  else area_dt$obs
  p_c <- if ("pred_count" %in% names(area_dt)) area_dt$pred_count else area_dt$pred
  E_a <- if ("total_E"    %in% names(area_dt)) pmax(area_dt$total_E, 1) else NULL
  clip_hi <- max(quantile(c(y_c, p_c), 0.999), max(y_c) * 2, 1)
  n_clip  <- sum(p_c > clip_hi)
  p_c_w   <- pmin(p_c, clip_hi)
  n <- length(y_c)
  RMSE_c <- sqrt(mean((y_c - p_c_w)^2))
  MAE_c  <- mean(abs(y_c - p_c_w))
  R2_c   <- if (n >= 3 && sd(y_c) > 0) cor(y_c, p_c)^2 else NA_real_
  rate <- !is.null(E_a)
  y_r <- if (rate) y_c   / E_a * 1e5 else NA_real_
  p_r <- if (rate) p_c_w / E_a * 1e5 else NA_real_
  data.table(
    n_areas    = n,
    # Legacy column names (count + rate).
    RMSE_count = RMSE_c, MAE_count = MAE_c, R2_count = R2_c,
    RMSE_rate  = if (rate) sqrt(mean((y_r - p_r)^2)) else NA_real_,
    MAE_rate   = if (rate) mean(abs(y_r - p_r))     else NA_real_,
    R2_rate    = if (rate && n >= 3 && sd(y_r) > 0) cor(y_r, p_r)^2 else NA_real_,
    Corr_rate  = if (rate && n >= 3) cor(y_r, p_r, use = "complete") else NA_real_,
    # Short aliases — 10b spatial CV averages these across folds.
    R2   = R2_c, RMSE = RMSE_c, MAE = MAE_c,
    n_clip     = n_clip)
}

#' Randomized PIT for discrete Poisson data (Czado et al. 2009)
randomized_pit_pois <- function(y, mu) {
  mu_safe <- pmax(mu, 1e-10)
  F_y    <- ppois(y,     lambda = mu_safe)
  F_y_m1 <- ppois(y - 1, lambda = mu_safe)
  F_y_m1[y == 0] <- 0
  u <- stats::runif(length(y))
  F_y_m1 + u * (F_y - F_y_m1)
}

#' Randomized PIT for discrete NB data (Czado et al. 2009)
randomized_pit_nb <- function(y, mu, nb_size) {
  mu_safe <- pmax(mu, 1e-10)
  F_y    <- pnbinom(y,     size = nb_size, mu = mu_safe)
  F_y_m1 <- pnbinom(y - 1, size = nb_size, mu = mu_safe)
  F_y_m1[y == 0] <- 0
  V <- runif(length(y))
  pmin(pmax(F_y_m1 + V * (F_y - F_y_m1), 0), 1)
}

# ═══════════════════════════════════════════════════════════════════════════════
# Stratification covariates for stratified CV scoring
# ═══════════════════════════════════════════════════════════════════════════════

#' NYC MODZCTA -> borough lookup (first three digits of the 5-digit code).
#' Returns a character vector of length(geo_id). Unknown -> "Other".
modzcta_to_borough <- function(geo_id) {
  pre <- substr(as.character(geo_id), 1, 3)
  out <- rep("Other", length(geo_id))
  out[pre %in% c("100", "101", "102")]               <- "Manhattan"
  out[pre == "103"]                                   <- "Staten Island"
  out[pre == "104"]                                   <- "Bronx"
  out[pre == "112"]                                   <- "Brooklyn"
  out[pre %in% c("110", "111", "113", "114", "116")] <- "Queens"
  out
}

#' Augment md$df with stratification covariates used for stratified scoring:
#'   - area_strat  : NY borough (from modzcta_to_borough) or UT county_name
#'   - urban_class : building-density tertile (Rural / Suburban / Urban)
#'   - wave_phase  : pre-peak / peak / post-peak from smoothed total daily Y
#'   - cell_E_q    : offset quintile (E_inla)
#' Idempotent (skips columns already present). UT county lookup uses md$grid_sf.
add_cv_strata <- function(md, st) {
  df <- as.data.table(md$df)

  # ── 1. Area stratum (borough / county) ────────────────────────────────────
  if (!"area_strat" %in% names(df)) {
    if (st == "NY" && "geo_id" %in% names(df)) {
      df[, area_strat := modzcta_to_borough(geo_id)]
    } else if (st == "UT" && !is.null(md$grid_sf)) {
      county_lkp <- data.table(
        grid_id     = as.character(md$grid_sf$grid_id),
        area_strat  = as.character(md$grid_sf$county_name))
      df[, grid_id := as.character(grid_id)]
      df <- merge(df, county_lkp, by = "grid_id", all.x = TRUE)
      df[is.na(area_strat), area_strat := "Other"]
    } else {
      df[, area_strat := "all"]
    }
  }

  # ── 2. Urban class (tertile of building_density, fallback log_pop_dens) ───
  if (!"urban_class" %in% names(df)) {
    bd_col <- if ("building_density" %in% names(df)) "building_density"
              else if ("log_pop_dens"     %in% names(df)) "log_pop_dens"
              else NULL
    if (!is.null(bd_col)) {
      bd <- df[[bd_col]]
      qs <- quantile(bd, c(1/3, 2/3), na.rm = TRUE, names = FALSE)
      df[, urban_class := fifelse(bd <= qs[1], "Rural",
                          fifelse(bd <= qs[2], "Suburban", "Urban"))]
    } else {
      df[, urban_class := "all"]
    }
  }

  # ── 3. Wave phase from smoothed total daily Y (pre / peak / post) ─────────
  if (!"wave_phase" %in% names(df)) {
    daily <- df[, .(Y_tot = sum(Y, na.rm = TRUE)), by = idx_temporal][order(idx_temporal)]
    if (nrow(daily) >= 14) {
      # 7-day centred moving mean (NA-padded edges) -> peak day
      k <- 7L
      filt <- stats::filter(daily$Y_tot, rep(1 / k, k), sides = 2)
      peak_t <- daily$idx_temporal[which.max(filt)]
      n_t    <- max(daily$idx_temporal)
      lo_t   <- pmax(1L, as.integer(round(0.7 * peak_t)))
      hi_t   <- pmin(n_t, as.integer(round(1.3 * peak_t)))
      df[, wave_phase := fifelse(idx_temporal < lo_t, "pre",
                          fifelse(idx_temporal > hi_t, "post", "peak"))]
    } else {
      df[, wave_phase := "all"]
    }
  }

  # ── 4. Offset quintile (light-touch diagnostic for offset-driven bias) ────
  if (!"cell_E_q" %in% names(df) && "E_inla" %in% names(df)) {
    qE <- quantile(df$E_inla, seq(0, 1, by = 0.2), na.rm = TRUE, names = FALSE)
    qE <- unique(qE)
    if (length(qE) >= 3) {
      df[, cell_E_q := paste0("Q", findInterval(E_inla, qE, all.inside = TRUE))]
    }
  }

  md$df <- as.data.frame(df)
  md
}

# ═══════════════════════════════════════════════════════════════════════════════
# Data Handling
# ═══════════════════════════════════════════════════════════════════════════════

#' Clean raw gt_dt: enforce minimum cases and exclude counties.
#' Note: no Y truncation — the NB likelihood handles tails; capping observations
#' would diverge the CV target from the in-sample target (09_fit_model does not
#' cap Y) and silently inflate apparent CV skill on Omicron-peak cells.
clean_cv_data <- function(gt_dt, min_cases = MIN_CASES,
                          outlier_q = OUTLIER_QUANTILE,
                          exclude_ids = character(0)) {
  dt <- copy(gt_dt)
  if (length(exclude_ids) > 0) dt <- dt[!grid_id %in% exclude_ids]
  total_by_grid <- dt[, .(total = sum(Y, na.rm = TRUE)), by = grid_id]
  keep <- total_by_grid[total >= min_cases, grid_id]
  dt <- dt[grid_id %in% keep]
  dt
}

#' Set Y = NA for held-out observations, preserving truth in Y_true
mask_holdout <- function(df, holdout_idx) {
  df$Y_true <- df$Y
  df$Y[holdout_idx] <- NA_integer_
  df
}

#' Identify spatial units (grid_id or geo_id) with enough case signal to be
#' meaningfully held out and scored. Sparse rural units (UT_wave1 was 15/25
#' counties with <50 cases over 92 days, >85% zero-days) destabilise the
#' NB hyperpar likelihood under CV masking AND contribute pure noise to
#' CRPS/Cov95 aggregates. We KEEP them in TRAINING (they remain observed
#' across the holdout window) but EXCLUDE them from the holdout indices,
#' so the model never has to predict them and the scoring is not diluted
#' by structural-zero cells.
#'
#' Thresholds default to total < 50 cases OR pct_zero > 0.80; both
#' overridable via options("cv_eval_total_thresh", "cv_eval_pct_zero_thresh")
#' or function args. Returns a character vector of unit IDs that pass.
compute_eval_units <- function(df,
                                unit_col = "grid_id",
                                y_col = "Y",
                                total_thresh = getOption("cv_eval_total_thresh", 50L),
                                pct_zero_thresh = getOption("cv_eval_pct_zero_thresh", 0.80)) {
  if (!unit_col %in% names(df)) return(character(0))
  dt <- as.data.table(df)
  # Use raw Y; df is pre-masking so Y == Y_true here.
  # Aggregate to (unit, day) before computing pct_zero so the metric
  # reflects days-without-cases rather than cell-day sparsity.
  day_col <- intersect(c("sim_day", "time_id", "t", "day"), names(dt))[1]
  if (!is.na(day_col)) {
    per_day <- dt[, .(y_sum = sum(get(y_col), na.rm = TRUE)),
                  by = c(unit_col, day_col)]
    per <- per_day[, .(total = sum(y_sum, na.rm = TRUE),
                       pct_zero = mean(y_sum == 0, na.rm = TRUE)),
                   by = c(unit_col)]
  } else {
    per <- dt[, .(total = sum(get(y_col), na.rm = TRUE),
                  pct_zero = mean(get(y_col) == 0, na.rm = TRUE)),
              by = c(unit_col)]
  }
  ok <- per[total >= total_thresh & pct_zero <= pct_zero_thresh][[unit_col]]
  as.character(ok)
}

#' NB-aware predictive variance for a Poisson-NB GLM with mean mu and rate sd.
#' Falls back to Gaussian sd_fit^2 when nb_size is missing/invalid.
nb_pred_var <- function(mu, sd_fit, nb_size) {
  if (!is.null(nb_size) && is.finite(nb_size) && nb_size > 0) {
    mu + mu^2 / nb_size
  } else {
    sd_fit^2
  }
}

# =============================================================================
# Out-of-fold attribution (cases attributable to mobility)
# =============================================================================

#' Compute out-of-fold mobility attribution (PCAtM) from a fitted INLA model.
#'
#' Returns cell-level attribution and fold-level summaries. Produces both the
#' exponential-counterfactual and a linearised (Taylor) alternative plus
#' winsorised/clip variants to guard against heavy-tailed eta draws.
#' Falls back to independent Gaussian draws when INLA config samples are
#' unavailable (flagged via `approx_indep`).
compute_attribution_oof <- function(fit_result, df_cv, holdout_idx,
                                    mob_pattern = "^(mob_lag|exposure_import)",
                                    n_draws = 200L,
                                    quantiles = c(0.05, 0.50, 0.95)) {
  fe <- fit_result$summary.fixed
  mob_covs <- grep(mob_pattern, rownames(fe), value = TRUE)
  if (length(mob_covs) == 0) return(NULL)

  # Held-out design matrix; NAs zeroed to mirror the INLA fit.
  X_oof <- as.matrix(df_cv[holdout_idx, mob_covs, drop = FALSE])
  X_oof[is.na(X_oof)] <- 0
  mu_hat  <- fit_result$summary.fitted.values$mean[holdout_idx]
  mu_safe <- pmax(mu_hat, 1e-12)

  # ── Posterior draws of beta_mob ──
  beta_draws  <- NULL
  approx_flag <- FALSE
  if (n_draws > 0L && requireNamespace("INLA", quietly = TRUE) &&
      !is.null(fit_result$misc$configs)) {
    sel <- setNames(as.list(rep(1L, length(mob_covs))), mob_covs)
    smps <- tryCatch(
      INLA::inla.posterior.sample(as.integer(n_draws), fit_result,
                                  selection = sel, num.threads = 1L,
                                  verbose = FALSE),
      error = function(e) { warning("inla.posterior.sample failed: ", conditionMessage(e)); NULL })
    if (!is.null(smps)) {
      lat_names <- rownames(smps[[1]]$latent)
      rows <- match(paste0(mob_covs, ":1"), lat_names)
      if (any(is.na(rows))) rows <- match(mob_covs, lat_names)
      if (!any(is.na(rows))) {
        beta_draws <- t(vapply(smps, function(s) s$latent[rows, 1],
                                numeric(length(mob_covs))))
        if (length(mob_covs) == 1L)
          beta_draws <- matrix(beta_draws, ncol = 1L)
        colnames(beta_draws) <- mob_covs
      }
    }
  }
  # Fallback: independent Gaussian (under-estimates joint correlation).
  if (is.null(beta_draws)) {
    bm <- fe[mob_covs, "mean"]; bs <- fe[mob_covs, "sd"]
    M  <- max(as.integer(n_draws), 200L)
    beta_draws <- mapply(function(m, s) rnorm(M, m, s), bm, bs)
    if (!is.matrix(beta_draws)) beta_draws <- matrix(beta_draws, nrow = M)
    colnames(beta_draws) <- mob_covs
    approx_flag <- TRUE
  }
  M <- nrow(beta_draws)

  # ── Per-draw counterfactual mu and attributable cases ──
  # eta_draws: n_cells x M  (X_oof %*% t(beta_draws))
  eta_draws      <- X_oof %*% t(beta_draws)
  mu_nomob_draws <- mu_safe * exp(-eta_draws)
  delta_draws    <- mu_safe - mu_nomob_draws
  # Linearised (Taylor) attribution: bounded, symmetric in eta.
  # Robust under heavy-tail standardised covariates (e.g. UT_w3 exposure_import_std).
  delta_lin_draws <- mu_safe * eta_draws

  # Winsorised (1%/99% per fold across cells, pooled over draws) eta. UT_w3
  # pathology fix: caps the cell-day eta tail that sends 1 - exp(-eta) into
  # non-physical regimes while preserving the bulk distribution and the
  # exponential link. Per-cell across-draw quantiles are computed once.
  q_lo <- as.numeric(quantile(eta_draws, 0.01, na.rm = TRUE))
  q_hi <- as.numeric(quantile(eta_draws, 0.99, na.rm = TRUE))
  eta_draws_w     <- pmin(pmax(eta_draws, q_lo), q_hi)
  delta_w_draws   <- mu_safe - mu_safe * exp(-eta_draws_w)

  # Fixed-bound clip eta in [-2, 2] (UT_w2 pathology fix). Bounds mobility's
  # multiplicative risk effect to e^2 ~ 7x, robust to per-fold quantile drift
  # under degenerate covariate columns where 1% and 99% quantiles can collapse.
  ETA_CLIP_LO     <- -2
  ETA_CLIP_HI     <-  2
  eta_draws_clip2 <- pmin(pmax(eta_draws, ETA_CLIP_LO), ETA_CLIP_HI)
  delta_clip2_draws <- mu_safe - mu_safe * exp(-eta_draws_clip2)

  # Degeneracy flag: >50% of mob_cov columns at near-zero variance means the
  # holdout cells carry no usable mobility signal (typically imputed-floor).
  # Such folds should be excluded from the robust wave-summary aggregation.
  DEGEN_VAR_TOL  <- 1e-6
  DEGEN_FRAC_TOL <- 0.50
  col_vars       <- apply(X_oof, 2L, stats::var, na.rm = TRUE)
  n_degen_cols   <- sum(col_vars < DEGEN_VAR_TOL, na.rm = TRUE)
  degen_frac     <- n_degen_cols / max(ncol(X_oof), 1L)
  is_degenerate  <- degen_frac > DEGEN_FRAC_TOL

  # Diagnostic: track how often eta gets large enough that exponential and
  # linearised forms diverge by > 50%. Surfaced into pcatm_summary so the
  # thesis can cite when to prefer linearised.
  eta_extreme_frac <- mean(abs(eta_draws) > 1, na.rm = TRUE)

  row_q <- function(mat, p)
    apply(mat, 1L, function(v) as.numeric(quantile(v, probs = p, na.rm = TRUE)))

  qlo <- quantiles[1]; qmed <- quantiles[2]; qhi <- quantiles[3]
  cell_attr <- data.table::data.table(
    holdout_pos     = seq_along(holdout_idx),
    mu              = mu_hat,
    mu_nomob_mean   = rowMeans(mu_nomob_draws),
    mu_nomob_lo     = row_q(mu_nomob_draws, qlo),
    mu_nomob_hi     = row_q(mu_nomob_draws, qhi),
    delta_mean      = rowMeans(delta_draws),
    delta_lo        = row_q(delta_draws, qlo),
    delta_hi        = row_q(delta_draws, qhi),
    # Linearised (bounded) alternative — symmetric in eta, no exp blow-up.
    delta_lin_mean  = rowMeans(delta_lin_draws),
    delta_lin_lo    = row_q(delta_lin_draws, qlo),
    delta_lin_hi    = row_q(delta_lin_draws, qhi))

  # ── Wave/fold-level PCAtM (draw-wise aggregation, preserves CrIs) ──
  total_mu            <- sum(mu_safe)
  total_delta_dr      <- colSums(delta_draws)
  total_delta_lin_dr  <- colSums(delta_lin_draws)
  total_delta_w_dr    <- colSums(delta_w_draws)
  total_delta_clip2_dr <- colSums(delta_clip2_draws)
  pcatm_draws         <- total_delta_dr      / max(total_mu, 1e-12)
  pcatm_lin_draws     <- total_delta_lin_dr  / max(total_mu, 1e-12)
  pcatm_w_draws       <- total_delta_w_dr    / max(total_mu, 1e-12)
  pcatm_clip2_draws   <- total_delta_clip2_dr / max(total_mu, 1e-12)
  pcatm_summary   <- data.table::data.table(
    PCAtM_mean         = mean(pcatm_draws),
    PCAtM_lo           = as.numeric(quantile(pcatm_draws, qlo, na.rm = TRUE)),
    PCAtM_med          = as.numeric(quantile(pcatm_draws, qmed, na.rm = TRUE)),
    PCAtM_hi           = as.numeric(quantile(pcatm_draws, qhi, na.rm = TRUE)),
    attr_cases_mean    = mean(total_delta_dr),
    attr_cases_lo      = as.numeric(quantile(total_delta_dr, qlo, na.rm = TRUE)),
    attr_cases_hi      = as.numeric(quantile(total_delta_dr, qhi, na.rm = TRUE)),
    # Linearised (Slater-style) PCAtM — preferred when eta_extreme_frac > 0.05.
    PCAtM_lin_mean        = mean(pcatm_lin_draws),
    PCAtM_lin_lo          = as.numeric(quantile(pcatm_lin_draws, qlo, na.rm = TRUE)),
    PCAtM_lin_med         = as.numeric(quantile(pcatm_lin_draws, qmed, na.rm = TRUE)),
    PCAtM_lin_hi          = as.numeric(quantile(pcatm_lin_draws, qhi, na.rm = TRUE)),
    attr_cases_lin_mean   = mean(total_delta_lin_dr),
    attr_cases_lin_lo     = as.numeric(quantile(total_delta_lin_dr, qlo, na.rm = TRUE)),
    attr_cases_lin_hi     = as.numeric(quantile(total_delta_lin_dr, qhi, na.rm = TRUE)),
    # Winsorised (1%/99%) PCAtM — UT_w3 pathology fix. Bounded variant of the
    # exponential link; preferred when eta_extreme_frac > 0.05 and the
    # linearised form is too conservative.
    PCAtM_w99_mean        = mean(pcatm_w_draws),
    PCAtM_w99_lo          = as.numeric(quantile(pcatm_w_draws, qlo,  na.rm = TRUE)),
    PCAtM_w99_med         = as.numeric(quantile(pcatm_w_draws, qmed, na.rm = TRUE)),
    PCAtM_w99_hi          = as.numeric(quantile(pcatm_w_draws, qhi,  na.rm = TRUE)),
    # Fixed-bound clip [-2, 2] PCAtM — robust to per-fold quantile drift.
    # Reported as the canonical OOF PCAtM in the thesis (Section 5.4).
    PCAtM_clip2_mean      = mean(pcatm_clip2_draws),
    PCAtM_clip2_lo        = as.numeric(quantile(pcatm_clip2_draws, qlo,  na.rm = TRUE)),
    PCAtM_clip2_med       = as.numeric(quantile(pcatm_clip2_draws, qmed, na.rm = TRUE)),
    PCAtM_clip2_hi        = as.numeric(quantile(pcatm_clip2_draws, qhi,  na.rm = TRUE)),
    eta_w99_lo            = q_lo,
    eta_w99_hi            = q_hi,
    eta_clip2_lo          = ETA_CLIP_LO,
    eta_clip2_hi          = ETA_CLIP_HI,
    eta_extreme_frac      = eta_extreme_frac,
    n_degen_cols          = n_degen_cols,
    degen_frac            = degen_frac,
    is_degenerate         = is_degenerate,
    total_mu              = total_mu,
    n_oof                 = length(holdout_idx),
    n_mob_covs            = length(mob_covs),
    n_draws               = M,
    approx_indep          = approx_flag)

  # ── Beta-mob summary (coefficient stability tracker) ──
  beta_summary <- data.table::data.table(
    covariate = mob_covs,
    mean = fe[mob_covs, "mean"],
    sd   = fe[mob_covs, "sd"],
    q05  = apply(beta_draws, 2L, quantile, probs = qlo,  na.rm = TRUE),
    q50  = apply(beta_draws, 2L, quantile, probs = qmed, na.rm = TRUE),
    q95  = apply(beta_draws, 2L, quantile, probs = qhi,  na.rm = TRUE))

  list(cell_attr     = cell_attr,
       pcatm_summary = pcatm_summary,
       beta_summary  = beta_summary,
       beta_draws    = beta_draws,   # M x p_mob, persisted for re-aggregation
       X_oof         = X_oof,        # n_oof x p_mob, persisted for re-aggregation
       mob_covs      = mob_covs)
}

#' Fit best-model (or supplied formula) with a held-out subset masked, return
#' predictions and NB size on the holdout. Common scaffold for temporal,
#' MODZCTA, and spatial CV. `model_name` selects the offset variant in get_E()
#' (Base/Temp use E_nomob; mobility models use E_inla).
#'
#' Also returns the posterior hyperparameter mode (`theta`) so the caller can
#' warm-start the next fold, and a `holdout_dt` carrying any stratification
#' columns present in `md$df` for stratified scoring downstream.
fit_holdout <- function(md, fm, holdout_idx, label, warmstart_theta = NULL,
                        model_name = NULL, n_attr_draws = 200L,
                        compute_cpo = TRUE) {
  if (is.null(model_name)) model_name <- md$best_model %||% "NLMob_ST_IV"
  df_cv <- mask_holdout(md$df, holdout_idx)
  E_vec <- get_E(model_name, df_cv)
  fit <- fit_model(fm, df_cv, E_vec,
                   model_name = label,
                   zero_frac = md$zero_frac, dispersion = md$dispersion,
                   cpo = compute_cpo, config = (n_attr_draws > 0L), cv_mode = TRUE,
                   n_grids = md$n_grids,
                   warmstart_theta = warmstart_theta)
  if (is.null(fit$result)) return(list(ok = FALSE, runtime = 0))
  mu     <- fit$result$summary.fitted.values$mean[holdout_idx]
  sd_fit <- fit$result$summary.fitted.values$sd[holdout_idx]
  y_true <- df_cv$Y_true[holdout_idx]
  hp <- fit$result$summary.hyperpar
  nb_size_row <- grep("size.*nbinomial", rownames(hp), ignore.case = TRUE)
  nb_size <- if (length(nb_size_row) > 0) hp[nb_size_row[1], "mean"] else NULL

  # ── CPO / LCPO over the TRAINING rows for this fold. Masked-holdout rows
  # carry NA Y so INLA never computes a CPO for them; the resulting vector
  # is the leave-one-out predictive density on the training half. We report
  # LCPO = -mean(log(cpo)) and a failure count for downstream diagnostics.
  cpo_train_lcpo     <- NA_real_
  cpo_train_failures <- NA_integer_
  cpo_train_n        <- NA_integer_
  if (!is.null(fit$result$cpo) && length(fit$result$cpo$cpo)) {
    train_mask <- rep(TRUE, length(fit$result$cpo$cpo))
    train_mask[holdout_idx] <- FALSE
    cpo_vec  <- fit$result$cpo$cpo[train_mask]
    fail_vec <- fit$result$cpo$failure[train_mask]
    ok_cpo   <- !is.na(cpo_vec) & cpo_vec > 0
    if (!is.null(fail_vec))
      ok_cpo <- ok_cpo & (is.na(fail_vec) | fail_vec == 0)
    if (any(ok_cpo)) {
      cpo_train_lcpo     <- -mean(log(pmax(cpo_vec[ok_cpo], 1e-300)))
      cpo_train_n        <- as.integer(sum(ok_cpo))
      cpo_train_failures <- as.integer(
        sum(!is.na(fail_vec) & fail_vec > 0, na.rm = TRUE))
    }
  }

  # Family used by this fold's fit_model call. Drives scoring branch in
  # cv_scores (Poisson vs NB vs Gaussian fallback). Pulled straight from
  # the INLA call to avoid second-guessing via nb_size presence alone.
  fold_family <- fit$result$.args$family %||%
                 (if (!is.null(nb_size)) "nbinomial" else "poisson")

  # Posterior hyperparameter mode for sequential cross-fold warmstart + C1
  # coefficient-stability logging.
  theta_post <- tryCatch(fit$result$mode$theta, error = function(e) NULL)
  if (is.null(theta_post) || !is.numeric(theta_post)) theta_post <- NULL
  fixed_full <- fit$result$summary.fixed   # full beta posterior summary (C1)

  # ── A1+A2: out-of-fold attribution (counterfactual mu + PCAtM with CrIs). ──
  attribution <- tryCatch(
    compute_attribution_oof(fit$result, df_cv, holdout_idx,
                            n_draws = n_attr_draws),
    error = function(e) {
      warning("compute_attribution_oof failed: ", conditionMessage(e)); NULL })

  # Holdout long table (one row per held-out cell-day) with strata columns
  # available in md$df. Used for stratified scoring; cheap to assemble here.
  strat_cols <- intersect(c("grid_id", "geo_id", "idx_temporal", "idx_spatial",
                            "area_strat", "urban_class", "wave_phase",
                            "cell_E_q"),
                          names(df_cv))
  holdout_dt <- as.data.table(df_cv[holdout_idx, strat_cols, drop = FALSE])
  holdout_dt[, `:=`(y_true = y_true, mu = mu, sd_fit = sd_fit)]

  # Memory hygiene: drop the heaviest INLA marginals after extraction. The
  # fit object is still referenced for `.args$family` above; remaining slots
  # (summary.fixed, summary.fitted.values, mode) are small.
  try({ fit$result$marginals.fitted.values <- NULL
        fit$result$marginals.random        <- NULL
        fit$result$marginals.linear.predictor <- NULL }, silent = TRUE)

  list(ok = TRUE, mu = mu, sd_fit = sd_fit, y_true = y_true,
       nb_size = nb_size, family = fold_family,
       runtime = fit$runtime,
       time_id = df_cv$idx_temporal[holdout_idx],
       fixed   = fixed_full,
       fixed_full = fixed_full,
       theta   = theta_post,
       cpo_train_lcpo     = cpo_train_lcpo,
       cpo_train_failures = cpo_train_failures,
       cpo_train_n        = cpo_train_n,
       attribution = attribution,
       holdout_dt = holdout_dt)
}

#' Rebuild unified data + formula for ST_IV, mirroring 09_fit_model.Rmd
rebuild_model_data <- function(gt_dt, out_dir, wave_id = NULL) {
  grid_sf <- st_read(file.path(CONFIG$grids_dir, sprintf("grid_100m_%s.gpkg", state)),
                     quiet = TRUE)
  grid_sf <- grid_sf[grid_sf$grid_id %in% unique(gt_dt$grid_id), ]

  aug <- build_augmented_adjacency(grid_sf, bridge_links, WORMHOLES, gt_dt = gt_dt)
  nb_aug <- aug$nb

  comp_check <- safe_n_comp_nb(nb_aug)
  if (comp_check$nc > 1) {
    cs <- table(comp_check$comp.id)
    keep <- which(comp_check$comp.id == as.integer(names(which.max(cs))))
    grid_sf <- grid_sf[keep, ]
    gt_dt <- gt_dt[grid_id %in% grid_sf$grid_id]
    nb_keep <- nb_aug[keep]
    old_to_new <- setNames(seq_along(keep), keep)
    nb_aug <- lapply(nb_keep, function(nbrs)
      as.integer(old_to_new[as.character(nbrs[nbrs %in% keep])]))
    attr(nb_aug, "class") <- "nb"
    attr(nb_aug, "region.id") <- as.character(seq_along(keep))
    nb_aug <- sanitize_nb(nb_aug)
  }

  graph_file <- file.path(out_dir, "models", "grid_adj_unified.graph")
  safe_nb2INLA(graph_file, nb_aug)
  n_grids <- length(unique(grid_sf$grid_id))

  df <- as.data.table(gt_dt)[order(grid_id, time_id)]
  gi <- data.table(grid_id = sort(unique(df$grid_id)),
                   idx_spatial = seq_len(uniqueN(df$grid_id)))
  if ("idx_spatial" %in% names(df)) df[, idx_spatial := NULL]
  df <- merge(df, gi, by = "grid_id", all.x = TRUE)
  df[, idx_temporal := as.integer(factor(time_id))]
  df[, idx_interaction := .I]

  # Match 09_fit_model exactly: K-clusters = min(N_ST_CLUSTERS, n_grids/100).
  # Without this, CV used min(N_ST_CLUSTERS, n_grids) (effectively
  # N_ST_CLUSTERS) and the cluster_graph dimension drifted from the in-sample
  # fit's. Also honour the same fit_n_st_clusters env var that 09 reads, so
  # SLURM jobs can pin K to the value used by the original fit.
  .cv_n_clusters_env <- suppressWarnings(
    as.integer(Sys.getenv("fit_n_st_clusters", unset = "")))
  .N_K <- if (length(.cv_n_clusters_env) && !is.na(.cv_n_clusters_env) &&
              .cv_n_clusters_env > 0)
            .cv_n_clusters_env else N_ST_CLUSTERS
  .k_match09 <- max(1L, min(.N_K, round(n_grids / 100)))
  st_res <- add_st_clusters(df, grid_sf, out_dir, n_clusters = .k_match09)
  df <- st_res$df
  df <- stabilize_offsets(df)

  has <- function(x) x %in% names(df)
  rho_col <- if (has("rho_norm")) "rho_norm" else if (has("rho")) "rho" else NULL
  if (!is.null(rho_col)) {
    df[[rho_col]] <- pmax(pmin(df[[rho_col]], RHO_TRIM_HIGH), RHO_TRIM_LOW)
    if (!has("log_rho")) df$log_rho <- log(df[[rho_col]])
    if (!has("log_rho_std")) {
      m <- mean(df$log_rho, na.rm = TRUE); s <- sd(df$log_rho, na.rm = TRUE)
      df$log_rho_std <- if (s > 1e-10) (df$log_rho - m) / s else 0
    }
  }
  if (has("contact_intensity") && !has("contact_intensity_std")) {
    s <- sd(df$contact_intensity, na.rm = TRUE)
    df$contact_intensity_std <- if (s > 1e-10) scale(df$contact_intensity)[, 1] else 0
  }
  df <- add_rho_bins(df)
  df <- as.data.frame(df)
  df <- df |>
    add_spatial_orthogonalization(grid_sf) |>
    add_temporal_lags() |>
    add_temporal_centering()

  avail <- select_covariates(df, POTENTIAL_COVARIATES)
  flags <- detect_flags(df)
  disp  <- var(df$Y, na.rm = TRUE) / max(mean(df$Y, na.rm = TRUE), 0.01)
  zero_frac <- mean(df$Y == 0, na.rm = TRUE)
  fms <- build_formulas(avail, graph_file, flags, disp,
                        st_res$cluster_graph, n_grids,
                        state = state, wave_id = wave_id)

  # Match 09_fit_model exactly: m1_extras (mob lags + neg control + import) are
  # attached to *every* MOB_MODELS formula, not only ST_IV. CV refits must use
  # the same formula structure as the in-sample fit.
  best_lags <- select_mob_lags(df)
  m1_extras <- c(best_lags, intersect(c("exposure_import_std", "mob_lag_neg"), names(df)))
  m1_extras <- m1_extras[sapply(m1_extras, function(v) sd(df[[v]], na.rm = TRUE) > 0.01)]
  add_extras <- function(flist) {
    if (length(m1_extras) == 0) return(flist)
    mob_models <- intersect(MOB_MODELS, names(flist))
    for (mn in mob_models)
      flist[[mn]] <- update(flist[[mn]],
                            as.formula(paste("~ . +", paste(m1_extras, collapse = " + "))))
    flist
  }
  fms <- add_extras(fms)

  if ("NLMob_ST_IV" %in% names(fms) && !is.null(st_res$cluster_graph)) {
    merge_res <- merge_degenerate_clusters(df, st_res$cluster_graph)
    if (merge_res$changed) {
      df <- merge_res$df
      fms <- build_formulas(avail, graph_file, flags, disp,
                            st_res$cluster_graph, n_grids,
                            state = state, wave_id = wave_id)
      fms <- add_extras(fms)
    }
  }

  list(df = df, grid_sf = grid_sf, graph_file = graph_file,
       fms = fms, dispersion = disp, zero_frac = zero_frac,
       n_grids = n_grids, n_days = max(df$idx_temporal),
       grid_index = gi, nb = nb_aug, cluster_graph = st_res$cluster_graph)
}

# ═══════════════════════════════════════════════════════════════════════════════
# Transport / Adjacency Setup
# ═══════════════════════════════════════════════════════════════════════════════

#' Configure transport globals for adjacency builder per state
setup_transport_globals <- function(st) {
  TIGER_BRIDGE_TYPES <<- c("S1100", "S1200")
  if (st == "NY") {
    BRIDGE_CSV <- CONFIG$bridge_csv
    SUBWAY_ROUTES_SHP   <<- CONFIG$subway_routes_shp
    SUBWAY_STATIONS_CSV <<- CONFIG$subway_stations_csv
    ROAD_SHP            <<- CONFIG$ny_roads_shp
    if (file.exists(BRIDGE_CSV)) {
      bridges <- fread(BRIDGE_CSV, na.strings = c("", "NA"))
      setnames(bridges, c("X - COORD  (LAT)", "Y - COORD   (LON)"),
               c("lat", "lon"), skip_absent = TRUE)
      bridge_links <<- bridges[nchar(BORO) >= 2 & !is.na(lat) & !is.na(lon)]
    } else bridge_links <<- data.table(lat = numeric(), lon = numeric())
  } else {
    SUBWAY_ROUTES_SHP   <<- ""
    SUBWAY_STATIONS_CSV <<- ""
    ROAD_SHP            <<- ""
    bridge_links <<- data.table(lat = numeric(), lon = numeric())
  }
  WORMHOLES <<- if (st == "NY") NYC_WORMHOLES else UT_WORMHOLES
  state <<- st
}

# ═══════════════════════════════════════════════════════════════════════════════
# Fold Creation (Spatial CV)
# ═══════════════════════════════════════════════════════════════════════════════

#' Create area-level folds for spatial CV
#' NYC: K-means on MODZCTA centroids. Utah: LOCO (each county = one fold)
create_area_folds <- function(df, grid_sf, st, k_ny = 10L) {
  if (st == "NY") {
    gm <- unique(data.table(grid_id = as.character(df$grid_id),
                            area_id = as.character(df$geo_id)))
    gm <- gm[area_id != "unknown" & area_id != "99999"]
    cents <- suppressWarnings(st_centroid(grid_sf))
    cents_dt <- data.table(grid_id = as.character(cents$grid_id),
                           x = st_coordinates(cents)[, 1],
                           y = st_coordinates(cents)[, 2])
    gm <- merge(gm, cents_dt, by = "grid_id")
    area_cents <- gm[, .(x = mean(x), y = mean(y), n_grids = .N), by = area_id]
    k_use <- min(k_ny, nrow(area_cents))
    set.seed(42)
    km <- kmeans(area_cents[, .(x, y)], centers = k_use, nstart = 25)
    area_cents[, cv_fold := km$cluster]
    fold_dt <- area_cents[, .(area_id, cv_fold)]
    n_folds <- max(fold_dt$cv_fold)
  } else {
    county_lkp <- data.table(grid_id = as.character(grid_sf$grid_id),
                             county_name = grid_sf$county_name)
    gm <- unique(data.table(grid_id = as.character(df$grid_id),
                            area_id = as.character(df$geo_id)))
    gm <- merge(gm, county_lkp, by = "grid_id", all.x = TRUE)
    gm_valid <- gm[!county_name %in% UT_EXCLUDE_TEST]
    gm_valid <- gm_valid[!is.na(county_name)]
    counties_sorted <- sort(unique(gm_valid$county_name))
    county_fold <- data.table(county_name = counties_sorted,
                              cv_fold = seq_along(counties_sorted))
    gm_valid <- merge(gm_valid, county_fold, by = "county_name")
    fold_dt <- unique(gm_valid[, .(area_id, cv_fold)])
    n_folds <- max(fold_dt$cv_fold, na.rm = TRUE)
  }
  grid_area <- unique(data.table(grid_id = as.character(df$grid_id),
                                 area_id = as.character(df$geo_id)))
  grid_area <- grid_area[area_id != "unknown" & area_id != "99999"]
  list(fold_dt = fold_dt, grid_area = grid_area, n_folds = n_folds)
}

# ═══════════════════════════════════════════════════════════════════════════════
# Utility
# ═══════════════════════════════════════════════════════════════════════════════

#' Extract NB size parameter from saved model results.
#' By default uses the per-wave best model from BEST_MODEL_BY_WAVE; pass
#' `model_name` or `wave_name` to override.
get_nb_size <- function(model_results, wave_name = NULL, model_name = NULL) {
  if (is.null(model_name)) {
    model_name <- (if (!is.null(wave_name)) BEST_MODEL_BY_WAVE[[wave_name]] else NULL) %||%
                  model_results$best_model %||% "NLMob_ST_IV"
  }
  hp <- model_results$model_outputs[[model_name]]$hyperparameters
  if (is.null(hp)) return(NULL)
  nb_row <- hp[grepl("size.*nbinomial", hp$parameter, ignore.case = TRUE), ]
  if (nrow(nb_row) > 0) nb_row$mean[1] else NULL
}

#' Load wave data and set up globals for CV. The returned `md` carries
#' `best_model` (selected by 09_fit_model via WAIC); CV callers use this
#' to look up the formula `md$fms[[best_model]]` and the warmstart file.
load_wave_for_cv <- function(wave_name) {
  st <- sub("_.*", "", wave_name)
  wave_id <- sub(".*_", "", wave_name)
  out_dir <- file.path(CONFIG$output_dir, wave_name)

  setup_transport_globals(st)

  model_data    <- readRDS(file.path(out_dir, "models", "model_data_gt.rds"))
  model_results <- readRDS(file.path(out_dir, "models", "model_results_unified.rds"))

  exc_key <- paste0(st, "_", wave_id)
  exc_ids <- EXCLUDE_COUNTIES[[exc_key]]
  gt_clean <- clean_cv_data(model_data$gt_dt,
                            exclude_ids = if (!is.null(exc_ids)) exc_ids else character(0))
  md <- rebuild_model_data(gt_clean, out_dir, wave_id = wave_id)

  nb_size <- get_nb_size(model_results, wave_name = wave_name)

  best_model <- BEST_MODEL_BY_WAVE[[wave_name]] %||%
                model_results$best_model %||% "NLMob_ST_IV"
  if (!best_model %in% names(md$fms)) {
    warning(sprintf("[%s] best_model '%s' not in rebuilt formulas; falling back to NLMob_ST_IV",
                    wave_name, best_model))
    best_model <- "NLMob_ST_IV"
  }
  md$best_model <- best_model

  # Attach stratification covariates (borough/county, urban_class, wave_phase,
  # cell_E_q) used by stratified CV scoring. Idempotent if already present.
  md <- add_cv_strata(md, st)

  list(md = md, model_results = model_results, nb_size = nb_size,
       state = st, wave_id = wave_id, out_dir = out_dir,
       best_model = best_model)
}
