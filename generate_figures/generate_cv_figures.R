#!/usr/bin/env Rscript
# generate_cv_figures.R
# Cross-validation diagnostic figure generator. Reads saved CV artifacts
# (no refitting) and writes publication-ready diagnostic PNGs.

suppressPackageStartupMessages({
  library(data.table)
  library(sf)
  library(ggplot2)
  library(patchwork)
  library(viridis)
  library(scales)
  library(dplyr)
  library(tidyr)
})

Sys.unsetenv("MallocStackLogging")
source("helpers/_utils.R")
source("helpers/_cv_helpers.R")  # provides BEST_MODEL_BY_WAVE
load_pipeline()

# Helper: per-wave best model (hard-coded in _cv_helpers.R), with fallbacks.
best_model_for <- function(wave, mr = NULL) {
  BEST_MODEL_BY_WAVE[[wave]] %||% (if (!is.null(mr)) mr$best_model else NULL) %||% "NLMob_ST_IV"
}

BASE      <- CONFIG$output_dir
FIG_DIR   <- file.path(BASE, "figures")
LATEX_FIG <- file.path("thesis_latex", "figures")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(LATEX_FIG, recursive = TRUE, showWarnings = FALSE)

WAVES  <- c("NY_wave1", "NY_wave2", "NY_wave3", "UT_wave1", "UT_wave2", "UT_wave3")
STATES <- setNames(sub("_.*", "", WAVES), WAVES)

# Pretty panel title: "NY_wave1" -> "NY Wave 1"
nice_wave <- function(w) {
  parts <- strsplit(w, "_wave", fixed = TRUE)[[1]]
  if (length(parts) == 2) paste0(parts[1], " Wave ", parts[2]) else w
}
NICE_WAVE <- setNames(vapply(WAVES, nice_wave, character(1)), WAVES)

STATE_COLORS <- c(NY = "#E41A1C", UT = "#377EB8")
WAVE_COLORS  <- c(NY_wave1 = "#E41A1C", NY_wave2 = "#FC8D62", NY_wave3 = "#E7298A",
                   UT_wave1 = "#377EB8", UT_wave2 = "#66C2A5", UT_wave3 = "#7570B3")

save_fig <- function(p, fname, w = 17/2.54, h = 12/2.54) {
  safe_save(FIG_DIR, fname, p, w = w, h = h)
  file.copy(file.path(FIG_DIR, fname),
            file.path(LATEX_FIG, fname), overwrite = TRUE)
}

THEME_CV <- theme_minimal(base_size = 10) +
  theme(plot.title    = element_text(face = "bold", hjust = 0.5, size = 12),
        plot.subtitle = element_text(hjust = 0.5, color = "grey40", size = 9),
        strip.text    = element_text(face = "bold", size = 9),
        legend.position = "bottom")

# ═══════════════════════════════════════════════════════════════════════════════
# Load data
# ═══════════════════════════════════════════════════════════════════════════════

cpo_list  <- list()
tcv_list  <- list()
model_list <- list()

for (w in WAVES) {
  f_cpo <- file.path(BASE, w, "validation", "cpo_diagnostics.rds")
  f_tcv <- file.path(BASE, w, "validation", "temporal_cv_results.rds")
  f_mod <- file.path(BASE, w, "models", "model_results_unified.rds")
  if (file.exists(f_cpo)) cpo_list[[w]] <- readRDS(f_cpo)
  if (file.exists(f_tcv)) tcv_list[[w]] <- readRDS(f_tcv)
  if (file.exists(f_mod)) model_list[[w]] <- readRDS(f_mod)
}

# Build tidy CPO summary table
if (length(cpo_list) > 0) {
  cpo_summary <- rbindlist(lapply(names(cpo_list), function(w) cpo_list[[w]]$summary))
  cpo_summary[, state := sub("_.*", "", wave)]
} else {
  message("[cv] no cpo_diagnostics.rds found — skipping CPO-based figures.")
  cpo_summary <- data.table()
}

# Build tidy temporal CV scores table
# NOTE: the temporal_cv_results.rds schema was reorganised; it no longer stores
# per-horizon `$scores`/`$daily` sublists but flat fields (cell_scores,
# window_scores, daily_agg, holdout_long, ...). The score-by-horizon and
# bias-timeseries figures below depend on the old schema; we build an empty
# placeholder so downstream guards can no-op without crashing the whole script.
tcv_has_old_schema <- length(tcv_list) > 0 &&
  is.list(tcv_list[[1]]) && length(tcv_list[[1]]) > 0 &&
  is.list(tcv_list[[1]][[1]]) &&
  all(c("scores", "daily") %in% names(tcv_list[[1]][[1]]))

if (tcv_has_old_schema) {
  tcv_scores <- rbindlist(lapply(names(tcv_list), function(w) {
    rbindlist(lapply(names(tcv_list[[w]]), function(h) {
      s <- copy(tcv_list[[w]][[h]]$scores)
      if (!"wave" %in% names(s))    s[, wave := w]
      if (!"horizon" %in% names(s)) s[, horizon := as.integer(h)]
      s
    }), fill = TRUE)
  }), fill = TRUE)
  tcv_scores[, state := sub("_.*", "", wave)]

  tcv_daily <- rbindlist(lapply(names(tcv_list), function(w) {
    rbindlist(lapply(names(tcv_list[[w]]), function(h) {
      d <- as.data.table(tcv_list[[w]][[h]]$daily)
      d[, `:=`(wave = w, horizon = as.integer(h))]
      d
    }), fill = TRUE)
  }), fill = TRUE)
  tcv_daily[, state := sub("_.*", "", wave)]
} else {
  message("[cv] temporal_cv_results.rds uses new schema — skipping per-horizon score/daily figures.")
  tcv_scores <- data.table()
  tcv_daily  <- data.table()
}

# ═══════════════════════════════════════════════════════════════════════════════
# A) CPO-BASED DIAGNOSTICS
# ═══════════════════════════════════════════════════════════════════════════════
if (nrow(cpo_summary) > 0) {

# ── A1: LCPO bar chart ──
p_lcpo <- ggplot(cpo_summary, aes(x = wave, y = LCPO, fill = state)) +
  geom_col(width = 0.65, alpha = 0.85) +
  geom_text(aes(label = sprintf("%.3f", LCPO)), vjust = -0.5, size = 3) +
  scale_fill_manual(values = STATE_COLORS) +
  labs(title = "In-Sample Predictive Quality (LCPO)",
       subtitle = "Lower = sharper per-observation predictive density",
       x = NULL, y = expression(-mean(log(CPO[i]))),
       fill = "State") +
  THEME_CV +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))
save_fig(p_lcpo, "cv_lcpo_barplot.png", w = 16/2.54, h = 10/2.54)

# ── A2: -log(CPO) density / violin per wave ──
cpo_violin_data <- rbindlist(lapply(names(cpo_list), function(w) {
  vals <- cpo_list[[w]]$cpo
  vals <- vals[!is.na(vals) & vals > 0 & is.finite(vals)]
  neg_log_cpo <- -log(vals)
  # Cap extreme tails for visibility
  cap <- quantile(neg_log_cpo, 0.995)
  neg_log_cpo <- pmin(neg_log_cpo, cap)
  # Subsample for speed
  if (length(neg_log_cpo) > 20000) neg_log_cpo <- sample(neg_log_cpo, 20000)
  data.table(wave = w, state = sub("_.*", "", w), neg_log_cpo = neg_log_cpo)
}))

p_cpo_violin <- ggplot(cpo_violin_data, aes(x = wave, y = neg_log_cpo, fill = state)) +
  geom_violin(alpha = 0.7, scale = "width", linewidth = 0.3) +
  geom_boxplot(width = 0.15, fill = "white", alpha = 0.6, outlier.shape = NA,
               linewidth = 0.3) +
  scale_fill_manual(values = STATE_COLORS) +
  labs(title = expression("Distribution of " * -log(CPO[i]) * " per Wave"),
       subtitle = "Heavier right tail = more poorly predicted observations",
       x = NULL, y = expression(-log(CPO[i])),
       fill = "State") +
  THEME_CV +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))
save_fig(p_cpo_violin, "cv_cpo_violin.png", w = 18/2.54, h = 12/2.54)

# ── A3: CPO spatial map (worst-predicted cells) ──
cpo_map_panels <- list()
for (w in WAVES) {
  st <- sub("_.*", "", w)
  grid_path <- file.path(BASE, "grids", sprintf("grid_100m_%s.gpkg", st))
  if (!file.exists(grid_path) || is.null(cpo_list[[w]]) || is.null(model_list[[w]])) next

  mr <- model_list[[w]]
  cpo_vals <- cpo_list[[w]]$cpo
  fv <- mr$model_outputs[[best_model_for(w, mr)]]$fitted_values

  if (is.null(fv) || nrow(fv) != length(cpo_vals)) next

  # Aggregate CPO to grid level: mean -log(CPO) per grid cell
  fv$neg_log_cpo <- ifelse(cpo_vals > 0 & !is.na(cpo_vals), -log(cpo_vals), NA_real_)
  grid_cpo <- as.data.table(fv)[, .(mean_neg_log_cpo = mean(neg_log_cpo, na.rm = TRUE)),
                                  by = grid_id]

  grid_sf <- st_read(grid_path, quiet = TRUE)
  grid_sf <- merge(grid_sf, grid_cpo, by = "grid_id")

  # Cap for color scale
  cap <- quantile(grid_sf$mean_neg_log_cpo, 0.98, na.rm = TRUE)
  grid_sf$mean_neg_log_cpo_cap <- pmin(grid_sf$mean_neg_log_cpo, cap)

  p <- ggplot(grid_sf) +
    geom_sf(aes(fill = mean_neg_log_cpo_cap), linewidth = 0, color = NA) +
    scale_fill_viridis_c(option = "inferno", name = expression(-log(CPO)),
                         labels = label_number(accuracy = 0.1)) +
    labs(title = NICE_WAVE[w] %||% w) +
    theme_void(base_size = 8) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 10),
          legend.position = "right",
          legend.key.height = unit(0.4, "cm"),
          legend.key.width = unit(0.3, "cm"))

  cpo_map_panels[[w]] <- p
}

if (length(cpo_map_panels) > 0) {
  p_cpo_maps <- wrap_plots(cpo_map_panels, ncol = 3) +
    plot_annotation(
      title = "Spatial Distribution of Predictive Difficulty",
      subtitle = expression("Grid-level mean " * -log(CPO[i]) * " — brighter = harder to predict"),
      theme = theme(plot.title    = element_text(face = "bold", hjust = 0.5, size = 13),
                    plot.subtitle = element_text(hjust = 0.5, color = "grey40", size = 9)))
  save_fig(p_cpo_maps, "cv_cpo_spatial.png",
           w = 24/2.54, h = ceiling(length(cpo_map_panels) / 3) * 10/2.54)
}

} else {
  message("[cv] Section A (CPO diagnostics) skipped — no cpo_list.")
}


# ═══════════════════════════════════════════════════════════════════════════════
# B) TEMPORAL FORWARD CV
# ═══════════════════════════════════════════════════════════════════════════════
# Section B depends on the legacy temporal_cv_results schema (per-horizon
# `$scores`/`$daily`). If absent, skip the whole section without erroring.
if (nrow(tcv_scores) > 0 && nrow(tcv_daily) > 0) {

# ── B4: Observed vs predicted scatter ──
# daily already has lo/hi from the CV pipeline; recalculate only if absent
if (!"lo" %in% names(tcv_daily)) tcv_daily[, lo := pmax(pred - 1.96 * sd, 0)]
if (!"hi" %in% names(tcv_daily)) tcv_daily[, hi := pred + 1.96 * sd]

p_scatter <- ggplot(tcv_daily, aes(x = obs, y = pred, color = wave)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
  geom_errorbar(aes(ymin = lo, ymax = hi), alpha = 0.3, width = 0, linewidth = 0.4) +
  geom_point(size = 2, alpha = 0.8) +
  scale_color_manual(values = WAVE_COLORS) +
  facet_wrap(~ state, scales = "free") +
  labs(title = "Observed vs. Predicted Daily Totals",
       subtitle = "Forecast horizons h = 3, 4, 5 days | error bars = 95% PI",
       x = "Observed daily cases", y = "Predicted daily cases",
       color = "Wave") +
  THEME_CV
save_fig(p_scatter, "cv_obs_vs_pred_scatter.png", w = 20/2.54, h = 10/2.54)

# ── B5: Bias time series ──
tcv_daily[, bias := pred - obs]
tcv_daily[, forecast_day := seq_len(.N), by = .(wave, horizon)]

p_bias <- ggplot(tcv_daily, aes(x = forecast_day, y = bias, color = wave)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_line(linewidth = 0.8, alpha = 0.85) +
  geom_point(size = 2) +
  scale_color_manual(values = WAVE_COLORS) +
  facet_grid(state ~ horizon, scales = "free_y",
             labeller = labeller(horizon = function(x) paste0("h = ", x))) +
  labs(title = "Forecast Bias Over Lead Time",
       subtitle = "pred - obs: positive = overprediction",
       x = "Forecast day index", y = "Bias (pred - obs)",
       color = "Wave") +
  THEME_CV
save_fig(p_bias, "cv_bias_timeseries.png", w = 22/2.54, h = 14/2.54)

# ── B6: CRPS by horizon line plot ──
p_crps <- ggplot(tcv_scores, aes(x = factor(horizon), y = CRPS,
                                  color = wave, group = wave)) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  scale_color_manual(values = WAVE_COLORS) +
  facet_wrap(~ state, scales = "free_y") +
  labs(title = "CRPS by Forecast Horizon",
       subtitle = "Lower = better calibrated probabilistic forecast",
       x = "Forecast horizon (days)", y = "CRPS",
       color = "Wave") +
  THEME_CV
save_fig(p_crps, "cv_crps_horizon.png", w = 18/2.54, h = 10/2.54)

# ── B7: Coverage grouped bar chart — NB-corrected ──
# Prefer NB daily coverage where available; fall back to Gaussian
if ("Cov95_NB_daily" %in% names(tcv_scores)) {
  tcv_scores[, Cov95_plot := fifelse(!is.na(Cov95_NB_daily), Cov95_NB_daily, Cov95)]
  cov_subtitle <- "NB predictive distribution | Red line = 95% nominal"
} else {
  tcv_scores[, Cov95_plot := Cov95]
  cov_subtitle <- "Red dashed line = nominal 95% | ST_IV model"
}

p_cov <- ggplot(tcv_scores, aes(x = wave, y = Cov95_plot * 100, fill = factor(horizon))) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6, alpha = 0.85) +
  geom_hline(yintercept = 95, linetype = "dashed", color = "red", linewidth = 0.6) +
  scale_fill_brewer(palette = "Set2") +
  labs(title = "95% Prediction Interval Coverage (Daily Aggregates)",
       subtitle = cov_subtitle,
       x = NULL, y = "Coverage (%)",
       fill = "Horizon") +
  THEME_CV +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))
save_fig(p_cov, "cv_coverage_barplot.png", w = 18/2.54, h = 10/2.54)

# ── B8: Score heatmap (wave × horizon) — includes NB-corrected metrics ──
heat_cols <- intersect(c("R2", "CRPS", "RMSE", "MAE", "Cov95",
                          "Cov95_NB_daily", "daily_bias", "daily_relRMSE"),
                        names(tcv_scores))
heat_dt <- tcv_scores[, c("wave", "horizon", heat_cols), with = FALSE]
heat_long <- melt(heat_dt, id.vars = c("wave", "horizon"),
                  variable.name = "metric", value.name = "value")

# Pretty labels for metrics
metric_labels <- c(R2 = "R\u00B2 (cell)", CRPS = "CRPS (cell)", RMSE = "RMSE (cell)",
                   MAE = "MAE (cell)", Cov95 = "Cov95 (cell)",
                   Cov95_NB_daily = "Cov95 NB\n(daily)", daily_bias = "Bias\n(daily)",
                   daily_relRMSE = "relRMSE\n(daily)")
heat_long[, metric_label := factor(metric_labels[as.character(metric)],
                                    levels = metric_labels[heat_cols])]

p_heatmap <- ggplot(heat_long, aes(x = factor(horizon), y = wave, fill = value)) +
  geom_tile(color = "white", linewidth = 0.8) +
  geom_text(aes(label = sprintf("%.3f", value)), size = 2.2, color = "white") +
  facet_wrap(~ metric_label, scales = "free", nrow = 1) +
  scale_fill_viridis_c(option = "magma") +
  labs(title = "Validation Score Heatmap",
       subtitle = "Wave (rows) \u00D7 Forecast Horizon (columns) | cell-level + daily aggregate metrics",
       x = "Horizon (days)", y = NULL, fill = "Value") +
  THEME_CV +
  theme(axis.text.y = element_text(size = 8),
        strip.text  = element_text(face = "bold", size = 8))
save_fig(p_heatmap, "cv_score_heatmap.png", w = 32/2.54, h = 10/2.54)

} else {
  message("[cv] Section B (temporal forward CV figures) skipped — new schema.")
}


# ═══════════════════════════════════════════════════════════════════════════════
# C) POSHAN-INSPIRED MODEL DIAGNOSTICS
# ═══════════════════════════════════════════════════════════════════════════════

# ── C9: Spatial random effects map (BYM2 posterior mean) ──
spatial_panels <- list()
for (w in WAVES) {
  st <- sub("_.*", "", w)
  mr <- model_list[[w]]
  if (is.null(mr)) next
  st4 <- mr$model_outputs[[best_model_for(w, mr)]]
  if (is.null(st4$spatial_effects)) next

  grid_path <- file.path(BASE, "grids", sprintf("grid_100m_%s.gpkg", st))
  if (!file.exists(grid_path)) next
  grid_sf <- st_read(grid_path, quiet = TRUE)

  n_grids <- mr$n_grids
  se <- st4$spatial_effects[1:n_grids, ]   # first n_grids = total (or structured)
  se_dt <- data.table(grid_id = se$grid_id, spatial_mean = se$mean)
  grid_sf <- merge(grid_sf, se_dt, by = "grid_id")

  cap_lo <- quantile(grid_sf$spatial_mean, 0.02, na.rm = TRUE)
  cap_hi <- quantile(grid_sf$spatial_mean, 0.98, na.rm = TRUE)
  grid_sf$spatial_cap <- pmax(pmin(grid_sf$spatial_mean, cap_hi), cap_lo)

  p <- ggplot(grid_sf) +
    geom_sf(aes(fill = spatial_cap), linewidth = 0, color = NA) +
    scale_fill_viridis_c(option = "viridis", name = "u\u1d62") +
    labs(title = NICE_WAVE[w] %||% w) +
    theme_void(base_size = 8) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 10),
          legend.position = "right",
          legend.key.height = unit(0.4, "cm"),
          legend.key.width = unit(0.3, "cm"))
  spatial_panels[[w]] <- p
}

if (length(spatial_panels) > 0) {
  p_spatial <- wrap_plots(spatial_panels, ncol = 3) +
    plot_annotation(
      title = "Spatial Random Effects — BYM2 Posterior Mean (ST_IV)",
      subtitle = "Positive = higher risk than expected | Negative = lower risk",
      theme = theme(plot.title    = element_text(face = "bold", hjust = 0.5, size = 13),
                    plot.subtitle = element_text(hjust = 0.5, color = "grey40", size = 9)))
  save_fig(p_spatial, "cv_spatial_random_effects.png",
           w = 24/2.54, h = ceiling(length(spatial_panels) / 3) * 10/2.54)
}

# ── C10: BYM2 component decomposition (structured vs unstructured) ──
bym2_panels <- list()
for (w in WAVES) {
  st <- sub("_.*", "", w)
  mr <- model_list[[w]]
  if (is.null(mr)) next
  st4 <- mr$model_outputs[[best_model_for(w, mr)]]
  if (is.null(st4$spatial_effects)) next

  n_grids <- mr$n_grids
  if (nrow(st4$spatial_effects) < 2 * n_grids) next

  grid_path <- file.path(BASE, "grids", sprintf("grid_100m_%s.gpkg", st))
  if (!file.exists(grid_path)) next
  grid_sf <- st_read(grid_path, quiet = TRUE)

  structured   <- st4$spatial_effects[1:n_grids, ]
  unstructured <- st4$spatial_effects[(n_grids + 1):(2 * n_grids), ]
  total_effect <- structured$mean + unstructured$mean

  bym_dt <- data.table(grid_id = structured$grid_id,
                       structured = structured$mean,
                       unstructured = unstructured$mean,
                       total = total_effect)
  grid_sf <- merge(grid_sf, bym_dt, by = "grid_id")

  # Cap for colour scales
  for (col in c("structured", "unstructured", "total")) {
    lo <- quantile(grid_sf[[col]], 0.02, na.rm = TRUE)
    hi <- quantile(grid_sf[[col]], 0.98, na.rm = TRUE)
    grid_sf[[paste0(col, "_cap")]] <- pmax(pmin(grid_sf[[col]], hi), lo)
  }

  # Mixing parameter phi
  hp <- st4$hyperparameters
  phi_row <- hp[grepl("Phi", hp$parameter), ]
  phi_val <- if (nrow(phi_row) > 0) round(phi_row$mean[1], 3) else NA

  p1 <- ggplot(grid_sf) +
    geom_sf(aes(fill = structured_cap), linewidth = 0, color = NA) +
    scale_fill_viridis_c(option = "plasma", name = "Effect") +
    labs(title = "Structured (ICAR)") +
    theme_void(base_size = 7) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 8),
          legend.key.height = unit(0.3, "cm"))

  p2 <- ggplot(grid_sf) +
    geom_sf(aes(fill = unstructured_cap), linewidth = 0, color = NA) +
    scale_fill_viridis_c(option = "viridis", name = "Effect") +
    labs(title = "Unstructured (IID)") +
    theme_void(base_size = 7) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 8),
          legend.key.height = unit(0.3, "cm"))

  p3 <- ggplot(grid_sf) +
    geom_sf(aes(fill = total_cap), linewidth = 0, color = NA) +
    scale_fill_viridis_c(option = "magma", name = "Effect") +
    labs(title = "Total") +
    theme_void(base_size = 7) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 8),
          legend.key.height = unit(0.3, "cm"))

  phi_label <- if (!is.na(phi_val)) bquote(phi == .(phi_val)) else ""
  bym2_panels[[w]] <- (p1 | p2 | p3) +
    plot_annotation(title = paste0(NICE_WAVE[w] %||% w, if (!is.na(phi_val)) sprintf("  (\u03C6 = %.3f)", phi_val) else ""),
                    theme = theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 10)))
}

if (length(bym2_panels) > 0) {
  p_bym2 <- wrap_plots(bym2_panels, ncol = 1) +
    plot_annotation(
      title = "BYM2 Decomposition — Structured (ICAR) vs Unstructured (IID)",
      theme = theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 13)))
  save_fig(p_bym2, "cv_bym2_decomposition.png",
           w = 24/2.54, h = length(bym2_panels) * 8/2.54)
}

# ── C11: Temporal random effect with credible band ──
# Overlay all three waves per state on a shared "Day" axis (matches the
# look of temporal_rr_{ny,ut}_waves.png). The RW1 anchor (idx_temporal == 1)
# is dropped so the smoothing-prior spike does not dominate the y-axis.
.te_wave_lab <- c(`1` = "Wave 1: Initial outbreak (Mar-May 2020)",
                  `2` = "Wave 2: Alpha winter surge (Dec 2020-Feb 2021)",
                  `3` = "Wave 3: Omicron surge (Dec 2021-Feb 2022)")
.te_wave_col <- c(`1` = "#1B9E77", `2` = "#D95F02", `3` = "#7570B3")
.te_state_lab <- c(NY = "New York City", UT = "Utah")

build_te_state_panel <- function(state) {
  rows <- list()
  for (wn in 1:3) {
    w  <- paste0(state, "_wave", wn)
    mr <- model_list[[w]]
    if (is.null(mr)) next
    st4 <- mr$model_outputs[[best_model_for(w, mr)]]
    if (is.null(st4$temporal_effects)) next
    te <- as.data.table(st4$temporal_effects)
    if ("q025" %in% names(te)) {
      setnames(te, c("q025", "q975"), c("lower", "upper"), skip_absent = TRUE)
    } else if ("0.025quant" %in% names(te)) {
      setnames(te, c("0.025quant", "0.975quant"), c("lower", "upper"), skip_absent = TRUE)
    }
    te <- te[idx_temporal > 1]
    te[, day  := idx_temporal - 1L]
    te[, wave := wn]
    rows[[as.character(wn)]] <- te[, .(day, mean, lower, upper, wave)]
  }
  if (length(rows) == 0) return(NULL)
  d <- rbindlist(rows)
  d[, wave_f := factor(wave, levels = 1:3, labels = .te_wave_lab[as.character(1:3)])]

  ggplot(d, aes(x = day, y = mean, colour = wave_f, fill = wave_f,
                ymin = lower, ymax = upper)) +
    geom_ribbon(alpha = 0.18, colour = NA) +
    geom_line(linewidth = 0.85) +
    geom_hline(yintercept = 0, linetype = "dashed",
               colour = "grey30", linewidth = 0.4) +
    scale_colour_manual(values = unname(.te_wave_col), name = NULL) +
    scale_fill_manual(values   = unname(.te_wave_col), name = NULL) +
    scale_x_continuous(breaks = seq(0, 90, by = 15), expand = c(0.01, 0)) +
    labs(x = "Day", y = "Temporal effect",
         title = paste0(.te_state_lab[state],
                        ": temporal random effect across waves")) +
    theme_bw(base_size = 12) +
    theme(legend.position    = "bottom",
          legend.box.spacing = unit(0, "pt"),
          plot.title         = element_text(face = "bold", size = 12),
          panel.grid.minor   = element_line(linewidth = 0.2, colour = "grey92"),
          panel.grid.major   = element_line(linewidth = 0.3, colour = "grey85"))
}

temp_panels <- Filter(Negate(is.null),
                      list(NY = build_te_state_panel("NY"),
                           UT = build_te_state_panel("UT")))

if (length(temp_panels) > 0) {
  p_temporal <- wrap_plots(temp_panels, ncol = 1) +
    plot_annotation(
      title    = "Temporal Random Effects — RW1 Posterior Mean with 95% CI (ST_IV)",
      subtitle = "Day 1 (random-walk anchor) omitted | Positive = elevated transmission",
      theme = theme(plot.title    = element_text(face = "bold", hjust = 0.5, size = 13),
                    plot.subtitle = element_text(hjust = 0.5, color = "grey40", size = 9)))
  save_fig(p_temporal, "cv_temporal_random_effects.png",
           w = 34/2.54, h = length(temp_panels) * 10/2.54)
}

# ── C12: Observed vs fitted temporal traces (top/bottom areas) ──
trace_panels <- list()
for (w in WAVES) {
  mr <- model_list[[w]]
  if (is.null(mr)) next
  st4 <- mr$model_outputs[[best_model_for(w, mr)]]
  if (is.null(st4$fitted_values)) next

  fv <- as.data.table(st4$fitted_values)
  if (!all(c("Y", "fitted_mean", "grid_id", "time_id") %in% names(fv))) next

  # Identify top 5 and bottom 5 grids by total observed cases
  grid_totals <- fv[, .(total_Y = sum(Y, na.rm = TRUE)), by = grid_id]
  active_grids <- grid_totals[total_Y > 0]
  if (nrow(active_grids) < 10) next

  top5 <- active_grids[order(-total_Y)][seq_len(min(5, nrow(active_grids))), grid_id]
  bot5 <- active_grids[order(total_Y)][seq_len(min(5, nrow(active_grids))), grid_id]
  sel_grids <- c(top5, bot5)

  fv_sel <- fv[grid_id %in% sel_grids]
  fv_sel[, group := ifelse(grid_id %in% top5, "Highest", "Lowest")]

  # Aggregate to time series per group
  ts_agg <- fv_sel[, .(obs = sum(Y, na.rm = TRUE),
                         fitted = sum(fitted_mean, na.rm = TRUE)),
                    by = .(time_id, group)]

  ts_long <- melt(ts_agg, id.vars = c("time_id", "group"),
                  measure.vars = c("obs", "fitted"),
                  variable.name = "type", value.name = "cases")

  p <- ggplot(ts_long, aes(x = time_id, y = cases, color = group, linetype = type)) +
    geom_line(linewidth = 0.7, alpha = 0.85) +
    scale_linetype_manual(values = c(obs = "dashed", fitted = "solid"),
                          labels = c(obs = "Observed", fitted = "Fitted")) +
    scale_color_manual(values = c(Highest = "#E41A1C", Lowest = "#377EB8")) +
    labs(title = NICE_WAVE[w] %||% w, x = "Day", y = "Cases",
         color = "Area group", linetype = "Type") +
    theme_minimal(base_size = 9) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 10),
          legend.position = "bottom",
          legend.key.size = unit(0.4, "cm"))

  trace_panels[[w]] <- p
}

if (length(trace_panels) > 0) {
  p_traces <- wrap_plots(trace_panels, ncol = 3) +
    plot_annotation(
      title = "Observed vs. Fitted — Top 5 vs Bottom 5 Areas by Case Volume",
      subtitle = "Dashed = observed | Solid = fitted (ST_IV)",
      theme = theme(plot.title    = element_text(face = "bold", hjust = 0.5, size = 13),
                    plot.subtitle = element_text(hjust = 0.5, color = "grey40", size = 9)))
  save_fig(p_traces, "cv_obs_vs_fitted_traces.png",
           w = 26/2.54, h = ceiling(length(trace_panels) / 3) * 10/2.54)
}

# ── C13: Observed vs fitted spatial facets (peak / trough snapshots) ──
snapshot_panels <- list()
for (w in WAVES) {
  st <- sub("_.*", "", w)
  mr <- model_list[[w]]
  if (is.null(mr)) next
  st4 <- mr$model_outputs[[best_model_for(w, mr)]]
  if (is.null(st4$fitted_values)) next

  grid_path <- file.path(BASE, "grids", sprintf("grid_100m_%s.gpkg", st))
  if (!file.exists(grid_path)) next
  grid_sf <- st_read(grid_path, quiet = TRUE)

  fv <- as.data.table(st4$fitted_values)
  if (!all(c("Y", "fitted_mean", "grid_id", "time_id") %in% names(fv))) next

  # Find peak day and trough day (by total observed)
  day_totals <- fv[, .(total = sum(Y, na.rm = TRUE)), by = time_id]
  peak_day  <- day_totals[which.max(total), time_id]
  # Trough: day with minimum positive total
  trough_day <- day_totals[total > 0][which.min(total), time_id]
  if (is.na(trough_day)) trough_day <- day_totals[which.min(total), time_id]

  for (day_info in list(list(d = peak_day, label = "Peak"),
                         list(d = trough_day, label = "Trough"))) {
    day <- day_info$d
    lab <- day_info$label

    fv_day <- fv[time_id == day, .(grid_id, Y, fitted_mean)]
    grid_day <- merge(grid_sf, fv_day, by = "grid_id")

    grid_long <- grid_day |>
      select(grid_id, Y, fitted_mean) |>
      pivot_longer(cols = c(Y, fitted_mean), names_to = "type", values_to = "cases") |>
      st_as_sf()
    grid_long$type <- factor(grid_long$type,
                              levels = c("Y", "fitted_mean"),
                              labels = c("Observed", "Fitted"))

    cap <- quantile(c(fv_day$Y, fv_day$fitted_mean), 0.98, na.rm = TRUE)
    cap <- max(cap, 1)

    p <- ggplot(grid_long) +
      geom_sf(aes(fill = pmin(cases, cap)), linewidth = 0, color = NA) +
      facet_wrap(~ type) +
      scale_fill_viridis_c(option = "plasma", name = "Cases",
                           labels = label_number(accuracy = 1),
                           limits = c(0, cap)) +
      labs(title = sprintf("%s — %s (day %d)", w, lab, day)) +
      theme_void(base_size = 7) +
      theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 9),
            strip.text = element_text(face = "bold", size = 8),
            legend.position = "right",
            legend.key.height = unit(0.3, "cm"),
            legend.key.width = unit(0.2, "cm"))

    snapshot_panels[[paste0(w, "_", lab)]] <- p
  }
}

if (length(snapshot_panels) > 0) {
  p_snapshots <- wrap_plots(snapshot_panels, ncol = 2) +
    plot_annotation(
      title = "Observed vs. Fitted — Spatial Snapshots at Peak and Trough Days",
      theme = theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 13)))
  save_fig(p_snapshots, "cv_spatial_snapshots.png",
           w = 24/2.54, h = length(snapshot_panels) / 2 * 7/2.54)
}

# ═══════════════════════════════════════════════════════════════════════════════
# PIT histograms — INLA PIT and Randomized PIT (Czado et al. 2009)
# ═══════════════════════════════════════════════════════════════════════════════

# Helper to build a PIT histogram panel
make_pit_panel <- function(pit_vals, title_txt) {
  pit_vals <- pit_vals[is.finite(pit_vals) & pit_vals >= 0 & pit_vals <= 1]
  if (length(pit_vals) < 100) return(NULL)
  dt <- data.table(pit = pit_vals)
  ggplot(dt, aes(x = pit)) +
    geom_histogram(aes(y = after_stat(density)), bins = 20,
                   fill = "steelblue", color = "white", alpha = 0.8) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "red", linewidth = 0.8) +
    labs(title = title_txt, x = "PIT value", y = "Density") +
    theme_minimal(base_size = 9) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 10))
}

# -- INLA PIT (legacy) --
pit_panels_inla <- list()
for (w in WAVES) {
  if (is.null(cpo_list[[w]])) next
  pit_panels_inla[[w]] <- make_pit_panel(cpo_list[[w]]$pit, w)
}
if (length(pit_panels_inla) > 0) {
  p_pit_inla <- wrap_plots(pit_panels_inla, ncol = 3) +
    plot_annotation(
      title = "INLA PIT Histograms — LOO Calibration (ST_IV)",
      subtitle = "Uniform = well calibrated | U-shape = underdispersed | Hump = overdispersed",
      theme = theme(plot.title    = element_text(face = "bold", hjust = 0.5, size = 13),
                    plot.subtitle = element_text(hjust = 0.5, color = "grey40", size = 9)))
  save_fig(p_pit_inla, "cv_pit_all.png",
           w = 22/2.54, h = ceiling(length(pit_panels_inla) / 3) * 8/2.54)
}

# -- Randomized PIT (correct for discrete NB) --
pit_panels_rand <- list()
for (w in WAVES) {
  if (is.null(cpo_list[[w]]) || is.null(cpo_list[[w]]$pit_randomized)) next
  pit_panels_rand[[w]] <- make_pit_panel(cpo_list[[w]]$pit_randomized, w)
}
if (length(pit_panels_rand) > 0) {
  p_pit_rand <- wrap_plots(pit_panels_rand, ncol = 3) +
    plot_annotation(
      title = "Randomized PIT Histograms — Discrete NB Calibration (Czado et al. 2009)",
      subtitle = "Uniform = well calibrated | correct for discrete count data",
      theme = theme(plot.title    = element_text(face = "bold", hjust = 0.5, size = 13),
                    plot.subtitle = element_text(hjust = 0.5, color = "grey40", size = 9)))
  save_fig(p_pit_rand, "cv_pit_randomized.png",
           w = 22/2.54, h = ceiling(length(pit_panels_rand) / 3) * 8/2.54)
}

# -- Combined 2-row comparison panel --
if (length(pit_panels_inla) > 0 && length(pit_panels_rand) > 0) {
  # Build two-row layout: top = INLA, bottom = Randomized
  row_inla <- wrap_plots(pit_panels_inla, ncol = 3) +
    plot_annotation(subtitle = "INLA PIT (continuous approximation)",
                    theme = theme(plot.subtitle = element_text(face = "bold", hjust = 0.5, size = 10)))
  row_rand <- wrap_plots(pit_panels_rand, ncol = 3) +
    plot_annotation(subtitle = "Randomized PIT (discrete NB — Czado et al. 2009)",
                    theme = theme(plot.subtitle = element_text(face = "bold", hjust = 0.5, size = 10)))
  p_pit_compare <- row_inla / row_rand +
    plot_annotation(
      title = "PIT Calibration Comparison",
      theme = theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 13)))
  n_rows <- ceiling(length(pit_panels_inla) / 3) + ceiling(length(pit_panels_rand) / 3)
  save_fig(p_pit_compare, "cv_pit_comparison.png",
           w = 22/2.54, h = n_rows * 8/2.54)
}

# ═══════════════════════════════════════════════════════════════════════════════
# D) IN-SAMPLE CALIBRATION SUMMARY (NB coverage + PIT KS comparison)
# ═══════════════════════════════════════════════════════════════════════════════
if ("Cov95_NB" %in% names(cpo_summary) &&
    "PIT_KS_INLA" %in% names(cpo_summary) &&
    "PIT_KS_randomized" %in% names(cpo_summary)) {

  # D1: NB coverage bar chart
  p_nb_cov <- ggplot(cpo_summary, aes(x = wave, y = Cov95_NB * 100, fill = state)) +
    geom_col(width = 0.65, alpha = 0.85) +
    geom_text(aes(label = sprintf("%.1f%%", Cov95_NB * 100)), vjust = -0.5, size = 3) +
    geom_hline(yintercept = 95, linetype = "dashed", color = "red", linewidth = 0.6) +
    scale_fill_manual(values = STATE_COLORS) +
    labs(title = "In-Sample 95% Coverage — NB Predictive Distribution",
         subtitle = "qnbinom quantiles | Red line = 95% nominal",
         x = NULL, y = "Coverage (%)", fill = "State") +
    THEME_CV +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))
  save_fig(p_nb_cov, "cv_insample_nb_coverage.png", w = 16/2.54, h = 10/2.54)

  # D2: PIT KS comparison (INLA vs Randomized)
  pit_ks_dt <- melt(cpo_summary[, .(wave, state, INLA = PIT_KS_INLA,
                                      Randomized = PIT_KS_randomized)],
                     id.vars = c("wave", "state"),
                     variable.name = "PIT_type", value.name = "KS_D")

  p_ks <- ggplot(pit_ks_dt, aes(x = wave, y = KS_D, fill = PIT_type)) +
    geom_col(position = position_dodge(width = 0.7), width = 0.6, alpha = 0.85) +
    geom_text(aes(label = sprintf("%.3f", KS_D)),
              position = position_dodge(width = 0.7), vjust = -0.5, size = 2.5) +
    scale_fill_manual(values = c(INLA = "#FC8D62", Randomized = "#66C2A5")) +
    labs(title = "PIT Uniformity — KS Distance (lower = better)",
         subtitle = "INLA PIT vs Randomized PIT (Czado et al. 2009)",
         x = NULL, y = "KS Distance (D)", fill = "PIT Method") +
    THEME_CV +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))
  save_fig(p_ks, "cv_pit_ks_comparison.png", w = 18/2.54, h = 10/2.54)
}

# ═══════════════════════════════════════════════════════════════════════════════
# E) MODZCTA Spatial CV Scatter (NYC) — true spatial-CV per-MODZCTA totals
# ═══════════════════════════════════════════════════════════════════════════════
ny_waves_modzcta <- list(
  list(id = "NY_wave1", label = "Wave 1: Mar\u2013May 2020 (Lockdown)"),
  list(id = "NY_wave2", label = "Wave 2: Dec 2020\u2013Feb 2021 (Winter Surge)"),
  list(id = "NY_wave3", label = "Wave 3: Dec 2021\u2013Feb 2022 (Omicron)")
)

scatter_list_modzcta <- list()
for (w in ny_waves_modzcta) {
  # Read spatial CV held-out predictions per cell-day, aggregate to per-MODZCTA totals.
  rds_path <- file.path(BASE, w$id, "validation", "spatial_cv_area_results.rds")
  if (!file.exists(rds_path)) next
  rv <- readRDS(rds_path)
  pp <- as.data.table(rv$predictions)
  if (!all(c("geo_id", "y_true", "mu") %in% names(pp))) next
  cv <- pp[, .(obs_rate = sum(y_true, na.rm = TRUE),
               pred_rate = sum(mu,     na.rm = TRUE)), by = geo_id]
  cv <- cv[is.finite(obs_rate) & is.finite(pred_rate) & obs_rate >= 0]
  if (nrow(cv) < 5) next
  cv[, log_obs  := log10(obs_rate + 1)]
  cv[, log_pred := log10(pred_rate + 1)]
  r2_val   <- round(cor(cv$obs_rate, cv$pred_rate)^2, 3)
  cor_val  <- round(cor(cv$obs_rate, cv$pred_rate), 3)
  cv[, ape := abs(pred_rate - obs_rate) / pmax(obs_rate, 1) * 100]
  mape_val <- round(median(cv$ape), 1)

  scatter_list_modzcta[[w$label]] <- ggplot(cv, aes(x = obs_rate, y = pred_rate)) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey50", linewidth = 0.6) +
    geom_point(color = "#2166AC", alpha = 0.7, size = 2.5) +
    geom_smooth(method = "lm", se = TRUE, color = "#B2182B", linewidth = 0.8, alpha = 0.15) +
    scale_x_log10() + scale_y_log10() +
    annotate("text", x = Inf, y = 0, hjust = 1.1, vjust = -0.5, size = 3.8,
             label = sprintf("R\u00b2 = %.3f\nr = %.3f\nMdAPE = %.1f%%", r2_val, cor_val, mape_val)) +
    labs(title = w$label, x = "Observed cases (MODZCTA total)", y = "Predicted cases (MODZCTA total)") +
    theme_minimal(base_size = 11) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 12),
          legend.position = "none")
}
if (length(scatter_list_modzcta) > 0) {
  p_modzcta <- wrap_plots(scatter_list_modzcta, ncol = 3) +
    plot_annotation(
      title = "Spatial Cross-Validation: Observed vs Predicted MODZCTA Case Rates",
      subtitle = "Out-of-sample predictions from 10-fold spatially blocked CV",
      theme = theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
                    plot.subtitle = element_text(hjust = 0.5, color = "grey40", size = 10)))
  save_fig(p_modzcta, "ny_modzcta_obs_vs_pred_scatter.png", w = 16, h = 5.5)
  file.copy(file.path(FIG_DIR, "ny_modzcta_obs_vs_pred_scatter.png"),
            file.path(LATEX_FIG, "ny_modzcta_obs_vs_pred_scatter.png"), overwrite = TRUE)
}

# ═══════════════════════════════════════════════════════════════════════════════
# F) NYC Temporal CV Scatter + Trajectory
# ═══════════════════════════════════════════════════════════════════════════════
ny_waves_tcv <- list(
  list(id = "NY_wave1", label = "Wave 1: Mar\u2013May 2020 (Lockdown)", start = as.Date("2020-03-01")),
  list(id = "NY_wave2", label = "Wave 2: Dec 2020\u2013Feb 2021 (Winter Surge)", start = as.Date("2020-12-01")),
  list(id = "NY_wave3", label = "Wave 3: Dec 2021\u2013Feb 2022 (Omicron)", start = as.Date("2021-12-01"))
)
tcv_scatter_list <- list(); tcv_ts_list <- list()
for (w in ny_waves_tcv) {
  raw_path <- file.path(BASE, w$id, "validation", "temporal_cv_daily_raw_h3.csv")
  if (!file.exists(raw_path)) next
  raw <- fread(raw_path)
  raw[, date := w$start + time_id - 1L]
  r2_val  <- round(cor(raw$obs, raw$pred)^2, 3)
  cor_val <- round(cor(raw$obs, raw$pred), 3)
  rmse_val <- round(sqrt(mean((raw$obs - raw$pred)^2)), 0)

  tcv_scatter_list[[w$label]] <- ggplot(raw, aes(x = obs, y = pred)) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey50", linewidth = 0.6) +
    geom_errorbar(aes(ymin = lo, ymax = hi), width = 0, alpha = 0.3, color = "steelblue") +
    geom_point(aes(color = factor(origin)), size = 3, alpha = 0.8) +
    geom_smooth(method = "lm", se = TRUE, color = "#B2182B", linewidth = 0.8, alpha = 0.15) +
    annotate("text", x = Inf, y = -Inf, hjust = 1.1, vjust = -0.5, size = 3.8,
             label = sprintf("R\u00b2 = %.3f\nr = %.3f\nRMSE = %d", r2_val, cor_val, rmse_val)) +
    labs(title = w$label, x = "Observed Daily Cases", y = "Predicted Daily Cases") +
    theme_minimal(base_size = 11) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 11), legend.position = "none")

  daily_path <- file.path(BASE, w$id, "validation", "temporal_cv_daily_h3.csv")
  if (!file.exists(daily_path)) next
  daily <- fread(daily_path)
  daily[, date := w$start + time_id - 1L]
  tcv_ts_list[[w$label]] <- ggplot(daily, aes(x = date)) +
    geom_ribbon(aes(ymin = lo, ymax = hi), fill = "steelblue", alpha = 0.25) +
    geom_line(aes(y = pred), color = "steelblue", linewidth = 1) +
    geom_point(aes(y = obs), color = "#B2182B", size = 2, alpha = 0.8) +
    labs(title = w$label, x = "Date", y = "Daily Cases") +
    scale_x_date(date_labels = "%b %d") +
    theme_minimal(base_size = 11) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 11))
}
if (length(tcv_scatter_list) > 0) {
  p_tcv_scatter <- wrap_plots(tcv_scatter_list, ncol = 3) +
    plot_annotation(title = "Temporal CV: Observed vs Predicted Daily Cases (NYC)",
                    subtitle = "Out-of-sample 3-day-ahead predictions from expanding-window forward CV",
                    theme = theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
                                  plot.subtitle = element_text(hjust = 0.5, color = "grey40", size = 10)))
  save_fig(p_tcv_scatter, "ny_temporal_cv_obs_vs_pred_scatter.png", w = 16, h = 5.5)
  file.copy(file.path(FIG_DIR, "ny_temporal_cv_obs_vs_pred_scatter.png"),
            file.path(LATEX_FIG, "ny_temporal_cv_obs_vs_pred_scatter.png"), overwrite = TRUE)
}
if (length(tcv_ts_list) > 0) {
  p_tcv_ts <- wrap_plots(tcv_ts_list, ncol = 3) +
    plot_annotation(title = "Temporal CV: Forecast Trajectories \u2014 NYC (3-day horizon)",
                    subtitle = "Blue ribbon = 95% PI; red dots = observed daily total cases",
                    theme = theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
                                  plot.subtitle = element_text(hjust = 0.5, color = "grey40", size = 10)))
  save_fig(p_tcv_ts, "ny_temporal_cv_trajectory.png", w = 16, h = 5)
  file.copy(file.path(FIG_DIR, "ny_temporal_cv_trajectory.png"),
            file.path(LATEX_FIG, "ny_temporal_cv_trajectory.png"), overwrite = TRUE)
}

# ═══════════════════════════════════════════════════════════════════════════════
# G) Utah Spatial + Temporal CV Figures
# ═══════════════════════════════════════════════════════════════════════════════
ut_waves <- list(
  list(id = "UT_wave1", label = "Wave 1: Mar\u2013May 2020 (Lockdown)", start = as.Date("2020-03-01")),
  list(id = "UT_wave2", label = "Wave 2: Dec 2020\u2013Feb 2021 (Winter Surge)", start = as.Date("2020-12-01")),
  list(id = "UT_wave3", label = "Wave 3: Dec 2021\u2013Feb 2022 (Omicron)", start = as.Date("2021-12-01"))
)

# Utah spatial scatter
ut_sp_scatter <- list()
for (w in ut_waves) {
  cv_path <- file.path(BASE, w$id, "validation", "spatial_cv_area_predictions.csv")
  if (!file.exists(cv_path)) next
  cv <- fread(cv_path)
  # New schema: y_true / mu (per-cell, per-time); aggregate to (geo_id, fold) totals
  if (!all(c("y_true", "mu") %in% names(cv))) next
  cv <- cv[, .(obs_rate = sum(y_true, na.rm = TRUE),
               pred_rate = sum(mu,     na.rm = TRUE)),
           by = .(geo_id, fold)]
  cv <- cv[is.finite(obs_rate) & is.finite(pred_rate) & obs_rate + pred_rate > 0]
  if (nrow(cv) < 5) next
  r2_val   <- round(cor(cv$obs_rate, cv$pred_rate)^2, 3)
  cor_val  <- round(cor(cv$obs_rate, cv$pred_rate), 3)
  cv[, ape := abs(pred_rate - obs_rate) / pmax(obs_rate, 1) * 100]
  mape_val <- round(median(cv$ape), 1)

  ut_sp_scatter[[w$label]] <- ggplot(cv, aes(x = obs_rate, y = pred_rate)) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey50", linewidth = 0.6) +
    geom_point(aes(color = factor(fold)), alpha = 0.8, size = 3) +
    geom_smooth(method = "lm", se = TRUE, color = "#B2182B", linewidth = 0.8, alpha = 0.15) +
    scale_x_log10() + scale_y_log10() +
    annotate("text", x = Inf, y = 0, hjust = 1.1, vjust = -0.5, size = 3.8,
             label = sprintf("R\u00b2 = %.3f\nr = %.3f\nMdAPE = %.1f%%", r2_val, cor_val, mape_val)) +
    labs(title = w$label, x = "Observed cases (county total)", y = "Predicted cases (county total)") +
    theme_minimal(base_size = 11) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 11), legend.position = "none")
}
if (length(ut_sp_scatter) > 0) {
  p_ut_sp <- wrap_plots(ut_sp_scatter, ncol = 3) +
    plot_annotation(title = "Spatial CV: Observed vs Predicted County Case Rates (Utah)",
                    subtitle = "Out-of-sample predictions from leave-one-county-out CV",
                    theme = theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
                                  plot.subtitle = element_text(hjust = 0.5, color = "grey40", size = 10)))
  save_fig(p_ut_sp, "ut_spatial_cv_scatter.png", w = 16, h = 5.5)
  file.copy(file.path(FIG_DIR, "ut_spatial_cv_scatter.png"),
            file.path(LATEX_FIG, "ut_spatial_cv_scatter.png"), overwrite = TRUE)
}

# Utah temporal scatter + trajectory
ut_tc_scatter <- list(); ut_tc_ts <- list()
for (w in ut_waves) {
  raw_path <- file.path(BASE, w$id, "validation", "temporal_cv_daily_raw_h3.csv")
  if (!file.exists(raw_path)) next
  raw <- fread(raw_path)
  raw[, date := w$start + time_id - 1L]
  r2_val  <- round(cor(raw$obs, raw$pred)^2, 3)
  cor_val <- round(cor(raw$obs, raw$pred), 3)
  rmse_val <- round(sqrt(mean((raw$obs - raw$pred)^2)), 0)

  ut_tc_scatter[[w$label]] <- ggplot(raw, aes(x = obs, y = pred)) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey50", linewidth = 0.6) +
    geom_errorbar(aes(ymin = lo, ymax = hi), width = 0, alpha = 0.3, color = "steelblue") +
    geom_point(aes(color = factor(origin)), size = 3, alpha = 0.8) +
    geom_smooth(method = "lm", se = TRUE, color = "#B2182B", linewidth = 0.8, alpha = 0.15) +
    annotate("text", x = Inf, y = -Inf, hjust = 1.1, vjust = -0.5, size = 3.8,
             label = sprintf("R\u00b2 = %.3f\nr = %.3f\nRMSE = %d", r2_val, cor_val, rmse_val)) +
    labs(title = w$label, x = "Observed Daily Cases", y = "Predicted Daily Cases") +
    theme_minimal(base_size = 11) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 11), legend.position = "none")

  daily_path <- file.path(BASE, w$id, "validation", "temporal_cv_daily_h3.csv")
  if (!file.exists(daily_path)) next
  daily <- fread(daily_path)
  daily[, date := w$start + time_id - 1L]
  ut_tc_ts[[w$label]] <- ggplot(daily, aes(x = date)) +
    geom_ribbon(aes(ymin = lo, ymax = hi), fill = "steelblue", alpha = 0.25) +
    geom_line(aes(y = pred), color = "steelblue", linewidth = 1) +
    geom_point(aes(y = obs), color = "#B2182B", size = 2, alpha = 0.8) +
    labs(title = w$label, x = "Date", y = "Daily Cases") +
    scale_x_date(date_labels = "%b %d") +
    theme_minimal(base_size = 11) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 11))
}
if (length(ut_tc_scatter) > 0) {
  p_ut_tc <- wrap_plots(ut_tc_scatter, ncol = 3) +
    plot_annotation(title = "Temporal CV: Observed vs Predicted Daily Cases (Utah)",
                    subtitle = "Out-of-sample 3-day-ahead predictions from expanding-window forward CV",
                    theme = theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
                                  plot.subtitle = element_text(hjust = 0.5, color = "grey40", size = 10)))
  save_fig(p_ut_tc, "ut_temporal_cv_scatter.png", w = 16, h = 5.5)
  file.copy(file.path(FIG_DIR, "ut_temporal_cv_scatter.png"),
            file.path(LATEX_FIG, "ut_temporal_cv_scatter.png"), overwrite = TRUE)
}
if (length(ut_tc_ts) > 0) {
  p_ut_traj <- wrap_plots(ut_tc_ts, ncol = 3) +
    plot_annotation(title = "Temporal CV: Forecast Trajectories \u2014 Utah (3-day horizon)",
                    subtitle = "Blue ribbon = 95% PI; red dots = observed daily total cases",
                    theme = theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
                                  plot.subtitle = element_text(hjust = 0.5, color = "grey40", size = 10)))
  save_fig(p_ut_traj, "ut_temporal_cv_trajectory.png", w = 16, h = 5)
  file.copy(file.path(FIG_DIR, "ut_temporal_cv_trajectory.png"),
            file.path(LATEX_FIG, "ut_temporal_cv_trajectory.png"), overwrite = TRUE)
}

