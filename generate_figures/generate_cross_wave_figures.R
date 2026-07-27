#!/usr/bin/env Rscript
# generate_cross_wave_figures.R
# Consolidated cross-wave figures (NY×3 + UT×3). Uses facet layouts and patchwork.

suppressPackageStartupMessages({
  library(data.table)
  library(sf)
  library(ggplot2)
  library(patchwork)
  library(viridis)
  library(ggspatial)
  library(scales)
  library(dplyr)
})

Sys.unsetenv("MallocStackLogging")
source("helpers/_utils.R")
load_pipeline()

LATEX_FIG <- file.path("thesis_latex", "figures")
dir.create(LATEX_FIG, recursive = TRUE, showWarnings = FALSE)

wave_ids   <- c("wave1", "wave2", "wave3")
states     <- c("NY", "UT")
wave_titles <- c(wave1 = "Wave 1: Mar\u2013May 2020",
                 wave2 = "Wave 2: Dec 2020\u2013Feb 2021",
                 wave3 = "Wave 3: Dec 2021\u2013Feb 2022")

# ══════════════════════════════════════════════════════════════════
# 1. Temporal RR — ALL 6 CONFIGS IN ONE FIGURE (facet_grid)
# ══════════════════════════════════════════════════════════════════

# despike_temporal() is defined in _utils.R; gated by USE_DESPIKE (FALSE
# since switching from monthly to weekly Advan inputs - May 2026).

trr_all <- rbindlist(lapply(states, function(st) {
  rbindlist(lapply(wave_ids, function(w) {
    wm <- CONFIG$waves[[w]]$months
    start_date   <- as.Date(paste0(wm[1], "-01"))
    month_breaks <- as.Date(paste0(wm[-1], "-01"))

    tp <- fread(file.path(CONFIG$output_dir, paste0(st, "_", w),
                          "tables", "unified_temporal_effects.csv"))
    tp[, date := start_date + idx_temporal - 1L]
    if (isTRUE(USE_DESPIKE) && length(month_breaks) > 0) {
      for (cc in c("mean", "q025", "q975"))
        tp[, (cc) := despike_temporal(tp, month_breaks, cc)]
    }
    tp[, `:=`(RR = exp(mean), rr_lo = exp(q025), rr_hi = exp(q975),
              state = st, wave = w)]
    tp[, .(date, RR, rr_lo, rr_hi, state, wave)]
  }))
}))

trr_all[, state := factor(state, levels = c("NY", "UT"),
                          labels = c("NYC", "Utah"))]
trr_all[, wave  := factor(wave, levels = wave_ids, labels = names(wave_titles))]

p_trr <- ggplot(trr_all, aes(x = date, y = RR)) +
  geom_ribbon(aes(ymin = rr_lo, ymax = rr_hi), fill = "#B2182B", alpha = 0.15) +
  geom_line(linewidth = 0.7, color = "#B2182B") +
  geom_point(size = 0.8, color = "#B2182B") +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey50") +
  facet_grid(state ~ wave, scales = "free") +
  scale_x_date(date_labels = "%b %d", date_breaks = "15 days") +
  labs(x = "Date", y = "Temporal Relative Risk",
       title = "Temporal Relative Risk Across Study Configurations (ST_IV)",
       subtitle = "Posterior mean with 95% credible interval") +
  theme_minimal(base_size = 10) +
  theme(plot.title    = element_text(face = "bold", hjust = 0.5, size = 13),
        plot.subtitle = element_text(hjust = 0.5, color = "grey40"),
        strip.text    = element_text(face = "bold", size = 10),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 7))

ggsave(file.path(LATEX_FIG, "temporal_rr_all.png"), p_trr,
       width = 14, height = 7, dpi = 300, bg = "white")

# ══════════════════════════════════════════════════════════════════
# 2. Fixed Effects RR — ALL 6 CONFIGS IN ONE FIGURE (facet_grid)
# ══════════════════════════════════════════════════════════════════

fe_all <- rbindlist(lapply(states, function(st) {
  rbindlist(lapply(wave_ids, function(w) {
    fe <- fread(file.path(CONFIG$output_dir, paste0(st, "_", w),
                          "tables", "unified_fixed_effects.csv"))
    fe <- fe[covariate != "(Intercept)"]
    fe[, `:=`(state = st, wave = w)]
    fe[, se := (log(RR_upper) - log(RR)) / qnorm(0.975)]
    fe[, RR_66lo := exp(log(RR) + qnorm(0.17) * se)]
    fe[, RR_66hi := exp(log(RR) + qnorm(0.83) * se)]
    fe
  }))
}))

fe_all[, state := factor(state, levels = c("NY", "UT"),
                         labels = c("NYC", "Utah"))]
fe_all[, wave  := factor(wave, levels = wave_ids, labels = names(wave_titles))]

# Rank covariates by average RR across configs for consistent ordering
cov_order <- fe_all[, .(avg_RR = mean(RR, na.rm = TRUE)), by = covariate]
setorder(cov_order, avg_RR)
fe_all[, covariate := factor(covariate, levels = cov_order$covariate)]

p_fe <- ggplot(fe_all, aes(x = RR, y = covariate, color = significant)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey50") +
  geom_linerange(aes(xmin = RR_lower, xmax = RR_upper),
                 linewidth = 0.3, alpha = 0.4) +
  geom_linerange(aes(xmin = RR_66lo, xmax = RR_66hi), linewidth = 1.2) +
  geom_point(size = 1.5) +
  scale_color_manual(values = c("TRUE" = "#B2182B", "FALSE" = "#999999"),
                     guide = "none") +
  facet_grid(state ~ wave, scales = "free_y", space = "free_y") +
  labs(x = "Relative Risk (exp(\u03b2))", y = NULL,
       title = "Fixed Effects Across Study Configurations (ST_IV)",
       subtitle = "Thick bar = 66% CI | Thin whisker = 95% CI | Red = significant") +
  theme_minimal(base_size = 9) +
  theme(plot.title    = element_text(face = "bold", hjust = 0.5, size = 13),
        plot.subtitle = element_text(hjust = 0.5, color = "grey40"),
        strip.text    = element_text(face = "bold", size = 9),
        axis.text.y   = element_text(size = 6),
        panel.spacing.y = unit(0.8, "lines"))

ggsave(file.path(LATEX_FIG, "fixed_effects_rr_all.png"), p_fe,
       width = 14, height = 12, dpi = 300, bg = "white")

# ══════════════════════════════════════════════════════════════════
# 3. WAIC Heatmap + Model Transition (side by side)
# ══════════════════════════════════════════════════════════════════

wave_dlabs <- c("NYC W1", "NYC W2", "NYC W3", "UT W1", "UT W2", "UT W3")
wave_dirs  <- c("NY_wave1", "NY_wave2", "NY_wave3",
                "UT_wave1", "UT_wave2", "UT_wave3")

waic_all <- rbindlist(lapply(seq_along(wave_dirs), function(i) {
  dt <- fread(file.path(CONFIG$output_dir, wave_dirs[i], "tables",
                        "model_comparison_unified.csv"))
  dt <- dt[Status == "OK"]
  dt[, dWAIC := WAIC - min(WAIC, na.rm = TRUE)]
  dt[, `:=`(wave = wave_dlabs[i], is_best = dWAIC == 0)]
  dt[, .(Model, WAIC, dWAIC, wave, is_best)]
}))

waic_all[, dWAIC := WAIC - min(WAIC, na.rm = TRUE), by = wave]
# WAIC-best per wave comes from the data; no manual override.
# Canonical: outputs/spatiotemporal_risk/<wave>/models/best_model.txt
# (NY_wave1: NLMob_ST_IV; NY_wave2/3: LinMob_ST_IV; UT_*: NLMob_ST_IV).

model_order <- c("Base", "Temp", "LinMob", "NLMob",
                 "LinMob_ST_IV", "NLMob_ST_IV")
model_labels <- display_model(model_order)
waic_all[, Model := factor(Model, levels = model_order, labels = model_labels)]
waic_all[, wave  := factor(wave, levels = wave_dlabs)]
waic_all[, label := ifelse(is_best, sprintf("%.0f*", dWAIC),
                           sprintf("%.0f", dWAIC))]

p_heat <- ggplot(waic_all, aes(x = wave, y = Model, fill = dWAIC)) +
  geom_tile(color = "white", linewidth = 0.8) +
  geom_text(aes(label = label, 
                color = ifelse(dWAIC > 10000, "white", "black")),  # threshold adjusted for log scale
            size = 12, fontface = "bold") +
  scale_fill_viridis_c(option = "magma", direction = -1,
                       trans = scales::log1p_trans(),              # log(1+x) handles zero and huge range
                       name = expression(Delta * WAIC),
                       breaks = c(0, 10, 100, 1000, 10000, 100000)) +
  scale_color_identity() +
  scale_y_discrete(limits = rev(model_labels)) +
  labs(title = expression(bold("a) " * Delta * "WAIC Heatmap")),
       x = NULL, y = NULL) +
  theme_minimal(base_size = 10) +
  theme(plot.title = element_text(hjust = 0.5, size = 11),
        axis.text  = element_text(size = 9),
        panel.grid = element_blank(),
        legend.position = "right")

trans_dt <- copy(waic_all)
trans_dt[, rank := frank(dWAIC, ties.method = "min"), by = wave]
best_dt  <- trans_dt[is_best == TRUE]

p_rank <- ggplot(trans_dt, aes(x = wave, y = rank,
                               group = Model, color = Model)) +
  geom_line(linewidth = 0.6, alpha = 0.4) +
  geom_point(size = 2, alpha = 0.4) +
  geom_line(data = best_dt, aes(x = wave, y = rank, group = 1),
            color = "black", linewidth = 1.2, linetype = "dashed",
            inherit.aes = FALSE) +
  geom_point(data = best_dt, aes(x = wave, y = rank),
             color = "black", size = 4, shape = 18, inherit.aes = FALSE) +
  geom_label(data = best_dt, aes(x = wave, y = rank, label = Model),
             color = "black", size = 2.5, fontface = "bold",
             nudge_y = -0.4, inherit.aes = FALSE, label.size = 0.3) +
  scale_y_reverse(breaks = 1:6, name = "WAIC Rank (1 = best)") +
  scale_color_brewer(palette = "Dark2") +
  labs(title = "b) Model Rank Transitions",
       x = NULL) +
  theme_minimal(base_size = 10) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 11),
        axis.text  = element_text(size = 9),
        legend.position = "right",
        panel.grid.minor = element_blank())

p_model_comp <- p_heat + p_rank +
  plot_layout(widths = c(1, 1)) +
  plot_annotation(
    title = "Model Comparison Across All Study Configurations",
    theme = theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 13)))

p_model_comp

ggsave(file.path(LATEX_FIG, "model_comparison_combined.png"), p_model_comp,
       width = 16, height = 6, dpi = 300, bg = "white")

# ══════════════════════════════════════════════════════════════════
# 4. Spatial RR — Cross-Wave (patchwork, different geometries)
# ══════════════════════════════════════════════════════════════════

grid_sf_ny <- st_read(file.path(CONFIG$grids_dir, "grid_100m_NY.gpkg"), quiet = TRUE)
grid_sf_ut <- st_read(file.path(CONFIG$grids_dir, "grid_100m_UT.gpkg"), quiet = TRUE)
ut_counties_sf <- st_read(CONFIG$ut_counties_path, quiet = TRUE) |>
  st_transform(st_crs(grid_sf_ut))

make_sp_panel <- function(st, w, grid_sf) {
  sp <- fread(file.path(CONFIG$output_dir, paste0(st, "_", w),
                        "tables", "unified_spatial_effects.csv"))
  n_areas <- nrow(sp) %/% 2L
  sp_s <- sp[idx_spatial <= n_areas]; setorder(sp_s, idx_spatial)
  sp_u <- sp[idx_spatial >  n_areas]; setorder(sp_u, idx_spatial)
  total_eff <- data.table(grid_id = sp_s$grid_id,
                          total_RR = exp(sp_s$mean + sp_u$mean))
  te_sf <- merge(grid_sf, total_eff, by = "grid_id", all.x = FALSE)
  q99 <- quantile(te_sf$total_RR, 0.99, na.rm = TRUE)
  q01 <- quantile(te_sf$total_RR, 0.01, na.rm = TRUE)

  p <- ggplot(te_sf) +
    geom_sf(aes(fill = pmin(pmax(total_RR, q01), q99)),
            color = NA, linewidth = 0) +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                         midpoint = 1, name = "RR",
                         limits = c(q01, q99)) +
    labs(title = wave_titles[w]) +
    theme_void(base_size = 8) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 9))

  if (st == "UT") {
    cents <- suppressWarnings(st_centroid(te_sf))
    p <- p +
      geom_sf(data = ut_counties_sf, fill = NA, color = "grey60",
              linewidth = 0.3, inherit.aes = FALSE) +
      geom_sf_text(data = cents, aes(label = round(total_RR, 2)),
                   size = 1.8, fontface = "bold", inherit.aes = FALSE)
  }
  p
}

# NYC row
sp_ny <- lapply(wave_ids, function(w) make_sp_panel("NY", w, grid_sf_ny))
# UT row
sp_ut <- lapply(wave_ids, function(w) make_sp_panel("UT", w, grid_sf_ut))

ny_row <- (sp_ny[[1]] + sp_ny[[2]] + sp_ny[[3]]) +
  plot_layout(ncol = 3, guides = "collect") &
  theme(legend.position = "bottom")

ut_row <- (sp_ut[[1]] + sp_ut[[2]] + sp_ut[[3]]) +
  plot_layout(ncol = 3, guides = "collect") &
  theme(legend.position = "bottom")

# Stack NYC on top of UT
p_spatial <- (wrap_elements(ny_row) / wrap_elements(ut_row)) +
  plot_annotation(
    title = "Spatial Relative Risk (BYM2) — All Study Configurations",
    subtitle = "Top: NYC (26,700 grid cells)  |  Bottom: Utah (29 counties)",
    theme = theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 13),
      plot.subtitle = element_text(hjust = 0.5, color = "grey40")))

ggsave(file.path(LATEX_FIG, "spatial_rr_all.png"), p_spatial,
       width = 16, height = 14, dpi = 300, bg = "white")

rm(grid_sf_ny, grid_sf_ut, ut_counties_sf); gc()

# ══════════════════════════════════════════════════════════════════
# 5. Observed vs Predicted Daily Cases — ALL CONFIGS (facet by state)
#    Following Poshan's comparison_plot_temporal pattern
# ══════════════════════════════════════════════════════════════════

wave_starts <- c(wave1 = as.Date("2020-03-01"),
                 wave2 = as.Date("2020-12-01"),
                 wave3 = as.Date("2021-12-01"))

daily_all <- rbindlist(lapply(states, function(st) {
  rbindlist(lapply(wave_ids, function(w) {
    f <- file.path(CONFIG$output_dir, paste0(st, "_", w),
                   "tables", "unified_fitted_values.csv")
    dt <- fread(f, select = c("time_id", "Y", "fitted_mean", "fitted_sd"))

    # Aggregate to daily city-wide totals
    daily <- dt[, .(observed  = sum(Y, na.rm = TRUE),
                     predicted = sum(fitted_mean, na.rm = TRUE),
                     agg_sd    = sqrt(sum(fitted_sd^2, na.rm = TRUE))),
                by = time_id]
    daily[, date := wave_starts[w] + time_id - 1L]
    daily[, `:=`(
      ci95_lo = pmax(0, predicted + qnorm(0.025) * agg_sd),
      ci95_hi = predicted + qnorm(0.975) * agg_sd
    )]
    daily[, state := fifelse(st == "NY", "New York City", "Utah")]
    daily[, wave  := wave_titles[w]]
    daily
  }))
}))

# Order wave factor
daily_all[, wave := factor(wave, levels = wave_titles)]

# Wave color palette
wave_colors <- c(
  "Wave 1: Mar\u2013May 2020"         = "#E64B35",
  "Wave 2: Dec 2020\u2013Feb 2021"    = "#4DBBD5",
  "Wave 3: Dec 2021\u2013Feb 2022"    = "#00A087"
)

p_obs_pred <- ggplot(daily_all, aes(x = date)) +
  # 95% CI ribbon for predicted
  geom_ribbon(aes(ymin = ci95_lo, ymax = ci95_hi, fill = wave),
              alpha = 0.2) +
  # Predicted (solid)
  geom_line(aes(y = predicted, color = wave, linetype = "Fitted"),
            linewidth = 0.7) +
  # Observed (dashed)
  geom_line(aes(y = observed, color = wave, linetype = "Observed"),
            linewidth = 0.5) +
  # Aesthetics
  scale_color_manual(values = wave_colors, name = "Wave") +
  scale_fill_manual(values = wave_colors, guide = "none") +
  scale_linetype_manual(values = c("Fitted" = "solid", "Observed" = "dashed"),
                        name = "Type") +
  scale_x_date(date_breaks = "2 weeks", date_labels = "%b %d") +
  scale_y_continuous(labels = label_comma()) +
  facet_wrap(~ state, ncol = 1, scales = "free") +
  labs(title = "Observed vs Fitted Daily Case Counts — All Study Configurations",
       subtitle = "Solid = model fitted values (best model, ST_IV) | Dashed = observed counts | Shaded = 95% CI",
       x = "Date", y = "Daily Cases") +
  theme_minimal(base_size = 11) +
  theme(
    plot.title    = element_text(face = "bold", hjust = 0.5, size = 13),
    plot.subtitle = element_text(hjust = 0.5, color = "grey40", size = 9),
    strip.text    = element_text(face = "bold", size = 12),
    legend.position = "bottom",
    legend.box    = "horizontal",
    axis.text.x   = element_text(angle = 45, hjust = 1, size = 8),
    panel.grid.minor = element_blank(),
    panel.spacing = unit(1, "lines")
  ) +
  guides(color = guide_legend(order = 1), linetype = guide_legend(order = 2))

ggsave(file.path(LATEX_FIG, "obs_vs_pred_all.png"), p_obs_pred,
       width = 14, height = 8, dpi = 300, bg = "white")

# ══════════════════════════════════════════════════════════════════
# 6. Observed vs Fitted — Location-Level Spaghetti (facet_grid)
#    Inspired by 12_visualization.Rmd §modzcta_temporal_rates
# ══════════════════════════════════════════════════════════════════

state_labels <- c(NY = "New York City (MODZCTAs)",
                  UT = "Utah (Counties)")

loc_all <- rbindlist(lapply(states, function(st) {
  rbindlist(lapply(wave_ids, function(w) {
    base <- file.path(CONFIG$output_dir, paste0(st, "_", w))
    ft <- fread(file.path(base, "tables", "unified_fitted_values.csv"),
                select = c("grid_id", "time_id", "Y", "fitted_mean"))
    md <- readRDS(file.path(base, "models", "model_data_gt.rds"))
    gt <- as.data.table(md$gt_dt)
    gm <- unique(gt[, .(grid_id, geo_id)])
    gm <- gm[geo_id != "unknown" & geo_id != "99999"]
    ft <- merge(ft, gm, by = "grid_id")

    loc_day <- ft[, .(obs  = sum(Y, na.rm = TRUE),
                       pred = sum(fitted_mean, na.rm = TRUE)),
                  by = .(geo_id, time_id)]
    loc_day[, date  := wave_starts[w] + time_id - 1L]
    loc_day[, state := state_labels[st]]
    loc_day[, wave  := wave_titles[w]]
    loc_day
  }))
}))

loc_all[, wave := factor(wave, levels = wave_titles)]

agg_all2 <- loc_all[, .(obs = sum(obs), pred = sum(pred)),
                     by = .(state, wave, date)]

loc_long <- rbind(
  loc_all[, .(geo_id, state, wave, date, cases = obs,  type = "Observed")],
  loc_all[, .(geo_id, state, wave, date, cases = pred, type = "Fitted")]
)
agg_long2 <- rbind(
  agg_all2[, .(state, wave, date, cases = obs,  type = "Observed")],
  agg_all2[, .(state, wave, date, cases = pred, type = "Fitted")]
)

p_spag <- ggplot() +
  geom_line(data = loc_long,
            aes(x = date, y = cases, group = interaction(geo_id, type),
                linetype = type),
            alpha = 0.12, linewidth = 0.25, color = "grey50") +
  geom_line(data = agg_long2,
            aes(x = date, y = cases, linetype = type, color = type),
            linewidth = 1.2) +
  facet_grid(state ~ wave, scales = "free") +
  scale_linetype_manual(values = c("Fitted" = "solid", "Observed" = "dashed"),
                        name = "Type") +
  scale_color_manual(values = c("Fitted" = "#E64B35", "Observed" = "#2166AC"),
                     name = "Type") +
  scale_x_date(date_breaks = "2 weeks", date_labels = "%b %d") +
  scale_y_continuous(labels = label_comma()) +
  labs(title = "Observed vs Fitted Daily Case Counts \u2014 All Study Configurations",
       subtitle = "Grey spaghetti = individual locations (MODZCTAs / counties)  |  Bold = aggregate  |  Red solid = fitted, Blue dashed = observed",
       x = "Date", y = "Daily Cases") +
  theme_minimal(base_size = 11) +
  theme(
    plot.title    = element_text(face = "bold", hjust = 0.5, size = 13),
    plot.subtitle = element_text(hjust = 0.5, color = "grey40", size = 9),
    strip.text    = element_text(face = "bold", size = 10),
    legend.position = "bottom",
    axis.text.x   = element_text(angle = 45, hjust = 1, size = 7),
    panel.grid.minor = element_blank(),
    panel.spacing = unit(0.8, "lines")
  )

ggsave(file.path(LATEX_FIG, "obs_vs_fitted_spaghetti.png"), p_spag,
       width = 16, height = 9, dpi = 300, bg = "white")
