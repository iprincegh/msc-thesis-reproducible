#!/usr/bin/env Rscript
# generate_nonlinearity_diagnostics.R
# Predicted vs Observed diagnostics (LOESS) using out-of-sample predictions.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

BASE <- "outputs/spatiotemporal_risk"
OUT  <- "thesis_latex/figures"

# ── Theme (matches publication figures) ──────────────────────────────────────
THEME_PUB <- theme_minimal(base_size = 10) +
  theme(
    plot.title    = element_text(face = "bold", hjust = 0.5, size = 11),
    plot.subtitle = element_text(hjust = 0.5, color = "grey40", size = 8),
    strip.text    = element_text(face = "bold", size = 9),
    panel.grid.minor = element_blank(),
    legend.position  = "bottom"
  )

WAVE_COLS <- c(
  "Wave 1" = "#2166AC", "Wave 2" = "#4393C3", "Wave 3" = "#92C5DE",  # NYC
  "UT Wave 1" = "#B2182B", "UT Wave 2" = "#D6604D", "UT Wave 3" = "#F4A582"
)

# ── Helper: build predicted-vs-observed panel ────────────────────────────────
make_pvo_panel <- function(dt, obs_col, pred_col, sd_col = NULL,
                           title, col_pal, facet_col = "wave_label") {
  dt <- copy(dt)
  setnames(dt, c(obs_col, pred_col), c("obs", "pred"), skip_absent = TRUE)
  if (!is.null(sd_col) && sd_col %in% names(dt)) {
    setnames(dt, sd_col, "sd_pred", skip_absent = TRUE)
    dt[, `:=`(lo = pred - 1.96 * sd_pred, hi = pred + 1.96 * sd_pred)]
  }

  max_val <- max(c(dt$obs, dt$pred), na.rm = TRUE) * 1.05

  p <- ggplot(dt, aes(obs, pred)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey50") +
    geom_point(aes(colour = get(facet_col)), alpha = 0.6, size = 2.5) +
    geom_smooth(method = "loess", se = TRUE, colour = "black",
                linewidth = 0.8, fill = "grey80", alpha = 0.3) +
    scale_colour_manual(values = col_pal, name = NULL) +
    labs(title = title,
         subtitle = "LOESS smoother (black) vs identity line (dashed)",
         x = "Observed counts", y = "Predicted counts") +
    THEME_PUB

  if (length(unique(dt[[facet_col]])) > 1) {
    p <- p + facet_wrap(as.formula(paste("~", facet_col)), scales = "free")
  }
  p
}

# ── Helper: build residuals-vs-fitted panel ──────────────────────────────────
make_resid_panel <- function(dt, obs_col, pred_col, title, col_pal,
                             facet_col = "wave_label") {
  dt <- copy(dt)
  setnames(dt, c(obs_col, pred_col), c("obs", "pred"), skip_absent = TRUE)
  dt[, resid := obs - pred]

  p <- ggplot(dt, aes(pred, resid)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
    geom_point(aes(colour = get(facet_col)), alpha = 0.6, size = 2.5) +
    geom_smooth(method = "loess", se = TRUE, colour = "black",
                linewidth = 0.8, fill = "grey80", alpha = 0.3) +
    scale_colour_manual(values = col_pal, name = NULL) +
    labs(title = title,
         subtitle = "LOESS smoother (black); flat at zero = no residual pattern",
         x = "Predicted counts", y = "Residual (Observed \u2212 Predicted)") +
    THEME_PUB

  if (length(unique(dt[[facet_col]])) > 1) {
    p <- p + facet_wrap(as.formula(paste("~", facet_col)), scales = "free")
  }
  p
}

# ══════════════════════════════════════════════════════════════════════════════
# NYC — MODZCTA-level temporal CV area totals (out-of-sample)
# ══════════════════════════════════════════════════════════════════════════════
cat("── NYC MODZCTA-level diagnostics ──\n")
ny_waves <- c("NY_wave1", "NY_wave2", "NY_wave3")
ny_labels <- c("Wave 1", "Wave 2", "Wave 3")

ny_area <- rbindlist(lapply(seq_along(ny_waves), function(i) {
  f <- file.path(BASE, ny_waves[i], "validation", "modzcta_temporal_cv_modzcta_totals.csv")
  d <- fread(f)
  d[, wave_label := ny_labels[i]]
  d
}))

p_ny1 <- make_pvo_panel(
  ny_area, "obs_total", "pred_total", "sd_total",
  title = "NYC — Predicted vs Observed (MODZCTA-level, out-of-sample)",
  col_pal = WAVE_COLS
)

p_ny2 <- make_resid_panel(
  ny_area, "obs_total", "pred_total",
  title = "NYC — Residuals vs Predicted (MODZCTA-level)",
  col_pal = WAVE_COLS
)

p_nyc <- p_ny1 / p_ny2 + plot_annotation(
  title = "Non-linearity Diagnostics — New York City",
  theme = theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 13))
)

ggsave(file.path(OUT, "nonlinearity_diag_nyc.png"), p_nyc,
       width = 12, height = 10, dpi = 300, bg = "white")
cat("  Saved nonlinearity_diag_nyc.png\n")

# ══════════════════════════════════════════════════════════════════════════════
# Utah — Daily temporal CV totals (out-of-sample, ~46 days per wave)
# ══════════════════════════════════════════════════════════════════════════════
cat("── Utah daily-level diagnostics ──\n")
ut_waves <- c("UT_wave1", "UT_wave2", "UT_wave3")
ut_labels <- c("Wave 1", "Wave 2", "Wave 3")

ut_daily <- rbindlist(lapply(seq_along(ut_waves), function(i) {
  # New CV schema emits horizon-3 daily forecasts; fall back to legacy h5 if present.
  cand <- c(
    file.path(BASE, ut_waves[i], "validation", "temporal_cv_daily_h3.csv"),
    file.path(BASE, ut_waves[i], "validation", "temporal_cv_daily_h5.csv")
  )
  f <- cand[file.exists(cand)][1]
  if (is.na(f)) {
    warning("[UT] no temporal_cv_daily_*.csv for ", ut_waves[i], " — skipping.")
    return(data.table())
  }
  d <- fread(f)
  d[, wave_label := ut_labels[i]]
  d
}), fill = TRUE)

p_ut1 <- make_pvo_panel(
  ut_daily, "obs", "pred", "sd",
  title = "Utah — Predicted vs Observed (Daily totals, out-of-sample)",
  col_pal = c("Wave 1" = "#B2182B", "Wave 2" = "#D6604D", "Wave 3" = "#F4A582")
)

p_ut2 <- make_resid_panel(
  ut_daily, "obs", "pred",
  title = "Utah — Residuals vs Predicted (Daily totals)",
  col_pal = c("Wave 1" = "#B2182B", "Wave 2" = "#D6604D", "Wave 3" = "#F4A582")
)

p_utah <- p_ut1 / p_ut2 + plot_annotation(
  title = "Non-linearity Diagnostics — Utah",
  theme = theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 13))
)

ggsave(file.path(OUT, "nonlinearity_diag_utah.png"), p_utah,
       width = 12, height = 10, dpi = 300, bg = "white")
cat("  Saved nonlinearity_diag_utah.png\n")

# ══════════════════════════════════════════════════════════════════════════════
# Combined compact figure (1 row per city, 2 cols: P-v-O + Residual)
# ══════════════════════════════════════════════════════════════════════════════
cat("── Combined figure ──\n")

# NYC — single panel (all waves pooled with colour)
p_cny1 <- make_pvo_panel(
  ny_area, "obs_total", "pred_total", "sd_total",
  title = "NYC — Predicted vs Observed",
  col_pal = WAVE_COLS
) + facet_null() + theme(legend.position = "right")

p_cny2 <- make_resid_panel(
  ny_area, "obs_total", "pred_total",
  title = "NYC — Residuals vs Predicted",
  col_pal = WAVE_COLS
) + facet_null() + theme(legend.position = "right")

p_cut1 <- make_pvo_panel(
  ut_daily, "obs", "pred", "sd",
  title = "Utah — Predicted vs Observed",
  col_pal = c("Wave 1" = "#B2182B", "Wave 2" = "#D6604D", "Wave 3" = "#F4A582")
) + facet_null() + theme(legend.position = "right")

p_cut2 <- make_resid_panel(
  ut_daily, "obs", "pred",
  title = "Utah — Residuals vs Predicted",
  col_pal = c("Wave 1" = "#B2182B", "Wave 2" = "#D6604D", "Wave 3" = "#F4A582")
) + facet_null() + theme(legend.position = "right")

p_combined <- (p_cny1 | p_cny2) / (p_cut1 | p_cut2) +
  plot_annotation(
    title = "Non-linearity Diagnostics — Predicted vs Observed",
    subtitle = "Out-of-sample area-level predictions with LOESS smoother",
    theme = theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 13),
      plot.subtitle = element_text(hjust = 0.5, color = "grey40", size = 10)
    )
  )

ggsave(file.path(OUT, "nonlinearity_diagnostics_combined.png"), p_combined,
       width = 14, height = 10, dpi = 300, bg = "white")
cat("  Saved nonlinearity_diagnostics_combined.png\n")

cat("── Done ──\n")
