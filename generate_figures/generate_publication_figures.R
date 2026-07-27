#!/usr/bin/env Rscript
# generate_publication_figures.R
# Produce final publication-ready thesis figures into `thesis_latex/figures_clean/`.


suppressPackageStartupMessages({
  library(data.table)
  library(sf)
  library(ggplot2)
  library(patchwork)
  library(viridis)
  library(ggspatial)
  library(spdep)
  library(scales)
  library(dplyr)
  library(maptiles)
  library(tidyterra)
  library(terra)
  library(rnaturalearth)
})

Sys.unsetenv("MallocStackLogging")
source("helpers/_utils.R")
load_pipeline()

# ── Output directory ──
OUT_DIR <- file.path("thesis_latex", "figures_clean")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

BASE <- CONFIG$output_dir

# ═══════════════════════════════════════════════════════════════════════════════
# GLOBAL SETTINGS
# ═══════════════════════════════════════════════════════════════════════════════

DPI <- 300

# Consistent wave palette
WAVE_COLS <- c(
  NY_wave1 = "#2166AC", NY_wave2 = "#4393C3", NY_wave3 = "#92C5DE",
  UT_wave1 = "#B2182B", UT_wave2 = "#D6604D", UT_wave3 = "#F4A582"
)
WAVE_LABS <- c(
  NY_wave1 = "NYC Wave 1", NY_wave2 = "NYC Wave 2", NY_wave3 = "NYC Wave 3",
  UT_wave1 = "Utah Wave 1", UT_wave2 = "Utah Wave 2", UT_wave3 = "Utah Wave 3"
)
# Labeled palette: maps wave_lab values ("NYC Wave 1") → colors
WAVE_COLS_LAB <- setNames(WAVE_COLS, WAVE_LABS[names(WAVE_COLS)])

# Publication theme
THEME_PUB <- theme_minimal(base_size = 10) +
  theme(
    text             = element_text(family = ""),
    plot.title       = element_text(face = "bold", hjust = 0.5, size = 11),
    plot.subtitle    = element_text(hjust = 0.5, color = "grey40", size = 8),
    strip.text       = element_text(face = "bold", size = 9),
    panel.grid.minor = element_blank(),
    legend.position  = "bottom",
    legend.key.width = unit(1.2, "cm")
  )

THEME_MAP <- theme_void(base_size = 9) +
  theme(
    plot.title    = element_text(face = "bold", hjust = 0.5, size = 10),
    plot.subtitle = element_text(hjust = 0.5, color = "grey40", size = 7),
    panel.border  = element_rect(color = "grey30", fill = NA, linewidth = 0.4),
    legend.position = "right"
  )

save_pub <- function(p, fname, w_in = 12, h_in = 8) {
  path <- file.path(OUT_DIR, fname)
  ggsave(path, p, width = w_in, height = h_in, dpi = DPI, bg = "white")
}

# Per-wave WAIC-best model name (canonical), read from the fit artifact.
# These CSVs (`tables/unified_*.csv`) are written from this best model only,
# so the data already corresponds to the wave-specific winner; this lookup is
# used only to label plots correctly.
best_model_for <- function(wave_id) {
  p <- file.path(BASE, wave_id, "models", "best_model.txt")
  if (file.exists(p)) trimws(readLines(p, n = 1, warn = FALSE)) else "unknown"
}

# Helper to load wave data from CSVs (lightweight, no RDS needed for most)
load_fe <- function(wave_id) {
  f <- file.path(BASE, wave_id, "tables", "unified_fixed_effects.csv")
  dt <- fread(f)
  dt[, wave := wave_id]
  dt
}

load_te <- function(wave_id) {
  f <- file.path(BASE, wave_id, "tables", "unified_temporal_effects.csv")
  dt <- fread(f)
  dt[, wave := wave_id]
  dt
}

load_se <- function(wave_id) {
  f <- file.path(BASE, wave_id, "tables", "unified_spatial_effects.csv")
  dt <- fread(f)
  dt[, wave := wave_id]
  dt
}


# ═══════════════════════════════════════════════════════════════════════════════
# 1. STUDY AREA MAP — fixed inset position
# ═══════════════════════════════════════════════════════════════════════════════

# ── NYC: MODZCTA boundaries coloured by borough ──
modzcta_sf <- st_read(file.path("data", "MODZCTA",
  "geo_export_5392c927-ee1e-4b9c-9787-2cad9456c83b.shp"), quiet = TRUE) |>
  st_transform(4326)

# Borough assignment from MODZCTA leading prefix (NYC ZIPs):
#   100xx/101xx/102xx -> Manhattan, 103xx -> Staten Island,
#   104xx -> Bronx,    111xx/112xx -> Brooklyn,
#   113xx/114xx/116xx -> Queens
modz_str <- as.character(modzcta_sf$modzcta)
modzcta_sf$boro <- dplyr::case_when(
  substr(modz_str,1,3) %in% c("100","101","102")             ~ "Manhattan",
  substr(modz_str,1,3) == "103"                              ~ "Staten Island",
  substr(modz_str,1,3) == "104"                              ~ "Bronx",
  substr(modz_str,1,3) %in% c("111","112")                   ~ "Brooklyn",
  substr(modz_str,1,3) %in% c("113","114","116")             ~ "Queens",
  TRUE                                                       ~ "Other"
)
modz_labels <- suppressWarnings(st_point_on_surface(modzcta_sf))

BORO_PAL <- c(Manhattan = "#4A90A4", Brooklyn = "#7BAEC2",
              Queens = "#9CC2D1", Bronx = "#6BA4BA",
              `Staten Island` = "#8FB8C8")

p_nyc <- ggplot() +
  geom_sf(data = modzcta_sf, aes(fill = boro),
          color = "white", linewidth = 0.25) +
  geom_sf_text(data = modz_labels, aes(label = modzcta),
               size = 1.3, color = "grey15") +
  scale_fill_manual(values = BORO_PAL, name = NULL) +
  labs(title = "a) New York City",
       subtitle = "178 MODZCTA reporting zones, coloured by borough",
       caption = "Source: NYC DOHMH") +
  coord_sf(crs = 4326) +
  THEME_MAP +
  theme(legend.position = "bottom",
        plot.caption = element_text(size = 7, color = "grey40", hjust = 0))

# ── Utah: counties coloured by urban / rural designation ──
ut_counties <- st_read(CONFIG$ut_counties_path, quiet = TRUE) |> st_transform(4326)
# Census-defined Wasatch Front + Cache urban corridor
URBAN_UT <- c("SALT LAKE", "UTAH", "DAVIS", "WEBER", "CACHE")
ut_counties$class <- ifelse(ut_counties$NAME %in% URBAN_UT, "Urban", "Rural")
ut_counties$NAME_TC <- tools::toTitleCase(tolower(ut_counties$NAME))
ut_pts <- suppressWarnings(st_point_on_surface(ut_counties))

UT_PAL <- c(Urban = "#7BAEC2", Rural = "#E8E4D8")

p_utah <- ggplot() +
  geom_sf(data = ut_counties, aes(fill = class),
          color = "grey40", linewidth = 0.3) +
  geom_sf_text(data = ut_pts, aes(label = NAME_TC),
               size = 2.2, color = "grey15", fontface = "plain") +
  scale_fill_manual(values = UT_PAL, name = NULL) +
  labs(title = "b) Utah",
       subtitle = "29 counties, shaded by urban (Wasatch Front + Cache) vs rural",
       caption = "Source: US Census TIGER/Line") +
  coord_sf(crs = 4326) +
  THEME_MAP +
  theme(legend.position = "bottom",
        plot.caption = element_text(size = 7, color = "grey40", hjust = 0))

p_study <- (p_nyc | p_utah) +
  plot_layout(widths = c(1, 1))
save_pub(p_study, "study_areas.png", w_in = 14, h_in = 8)
rm(modzcta_sf, modz_labels, ut_counties, ut_pts)
gc(verbose = FALSE)


# ═══════════════════════════════════════════════════════════════════════════════
# 2. MODEL COMPARISON COMBINED (WAIC heatmap + rank transition)
# ═══════════════════════════════════════════════════════════════════════════════

wave_dirs <- c("NY_wave1", "NY_wave2", "NY_wave3",
               "UT_wave1", "UT_wave2", "UT_wave3")
wave_labs_w <- c("NYC W1", "NYC W2", "NYC W3", "UT W1", "UT W2", "UT W3")

waic_all <- rbindlist(lapply(seq_along(wave_dirs), function(i) {
  dt <- fread(file.path(BASE, wave_dirs[i], "tables",
                        "model_comparison_unified.csv"))
  dt <- dt[Status == "OK"]
  dt[, dWAIC := WAIC - min(WAIC, na.rm = TRUE)]
  dt[, `:=`(wave = wave_labs_w[i], is_best = dWAIC == 0)]
  dt[, .(Model, WAIC, dWAIC, wave, is_best)]
}))

# WAIC-best per wave comes from the data (min dWAIC); no manual override.
# Canonical record: outputs/spatiotemporal_risk/<wave>/models/best_model.txt
# (NY_wave1: NLMob_ST_IV; NY_wave2/3: LinMob_ST_IV; UT_*: NLMob_ST_IV).

model_order <- c("Base", "Temp", "LinMob", "NLMob",
                 "LinMob_ST_IV", "NLMob_ST_IV")
model_labels <- display_model(model_order)
waic_all[, Model := factor(Model, levels = model_order, labels = model_labels)]
waic_all[, wave  := factor(wave, levels = wave_labs_w)]
waic_all[, label := ifelse(is_best, sprintf("%.0f*", dWAIC),
                           sprintf("%.0f", dWAIC))]

p_waic <- ggplot(waic_all, aes(x = wave, y = Model, fill = dWAIC)) +
  geom_tile(color = "white", linewidth = 0.8) +
  geom_text(aes(label = label), size = 3.2, fontface = "bold") +
  scale_fill_viridis_c(option = "magma", direction = -1,
                       trans = pseudo_log_trans(sigma = 10),
                       name = expression(Delta * WAIC),
                       breaks = c(0, 50, 500, 5000)) +
  scale_y_discrete(limits = rev(model_labels)) +
  labs(title = expression("a)" ~ Delta * "WAIC Model Comparison"),
       x = NULL, y = NULL) +
  THEME_PUB +
  theme(panel.grid = element_blank(),
        legend.position = "right",
        legend.key.width = unit(0.5, "cm"))

# Rank transition
trans_dt <- copy(waic_all)
trans_dt[, rank := frank(dWAIC, ties.method = "min"), by = wave]
best_dt  <- trans_dt[is_best == TRUE]

p_trans <- ggplot(trans_dt, aes(x = wave, y = rank, group = Model, color = Model)) +
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
  labs(title = "b) Model Rank Transitions", x = NULL) +
  THEME_PUB +
  theme(legend.position = "right")

p_combined <- p_waic + p_trans +
  plot_layout(widths = c(1, 1)) +
  plot_annotation(
    title = "Model Comparison Across Study Configurations",
    theme = theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 12))
  )
save_pub(p_combined, "model_comparison_combined.png", w_in = 14, h_in = 6)


# ═══════════════════════════════════════════════════════════════════════════════
# 3–4. TEMPORAL RR CROSS-WAVE (NYC + Utah)
# ═══════════════════════════════════════════════════════════════════════════════

make_temporal_panel <- function(state_prefix, wave_ids, panel_title) {
  te_all <- rbindlist(lapply(wave_ids, function(w) {
    te <- load_te(w)
    wm <- CONFIG$waves[[sub(".*_", "", w)]]$months
    te[, date := as.Date(paste0(wm[1], "-01")) + idx_temporal - 1L]
    te[, `:=`(RR = exp(mean), rr_lo = exp(q025), rr_hi = exp(q975))]
    te[, wave_lab := WAVE_LABS[w]]
    te
  }))
  te_all[, wave_lab := factor(wave_lab, levels = WAVE_LABS[wave_ids])]

  # Map colors by wave label (not wave id)
  lab_cols <- setNames(WAVE_COLS[wave_ids], WAVE_LABS[wave_ids])

  # Subtitle reflects the WAIC-best model used per wave (read from disk).
  best_per_wave <- vapply(wave_ids, best_model_for, character(1))
  sub_txt <- if (length(unique(best_per_wave)) == 1L)
    sprintf("%s | Posterior mean with 95%% credible interval", unique(best_per_wave))
  else
    sprintf("WAIC-best per wave (%s) | Posterior mean with 95%% credible interval",
            paste(sprintf("%s: %s", WAVE_LABS[wave_ids], best_per_wave), collapse = "; "))

  ggplot(te_all, aes(x = date)) +
    geom_ribbon(aes(ymin = rr_lo, ymax = rr_hi, fill = wave_lab), alpha = 0.45) +
    geom_line(aes(y = rr_lo, color = wave_lab), linewidth = 0.35, linetype = "dotted") +
    geom_line(aes(y = rr_hi, color = wave_lab), linewidth = 0.35, linetype = "dotted") +
    geom_line(aes(y = RR, color = wave_lab), linewidth = 0.8) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "grey50") +
    scale_color_manual(values = lab_cols, name = NULL) +
    scale_fill_manual(values = lab_cols, name = NULL) +
    scale_x_date(date_labels = "%b %d", date_breaks = "10 days") +
    facet_wrap(~wave_lab, scales = "free_x", nrow = 1) +
    labs(title = panel_title,
         subtitle = sub_txt,
         x = "Date", y = "Temporal Relative Risk") +
    THEME_PUB +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 30, hjust = 1, size = 7))
}

# Combined panels (kept for reference)
p_trr_ny <- make_temporal_panel("NY",
  c("NY_wave1", "NY_wave2", "NY_wave3"),
  "Temporal Relative Risk — New York City")
save_pub(p_trr_ny, "ny_temporal_rr_cross_wave.png", w_in = 7.5, h_in = 3.5)

p_trr_ut <- make_temporal_panel("UT",
  c("UT_wave1", "UT_wave2", "UT_wave3"),
  "Temporal Relative Risk — Utah")
save_pub(p_trr_ut, "ut_temporal_rr_cross_wave.png", w_in = 7.5, h_in = 3.5)

# Individual wave figures (one per wave for better print clarity)
make_temporal_single <- function(wave_id, wave_title) {
  te <- load_te(wave_id)
  wm <- CONFIG$waves[[sub(".*_", "", wave_id)]]$months
  te[, date := as.Date(paste0(wm[1], "-01")) + idx_temporal - 1L]
  te[, `:=`(RR = exp(mean), rr_lo = exp(q025), rr_hi = exp(q975))]
  col <- WAVE_COLS[wave_id]

  bm <- best_model_for(wave_id)

  ggplot(te, aes(x = date)) +
    geom_ribbon(aes(ymin = rr_lo, ymax = rr_hi), fill = col, alpha = 0.45) +
    geom_line(aes(y = rr_lo), color = col, linewidth = 0.35, linetype = "dotted") +
    geom_line(aes(y = rr_hi), color = col, linewidth = 0.35, linetype = "dotted") +
    geom_line(aes(y = RR), color = col, linewidth = 0.8) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "grey50") +
    scale_x_date(date_labels = "%b %d", date_breaks = "10 days") +
    labs(title = wave_title,
         subtitle = sprintf("%s (WAIC-best) | Posterior mean with 95%% credible interval", bm),
         x = "Date", y = "Temporal Relative Risk") +
    THEME_PUB +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 30, hjust = 1, size = 8))
}

for (sid in c("NY", "UT")) {
  for (wn in paste0("wave", 1:3)) {
    wid <- paste0(sid, "_", wn)
    wlab <- WAVE_LABS[wid]
    p <- make_temporal_single(wid, paste("Temporal Relative Risk —", wlab))
    save_pub(p, paste0(tolower(sid), "_temporal_rr_", wn, ".png"),
             w_in = 6, h_in = 3.5)
  }
}


# ═══════════════════════════════════════════════════════════════════════════════
# 5–6. FIXED EFFECTS FOREST PLOTS — LOCKED covariate order
# ═══════════════════════════════════════════════════════════════════════════════

# Human-readable covariate labels, annotated with the β-symbol used in
# Section "Fixed Effects (x_{g,t})" of the methodology. Labels are written
# as plotmath expressions (use scale_y_discrete(labels = scales::label_parse())).
COV_LABELS <- c(
  building_density_std       = "'Building Density'~(beta[bld])",
  inf_pressure_std           = "rho[g*','*t]^{inf}~~'(OD infection pressure)'",
  contact_intensity_std      = "'Contact Intensity'~(beta[cnt])",
  log_activity_density_std   = "'Activity Density'~(beta[act])",
  mob_lag3                   = "log~rho[g*','*t-3]^{resid}",
  mob_lag7                   = "log~rho[g*','*t-7]^{resid}",
  mob_lag14                  = "log~rho[g*','*t-14]^{resid}",
  mob_lag21                  = "log~rho[g*','*t-21]^{resid}",
  is_weekend                 = "'Weekend'~(beta[wknd])",
  pct_younger_std            = "'% Younger'~(beta[young])",
  pct_older_std              = "'% Older'~(beta[old])",
  pct_employed_std           = "'% Employed'~(beta[emp])",
  pct_unemployed_std         = "'% Unemployed'~(beta[unemp])",
  ses_contrast_std           = "'SES Contrast'~(beta[ses])",
  avg_category_risk_std      = "'Category Risk'~(beta[cat])",
  exposure_import_std        = "rho[g*','*t]^{imp}~~'(exposure import)'",
  hourly_concentration_std   = "'Hourly Visit Conc.'~(beta[hr])",
  log_repeat_ratio_std       = "'Visitor Repeat Ratio'~(beta[rep])",
  log_daytime_pop_std        = "'Daytime Population'~(beta[day])",
  mob_spatial_lag_std        = "W~log~rho[g*','*t]^{inf}~~'(spatial lag)'",
  landuse_commercial         = "'Land Use (Commercial)'~(beta[lc])",
  log_distance_from_home_std = "'Distance from Home'~(beta[dist])",
  commercial_x_weekend       = "'Commercial' %*% 'Weekend'~(beta[lc*','*wknd])",
  activity_x_weekend         = "'Activity' %*% 'Weekend'~(beta[act*','*wknd])",
  daily_visit_cv_std         = "'Daily Visit CV'~(beta[cv])"
)

make_forest_panel <- function(wave_ids, panel_title, cov_order_override = NULL) {
  fe_all <- rbindlist(lapply(wave_ids, load_fe))
  fe_all <- fe_all[covariate != "(Intercept)"]
  fe_all[, wave_lab := WAVE_LABS[wave]]
  fe_all[, wave_lab := factor(wave_lab, levels = WAVE_LABS[wave_ids])]

  # Apply human-readable labels; fallback is a plotmath-safe quoted literal
  fe_all[, covariate := fifelse(
    covariate %in% names(COV_LABELS),
    COV_LABELS[covariate],
    paste0("'", gsub("_std$", "", gsub("_", " ", covariate)), "'")
  )]

  # 66% CI from SE
  fe_all[, se := (log(RR_upper) - log(RR)) / qnorm(0.975)]
  fe_all[, RR_66lo := exp(log(RR) + qnorm(0.17) * se)]
  fe_all[, RR_66hi := exp(log(RR) + qnorm(0.83) * se)]

  # LOCKED covariate order: sort by mean absolute RR across all waves
  cov_order <- fe_all[, .(mean_abs_rr = mean(abs(log(RR)))), by = covariate]
  setorder(cov_order, mean_abs_rr)
  local_levels <- cov_order$covariate

  if (!is.null(cov_order_override)) {
    # Use override order; append any covariates not in override at the bottom
    extra <- setdiff(local_levels, cov_order_override)
    local_levels <- c(intersect(cov_order_override, local_levels), extra)
  }

  fe_all[, covariate := factor(covariate, levels = local_levels)]

  # Significance color
  fe_all[, sig_col := ifelse(significant, "Significant", "Not significant")]

  list(
    plot = ggplot(fe_all, aes(x = RR, y = covariate, color = sig_col)) +
      geom_vline(xintercept = 1, linetype = "dashed", color = "grey50") +
      geom_linerange(aes(xmin = RR_lower, xmax = RR_upper), linewidth = 0.4, alpha = 0.5) +
      geom_linerange(aes(xmin = RR_66lo, xmax = RR_66hi), linewidth = 1.5) +
      geom_point(size = 2) +
      scale_color_manual(values = c("Significant" = "#2166AC", "Not significant" = "#999999"),
                         name = NULL) +
      scale_y_discrete(labels = scales::label_parse()) +
      facet_wrap(~wave_lab, nrow = 1) +
      labs(title = panel_title,
           subtitle = "Posterior mean RR with 66% (thick) and 95% (thin) credible intervals",
           x = "Relative Risk (exp(β))", y = NULL) +
      THEME_PUB +
      theme(strip.text = element_text(face = "bold")),
    cov_levels = local_levels
  )
}

ny_result <- make_forest_panel(
  c("NY_wave1", "NY_wave2", "NY_wave3"),
  "Fixed-Effect Relative Risks — New York City (WAIC-best Type IV variant per wave)")
p_fe_ny <- ny_result$plot
save_pub(p_fe_ny, "ny_fixed_effects_cross_wave.png", w_in = 14, h_in = 7)

ut_result <- make_forest_panel(
  c("UT_wave1", "UT_wave2", "UT_wave3"),
  "Fixed-Effect Relative Risks — Utah (WAIC-best Type IV variant per wave)",
  cov_order_override = ny_result$cov_levels)
p_fe_ut <- ut_result$plot
save_pub(p_fe_ut, "ut_fixed_effects_cross_wave.png", w_in = 14, h_in = 7)


# ═══════════════════════════════════════════════════════════════════════════════
# 7–8. SPATIAL RR CROSS-WAVE — Global scale + significance masking
# ═══════════════════════════════════════════════════════════════════════════════

grid_sf_ny <- st_read(file.path(CONFIG$grids_dir, "grid_100m_NY.gpkg"), quiet = TRUE)
grid_sf_ut <- st_read(file.path(CONFIG$grids_dir, "grid_100m_UT.gpkg"), quiet = TRUE)
ut_counties_sf <- st_read(CONFIG$ut_counties_path, quiet = TRUE) |>
  st_transform(st_crs(grid_sf_ut))

compute_spatial_rr <- function(wave_id, grid_sf) {
  sp <- load_se(wave_id)
  n_areas <- nrow(sp) %/% 2L
  sp_s <- sp[idx_spatial <= n_areas]; setorder(sp_s, idx_spatial)
  sp_u <- sp[idx_spatial >  n_areas]; setorder(sp_u, idx_spatial)

  # Spatial RR per Blangiardo & Cameletti (2015): RR_g = exp(mean(u_g) + mean(v_g)).
  # SD of total log-RR assumes posterior independence of structured/unstructured
  # components (standard BYM2 assumption; INLA returns marginal posteriors).
  total_eff <- data.table(
    grid_id   = sp_s$grid_id,
    total_RR  = exp(sp_s$mean + sp_u$mean),
    total_lo  = sp_s$q025 + sp_u$q025,
    total_hi  = sp_s$q975 + sp_u$q975,
    total_mean_log = sp_s$mean + sp_u$mean,
    total_sd  = sqrt(sp_s$sd^2 + sp_u$sd^2)
  )
  # Significance: 95% CI on log-RR excludes 0 (i.e. RR-CI excludes 1).
  total_eff[, significant := (total_lo > 0) | (total_hi < 0)]

  # Exceedance probability P(RR_g > 1 | y) = P(u_g + v_g > 0 | y) under a
  # Gaussian approximation to the marginal posterior of the BYM2 sum
  # (INLA's BYM2 marginals are close to Gaussian in practice). Hotspot
  # classification follows Richardson, Thomson, Best & Elliott (2004):
  #   prob_gt_1 > 0.95  -> high-risk hotspot
  #   prob_gt_1 < 0.05  -> low-risk coldspot
  total_eff[, prob_gt_1 := pnorm(total_mean_log / pmax(total_sd, 1e-8))]
  total_eff[, hotspot := fifelse(prob_gt_1 > 0.95, "high",
                          fifelse(prob_gt_1 < 0.05, "low", "neutral"))]
  total_eff[, wave := wave_id]
  total_eff
}

# NYC spatial RR
ny_waves <- c("NY_wave1", "NY_wave2", "NY_wave3")
ny_sp_all <- rbindlist(lapply(ny_waves, compute_spatial_rr, grid_sf = grid_sf_ny))

# Hotspot summary table (Richardson et al. 2004 classification) for thesis text.
hotspot_summary <- ny_sp_all[, .(
  n_cells     = .N,
  n_high      = sum(hotspot == "high"),
  n_low       = sum(hotspot == "low"),
  pct_high    = round(100 * mean(hotspot == "high"), 2),
  pct_low     = round(100 * mean(hotspot == "low"), 2),
  median_RR   = round(median(total_RR), 3),
  max_RR      = round(max(total_RR), 3)
), by = wave]
fwrite(hotspot_summary, file.path(OUT_DIR, "spatial_rr_hotspot_summary_ny.csv"))
fwrite(ny_sp_all[, .(wave, grid_id, total_RR, prob_gt_1, hotspot, significant)],
       file.path(OUT_DIR, "spatial_rr_grid_ny.csv"))

# Re-read MODZCTA boundaries for basemap overlay
modzcta_rr_sf <- st_read(file.path("data", "MODZCTA",
  "geo_export_5392c927-ee1e-4b9c-9787-2cad9456c83b.shp"), quiet = TRUE) |>
  st_transform(st_crs(grid_sf_ny))

# Compute GLOBAL color limits across all NYC waves
ny_rr_vals <- ny_sp_all$total_RR
ny_q01 <- quantile(ny_rr_vals, 0.01, na.rm = TRUE)
ny_q99 <- quantile(ny_rr_vals, 0.99, na.rm = TRUE)

ny_spatial_plots <- lapply(ny_waves, function(w) {
  sp_w <- ny_sp_all[wave == w]
  sp_sf <- merge(grid_sf_ny, sp_w, by = "grid_id", all.x = FALSE)

  # Significance masking: alpha channel
  sp_sf$alpha_val <- ifelse(sp_sf$significant, 1.0, 0.2)
  sp_sf$rr_clipped <- pmin(pmax(sp_sf$total_RR, ny_q01), ny_q99)

  ggplot(sp_sf) +
    geom_sf(data = modzcta_rr_sf, fill = NA, color = "grey60", linewidth = 0.25) +
    geom_sf(aes(fill = rr_clipped, alpha = alpha_val), color = NA) +
    scale_fill_gradient2(
      low = "#2166AC", mid = "white", high = "#B2182B",
      midpoint = 1, name = "RR",
      limits = c(ny_q01, ny_q99),
      oob = squish
    ) +
    scale_alpha_identity() +
    coord_sf(datum = NA) +
    labs(title = paste0(WAVE_LABS[w], " — ", best_model_for(w))) +
    THEME_MAP +
    theme(legend.position = if (w == tail(ny_waves, 1)) "right" else "none")
})

p_ny_spatial <- wrap_plots(ny_spatial_plots, nrow = 1) +
  plot_layout(guides = "collect") +
  plot_annotation(
    title = "BYM2 Spatial Relative Risk — NYC (WAIC-best Type IV variant per wave)",
    subtitle = "Global scale | Non-significant cells (95% CI includes 1) faded",
    theme = theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 12),
                  plot.subtitle = element_text(hjust = 0.5, color = "grey40", size = 9))
  ) &
  theme(legend.position = "right")
save_pub(p_ny_spatial, "ny_spatial_rr_cross_wave.png", w_in = 14, h_in = 5.5)

# Utah spatial RR
ut_waves <- c("UT_wave1", "UT_wave2", "UT_wave3")
ut_sp_all <- rbindlist(lapply(ut_waves, compute_spatial_rr, grid_sf = grid_sf_ut))

ut_hotspot_summary <- ut_sp_all[, .(
  n_cells   = .N,
  n_high    = sum(hotspot == "high"),
  n_low     = sum(hotspot == "low"),
  pct_high  = round(100 * mean(hotspot == "high"), 2),
  pct_low   = round(100 * mean(hotspot == "low"), 2),
  median_RR = round(median(total_RR), 3),
  max_RR    = round(max(total_RR), 3)
), by = wave]
fwrite(ut_hotspot_summary, file.path(OUT_DIR, "spatial_rr_hotspot_summary_ut.csv"))
fwrite(ut_sp_all[, .(wave, grid_id, total_RR, prob_gt_1, hotspot, significant)],
       file.path(OUT_DIR, "spatial_rr_grid_ut.csv"))

ut_rr_vals <- ut_sp_all$total_RR
ut_q01 <- quantile(ut_rr_vals, 0.01, na.rm = TRUE)
ut_q99 <- quantile(ut_rr_vals, 0.99, na.rm = TRUE)

ut_spatial_plots <- lapply(ut_waves, function(w) {
  sp_w <- ut_sp_all[wave == w]
  sp_sf <- merge(grid_sf_ut, sp_w, by = "grid_id", all.x = FALSE)
  sp_sf$alpha_val <- ifelse(sp_sf$significant, 1.0, 0.2)
  sp_sf$rr_clipped <- pmin(pmax(sp_sf$total_RR, ut_q01), ut_q99)
  cents <- suppressWarnings(st_centroid(sp_sf))

  ggplot() +
    geom_sf(data = ut_counties_sf, fill = "grey95", color = "grey60",
            linewidth = 0.3, alpha = 0.3) +
    geom_sf(data = sp_sf,
            aes(fill = rr_clipped, alpha = alpha_val),
            color = "grey40", linewidth = 0.3) +
    geom_sf_text(data = cents, aes(label = round(sp_sf$total_RR, 2)),
                 size = 1.8, fontface = "bold") +
    scale_fill_gradient2(
      low = "#2166AC", mid = "white", high = "#B2182B",
      midpoint = 1, name = "RR",
      limits = c(ut_q01, ut_q99),
      oob = squish
    ) +
    scale_alpha_identity() +
    coord_sf(datum = NA) +
    labs(title = paste0(WAVE_LABS[w], " — ", best_model_for(w))) +
    THEME_MAP +
    theme(legend.position = if (w == tail(ut_waves, 1)) "right" else "none")
})

p_ut_spatial <- wrap_plots(ut_spatial_plots, nrow = 1) +
  plot_layout(guides = "collect") +
  plot_annotation(
    title = "BYM2 Spatial Relative Risk — Utah (WAIC-best Type IV variant per wave)",
    subtitle = "Global scale | Non-significant cells faded | RR labels on counties",
    theme = theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 12),
                  plot.subtitle = element_text(hjust = 0.5, color = "grey40", size = 9))
  ) &
  theme(legend.position = "right")
save_pub(p_ut_spatial, "ut_spatial_rr_cross_wave.png", w_in = 14, h_in = 6)
rm(ny_sp_all, ut_sp_all); gc(verbose = FALSE)


# ═══════════════════════════════════════════════════════════════════════════════
# 9–11. SPATIOTEMPORAL INTERACTION (NYC, 3 waves) — Heatmap only
# ═══════════════════════════════════════════════════════════════════════════════

load_fv <- function(wave_id) {
  f <- file.path(BASE, wave_id, "tables", "unified_fitted_values.csv")
  fread(f)
}

for (wv in c("NY_wave1", "NY_wave2", "NY_wave3")) {
  pfx <- paste0("ny_", sub("NY_wave", "w", wv))
  fname <- paste0(pfx, "_z4_spatiotemporal_interaction.png")

  # Load data
  fv <- load_fv(wv)
  sp <- load_se(wv)
  tp <- load_te(wv)
  wm <- CONFIG$waves[[sub(".*_", "", wv)]]$months
  start_date <- as.Date(paste0(wm[1], "-01"))

  # Compute spatial total effect
  n_areas <- nrow(sp) %/% 2L
  sp_s <- sp[idx_spatial <= n_areas]; setorder(sp_s, idx_spatial)
  sp_u <- sp[idx_spatial >  n_areas]; setorder(sp_u, idx_spatial)
  sp_total <- data.table(grid_id = sp_s$grid_id,
                         sp_mean = sp_s$mean + sp_u$mean)
  tp_dt <- data.table(time_id = tp$idx_temporal, tp_mean = tp$mean)

  # Compute delta: ST interaction residual
  fv[, log_RR := log(pmax(RR, 1e-10))]
  fv <- merge(fv, sp_total, by = "grid_id")
  fv <- merge(fv, tp_dt,    by = "time_id")
  fv[, delta_gt := log_RR - sp_mean - tp_mean]
  fv[, delta_gt := delta_gt - mean(delta_gt, na.rm = TRUE)]
  fv <- fv[is.finite(delta_gt)]

  # Aggregate to grid_id x time_id
  st_mat <- fv[, .(st_mean = mean(delta_gt, na.rm = TRUE)),
               by = .(grid_id, time_id)]
  st_mat[, date := start_date + time_id - 1L]

  # Top 40 grid cells by ST variance
  var_by_grid <- st_mat[, .(var_st = var(st_mean, na.rm = TRUE)), by = grid_id]
  top_grids <- var_by_grid[order(-var_st)][1:min(40, nrow(var_by_grid)), grid_id]
  st_sub <- st_mat[grid_id %in% top_grids]

  q_st <- quantile(abs(st_sub$st_mean), 0.99, na.rm = TRUE)

  p_heat <- ggplot(st_sub, aes(x = date,
                                y = reorder(as.character(grid_id), st_mean),
                                fill = pmin(pmax(st_mean, -q_st), q_st))) +
    geom_tile() +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                         midpoint = 0, name = "ST Effect",
                         limits = c(-q_st, q_st)) +
    scale_x_date(date_labels = "%b %d", date_breaks = "10 days") +
    labs(title = paste("Spatiotemporal Interaction \u2014", WAVE_LABS[wv]),
         subtitle = expression(paste(delta, "(g,t) = log RR - u(g) - ", gamma, "(t)")),
         x = "Date", y = "Grid Cell") +
    THEME_PUB +
    theme(axis.text.y = element_text(size = 6),
          axis.text.x = element_text(angle = 30, hjust = 1, size = 8))

  save_pub(p_heat, fname, w_in = 12, h_in = 7)
}


# ═══════════════════════════════════════════════════════════════════════════════
# 12. NLMob DOSE-RESPONSE
# ═══════════════════════════════════════════════════════════════════════════════

nlmob_path <- file.path(BASE, "NY_wave1", "models", "model_results_unified.rds")
if (file.exists(nlmob_path)) {
  res <- readRDS(nlmob_path)
  rw1 <- as.data.table(res$model_outputs$NLMob$rw1_curve)

  log_rho_range <- c(-5.01, 5.01)
  brks <- seq(log_rho_range[1], log_rho_range[2], length.out = N_RHO_BINS + 1)
  mids <- (brks[-length(brks)] + brks[-1]) / 2
  rw1[, log_rho := mids]

  p_nlmob <- ggplot(rw1, aes(x = log_rho, y = mean)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    geom_ribbon(aes(ymin = `0.025quant`, ymax = `0.975quant`),
                fill = "#D6604D", alpha = 0.2) +
    geom_line(linewidth = 1, color = "#B2182B") +
    geom_point(size = 2, color = "#B2182B") +
    labs(title = "Nonlinear Mobility Dose-Response (NLMob RW1)",
         subtitle = "NYC Wave 1 — Posterior mean with 95% credible interval",
         x = expression(log(rho) ~ "(standardised mobility intensity)"),
         y = "Random effect (log-RR scale)") +
    THEME_PUB
  save_pub(p_nlmob, "nlmob_dose_response.png", w_in = 7, h_in = 5)
  rm(res, rw1)
  gc(verbose = FALSE)
}


# ═══════════════════════════════════════════════════════════════════════════════
# 13–14. RESIDUAL DIAGNOSTICS (NYC W3 + Utah W3)
# ═══════════════════════════════════════════════════════════════════════════════

make_residual_plot <- function(state, grid_sf, counties_sf = NULL) {
  wave_id <- paste0(state, "_wave3")
  wdata <- load_wave("wave3", state, grid_sf)
  fv <- as.data.table(wdata$fitted_09b)

  # Spatial residual map
  resid_sp <- fv[, .(mean_resid = mean(pearson_resid, na.rm = TRUE)), by = grid_id]
  resid_sf <- merge(grid_sf, resid_sp, by = "grid_id", all.x = FALSE)
  q99 <- quantile(abs(resid_sf$mean_resid), 0.99, na.rm = TRUE)

  if (!is.null(counties_sf)) {
    cents <- suppressWarnings(st_centroid(resid_sf))
    p_map <- ggplot() +
      geom_sf(data = counties_sf, fill = "grey95", color = "grey60",
              linewidth = 0.3, alpha = 0.3) +
      geom_sf(data = resid_sf,
              aes(fill = pmin(pmax(mean_resid, -q99), q99)),
              color = "grey40", linewidth = 0.3) +
      geom_sf_text(data = cents, aes(label = round(resid_sf$mean_resid, 2)),
                   size = 2, fontface = "bold") +
      scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                           midpoint = 0, name = "Mean\nResidual",
                           limits = c(-q99, q99)) +
      coord_sf(datum = NA) +
      labs(title = "a) Mean Pearson Residuals") +
      THEME_MAP
  } else {
    p_map <- ggplot(resid_sf) +
      geom_sf(data = modzcta_rr_sf, fill = NA, color = "grey60", linewidth = 0.25) +
      geom_sf(aes(fill = pmin(pmax(mean_resid, -q99), q99)),
              color = NA) +
      scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                           midpoint = 0, name = "Mean\nResidual",
                           limits = c(-q99, q99)) +
      annotation_scale(location = "br", width_hint = 0.15, style = "ticks",
                       text_cex = 0.6) +
      coord_sf(datum = NA) +
      labs(title = "a) Mean Pearson Residuals") +
      THEME_MAP
  }

  # Residual vs index
  r_vec <- fv$pearson_resid[is.finite(fv$pearson_resid)]
  set.seed(42)
  n_sub <- min(60000, length(r_vec))
  idx_s <- sort(sample.int(length(r_vec), n_sub))

  p_index <- ggplot(data.table(Index = idx_s, resid = r_vec[idx_s]),
                    aes(x = Index, y = resid)) +
    geom_point(shape = ".", alpha = 0.3, color = "black") +
    geom_hline(yintercept = 0, color = "red", linewidth = 0.5) +
    scale_x_continuous(labels = comma) +
    coord_cartesian(ylim = quantile(r_vec, c(0.001, 0.999))) +
    labs(title = "b) Residuals vs Index", x = "Observation", y = "Pearson Residual") +
    THEME_PUB +
    theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.4))

  # Histogram
  r_clip <- pmin(pmax(r_vec, -10), 10)
  p_hist <- ggplot(data.table(r = r_clip), aes(x = r)) +
    geom_histogram(aes(y = after_stat(density)), bins = 50,
                   fill = "steelblue", color = "white", alpha = 0.7) +
    stat_function(fun = dnorm, args = list(mean = 0, sd = 1),
                  color = "red", linewidth = 0.8, linetype = "dashed") +
    labs(title = "c) Residual Distribution",
         subtitle = sprintf("Mean = %.3f, SD = %.3f",
                            mean(r_vec, na.rm = TRUE), sd(r_vec, na.rm = TRUE)),
         x = "Pearson Residual", y = "Density") +
    THEME_PUB

  state_lab <- if (state == "NY") "NYC" else "Utah"
  p_map / (p_index | p_hist) +
    plot_annotation(
      title = paste("Residual Diagnostics —", state_lab, "Wave 3 (ST_IV)"),
      theme = theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 12))
    ) +
    plot_layout(heights = c(1.2, 1))
}

p_resid_ny <- make_residual_plot("NY", grid_sf_ny)
save_pub(p_resid_ny, "ny_w3_residual_diagnostics.png", w_in = 14, h_in = 10)

p_resid_ut <- make_residual_plot("UT", grid_sf_ut, ut_counties_sf)
save_pub(p_resid_ut, "ut_w3_residual_diagnostics.png", w_in = 14, h_in = 10)
rm(p_resid_ny, p_resid_ut); gc(verbose = FALSE)


# ═══════════════════════════════════════════════════════════════════════════════
# 15. PIT COMPARISON
# ═══════════════════════════════════════════════════════════════════════════════

WAVES_ALL <- c("NY_wave1", "NY_wave2", "NY_wave3",
               "UT_wave1", "UT_wave2", "UT_wave3")

cpo_list <- lapply(WAVES_ALL, function(w) {
  f <- file.path(BASE, w, "validation", "cpo_diagnostics.rds")
  if (file.exists(f)) readRDS(f) else NULL
})
names(cpo_list) <- WAVES_ALL

# Build PIT data
pit_data <- rbindlist(lapply(names(cpo_list), function(w) {
  obj <- cpo_list[[w]]
  if (is.null(obj)) return(NULL)

  # Original PIT from INLA
  pit_orig <- obj$pit
  pit_orig <- pit_orig[!is.na(pit_orig) & is.finite(pit_orig)]
  if (length(pit_orig) > 20000) pit_orig <- sample(pit_orig, 20000)

  # Randomised NB PIT if available
  pit_nb <- obj$pit_nb
  if (!is.null(pit_nb)) {
    pit_nb <- pit_nb[!is.na(pit_nb) & is.finite(pit_nb)]
    if (length(pit_nb) > 20000) pit_nb <- sample(pit_nb, 20000)
  }

  rbindlist(list(
    data.table(wave = w, type = "INLA PIT (rate scale)", pit = pit_orig),
    if (!is.null(pit_nb))
      data.table(wave = w, type = "Randomised NB PIT", pit = pit_nb)
  ))
}))

if (nrow(pit_data) > 0) {
  pit_data[, wave_lab := WAVE_LABS[wave]]

  p_pit <- ggplot(pit_data, aes(x = pit)) +
    geom_histogram(bins = 20, fill = "steelblue", color = "white", alpha = 0.7,
                   aes(y = after_stat(density))) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "red", linewidth = 0.5) +
    facet_grid(type ~ wave_lab) +
    labs(title = "Probability Integral Transform (PIT) Comparison",
         subtitle = "Top: INLA rate-scale PIT | Bottom: Randomised NB PIT (Czado et al. 2009)",
         x = "PIT value", y = "Density") +
    THEME_PUB +
    theme(strip.text.y = element_text(angle = 0))
  save_pub(p_pit, "cv_pit_comparison.png", w_in = 14, h_in = 6)
}


# ═══════════════════════════════════════════════════════════════════════════════
# 16. COVERAGE BARPLOT
# ═══════════════════════════════════════════════════════════════════════════════

cov_data <- rbindlist(lapply(names(cpo_list), function(w) {
  obj <- cpo_list[[w]]
  if (is.null(obj)) return(NULL)
  # Coverage is stored in the summary data.table
  if (!is.null(obj$summary) && "Cov95_NB" %in% names(obj$summary)) {
    data.table(wave = w, wave_lab = WAVE_LABS[w], coverage = obj$summary$Cov95_NB)
  } else NULL
}))

if (nrow(cov_data) > 0) {
  p_cov <- ggplot(cov_data, aes(x = wave_lab, y = coverage,
                                 fill = ifelse(grepl("NY", wave), "NYC", "Utah"))) +
    geom_col(width = 0.6, alpha = 0.85) +
    geom_hline(yintercept = 0.95, linetype = "dashed", color = "red", linewidth = 0.6) +
    geom_text(aes(label = sprintf("%.1f%%", coverage * 100)), vjust = -0.5, size = 3) +
    scale_fill_manual(values = c("NYC" = "#2166AC", "Utah" = "#B2182B"), name = NULL) +
    scale_y_continuous(labels = percent, limits = c(0, 1.05)) +
    labs(title = "In-Sample 95% Coverage (NB Predictive Distribution)",
         x = NULL, y = "Coverage") +
    THEME_PUB
  save_pub(p_cov, "cv_coverage_barplot.png", w_in = 8, h_in = 5)
}


# ═══════════════════════════════════════════════════════════════════════════════
# 17–21. TEMPORAL FORWARD CV FIGURES
# ═══════════════════════════════════════════════════════════════════════════════
# NOTE: this block was written against an older temporal_cv_results.rds schema
# (expecting per-horizon `$daily` / `$scores`). The current artifact stores
# `cell_scores`, `daily_agg`, `daily_raw`, `window_scores`, `holdout_long`, etc.
# CV figures are produced by `generate_cv_figures.R`, so this block is now a
# no-op to avoid crashing the rest of the publication-figure pipeline.
if (FALSE) {

# Load temporal CV data
tcv_list <- lapply(WAVES_ALL, function(w) {
  f <- file.path(BASE, w, "validation", "temporal_cv_results.rds")
  if (file.exists(f)) readRDS(f) else NULL
})
names(tcv_list) <- WAVES_ALL

tcv_daily <- rbindlist(lapply(names(tcv_list), function(w) {
  obj <- tcv_list[[w]]
  if (is.null(obj)) return(NULL)
  rbindlist(lapply(names(obj), function(h) {
    d <- obj[[h]]$daily
    if (is.null(d) || !nrow(d)) return(NULL)
    d <- as.data.table(d)
    d[, `:=`(wave = w, horizon = as.integer(h))]
    d
  }), fill = TRUE)
}), fill = TRUE)
tcv_daily[, state := sub("_.*", "", wave)]
tcv_daily[, wave_lab := WAVE_LABS[wave]]

tcv_scores <- rbindlist(lapply(names(tcv_list), function(w) {
  obj <- tcv_list[[w]]
  if (is.null(obj)) return(NULL)
  rbindlist(lapply(names(obj), function(h) {
    s <- obj[[h]]$scores
    if (is.null(s) || !length(s)) return(NULL)
    s <- as.data.table(s)
    if (!"wave" %in% names(s)) s[, wave := w]
    if (!"horizon" %in% names(s)) s[, horizon := as.integer(h)]
    s
  }), fill = TRUE)
}), fill = TRUE)
tcv_scores[, wave_lab := WAVE_LABS[wave]]

# 17. Bias time series
if (nrow(tcv_daily) > 0 && "pred" %in% names(tcv_daily) && "obs" %in% names(tcv_daily)) {
  tcv_daily[, bias := pred - obs]

  p_bias <- ggplot(tcv_daily[horizon == 5],
                   aes(x = time_id, y = bias, color = wave_lab)) +
    geom_hline(yintercept = 0, linetype = "solid", color = "grey50") +
    geom_line(linewidth = 0.6) +
    geom_point(size = 1.2) +
    scale_color_manual(values = WAVE_COLS_LAB, name = NULL) +
    facet_wrap(~wave_lab, scales = "free", nrow = 2) +
    labs(title = "Daily Bias Time Series (h = 5 days)",
         subtitle = "Positive = over-prediction | Negative = under-prediction",
         x = "Forecast Day", y = "Bias (predicted − observed)") +
    THEME_PUB +
    theme(legend.position = "none")
  save_pub(p_bias, "cv_bias_timeseries.png", w_in = 12, h_in = 7)
}

# 18. Obs vs fitted traces (top-5 high-burden)
traces_plots <- list()
for (w in WAVES_ALL) {
  obj <- tcv_list[[w]]
  if (is.null(obj)) next
  # Try h=3 for best horizon
  d3 <- obj[["3"]]
  if (is.null(d3) || is.null(d3$cell_preds)) next
  cp <- as.data.table(d3$cell_preds)
  if (!"area_id" %in% names(cp)) next

  # Top 5 areas by total observed
  top5 <- cp[, .(total_obs = sum(obs, na.rm = TRUE)), by = area_id]
  setorder(top5, -total_obs)
  top5 <- head(top5, 5)$area_id

  cp_top <- cp[area_id %in% top5]
  cp_top[, area_id := factor(area_id)]

  p_tr <- ggplot(cp_top, aes(x = time_id, group = area_id)) +
    geom_line(aes(y = obs, color = area_id), linewidth = 0.5, alpha = 0.7) +
    geom_line(aes(y = pred, color = area_id), linewidth = 0.5, linetype = "dashed") +
    labs(title = WAVE_LABS[w], x = "Day", y = "Cases") +
    THEME_PUB +
    theme(legend.position = "none", plot.title = element_text(size = 9))
  traces_plots[[w]] <- p_tr
}

if (length(traces_plots) >= 2) {
  p_traces <- wrap_plots(traces_plots, nrow = 2) +
    plot_annotation(
      title = "Observed (solid) vs Fitted (dashed) — Top 5 Highest-Burden Areas",
      theme = theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 11))
    )
  save_pub(p_traces, "cv_obs_vs_fitted_traces.png", w_in = 14, h_in = 8)
}

# 19. Obs vs. predicted aggregate (all 6 configs)
if (nrow(tcv_daily) > 0) {
  tcv_h5 <- tcv_daily[horizon == 5]

  p_ovp <- ggplot(tcv_h5, aes(x = time_id)) +
    geom_line(aes(y = obs), color = "#2166AC", linewidth = 0.6) +
    geom_line(aes(y = pred), color = "#B2182B", linewidth = 0.6, alpha = 0.7) +
    facet_wrap(~wave_lab, scales = "free", nrow = 2) +
    labs(title = "Observed (blue) vs Fitted (red) Daily Aggregate Cases",
         subtitle = "h = 5 day forecast horizon",
         x = "Forecast Day", y = "Daily Cases") +
    THEME_PUB +
    theme(legend.position = "none")
  save_pub(p_ovp, "obs_vs_pred_all.png", w_in = 14, h_in = 7)
}

# 20. Spaghetti traces (unit-level)
spaghetti_plots <- list()
for (w in WAVES_ALL) {
  obj <- tcv_list[[w]]
  if (is.null(obj)) next
  d3 <- obj[["3"]]
  if (is.null(d3) || is.null(d3$cell_preds)) next
  cp <- as.data.table(d3$cell_preds)
  if (!"area_id" %in% names(cp)) next

  # Area-level aggregate
  area_agg <- cp[, .(obs = sum(obs, na.rm = TRUE),
                      pred = sum(pred, na.rm = TRUE)), by = .(area_id, time_id)]

  p_sp <- ggplot(area_agg, aes(x = time_id, group = area_id)) +
    geom_line(aes(y = obs), color = alpha("#2166AC", 0.3), linewidth = 0.3) +
    geom_line(aes(y = pred), color = alpha("#B2182B", 0.3), linewidth = 0.3) +
    labs(title = WAVE_LABS[w], x = "Day", y = "Cases") +
    THEME_PUB +
    theme(legend.position = "none", plot.title = element_text(size = 9))
  spaghetti_plots[[w]] <- p_sp
}

if (length(spaghetti_plots) >= 2) {
  p_spaghetti <- wrap_plots(spaghetti_plots, nrow = 2) +
    plot_annotation(
      title = "Unit-Level Model Fit — Observed (blue) vs Fitted (red)",
      subtitle = "Each trace = one MODZCTA (NYC) or county (Utah)",
      theme = theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 11),
                    plot.subtitle = element_text(hjust = 0.5, color = "grey40", size = 9))
    )
  save_pub(p_spaghetti, "obs_vs_fitted_spaghetti.png", w_in = 14, h_in = 8)
}

# 21. CRPS by horizon
if (nrow(tcv_scores) > 0 && "CRPS" %in% names(tcv_scores)) {
  p_crps <- ggplot(tcv_scores, aes(x = factor(horizon), y = CRPS,
                                    color = wave_lab, group = wave_lab)) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 2.5) +
    scale_color_manual(values = WAVE_COLS_LAB, name = NULL) +
    labs(title = "CRPS by Forecast Horizon",
         subtitle = "NYC values stable; Utah larger due to absolute count magnitude",
         x = "Forecast Horizon (days)", y = "CRPS") +
    THEME_PUB
  save_pub(p_crps, "cv_crps_horizon.png", w_in = 8, h_in = 5)
}

} # end if(FALSE) — legacy CV figure block (schema mismatch)


# ═══════════════════════════════════════════════════════════════════════════════
# 22. LISA HOTSPOT CROSS-WAVE (NYC)
# ═══════════════════════════════════════════════════════════════════════════════

# The LISA maps require spatial weights and local Moran's I computation
lisa_plots <- lapply(ny_waves, function(w) {
  sp <- load_se(w)
  n_areas <- nrow(sp) %/% 2L
  sp_s <- sp[idx_spatial <= n_areas]; setorder(sp_s, idx_spatial)
  sp_u <- sp[idx_spatial >  n_areas]; setorder(sp_u, idx_spatial)

  total_eff <- data.table(
    grid_id  = sp_s$grid_id,
    total_RR = exp(sp_s$mean + sp_u$mean)
  )
  rr_sf <- merge(grid_sf_ny, total_eff, by = "grid_id", all.x = FALSE)

  # Build spatial weights
  nb <- poly2nb(rr_sf, queen = TRUE)
  lw <- nb2listw(nb, style = "W", zero.policy = TRUE)

  # Local Moran's I
  lisa <- localmoran(rr_sf$total_RR, lw, zero.policy = TRUE)
  rr_sf$lisa_I <- lisa[, 1]
  rr_sf$lisa_p <- lisa[, 5]

  # Classify clusters (FDR adjusted)
  p_adj <- p.adjust(rr_sf$lisa_p, method = "fdr")
  sig <- p_adj < 0.005
  z_rr <- scale(rr_sf$total_RR)[, 1]
  lag_z <- lag.listw(lw, z_rr, zero.policy = TRUE)

  cluster <- rep("Not significant", length(z_rr))
  cluster[sig & z_rr > 0 & lag_z > 0] <- "High-High"
  cluster[sig & z_rr < 0 & lag_z < 0] <- "Low-Low"
  cluster[sig & z_rr > 0 & lag_z < 0] <- "High-Low"
  cluster[sig & z_rr < 0 & lag_z > 0] <- "Low-High"
  rr_sf$cluster <- factor(cluster, levels = c("High-High", "Low-Low",
                                               "High-Low", "Low-High",
                                               "Not significant"))

  lisa_cols <- c("High-High" = "#B2182B", "Low-Low" = "#2166AC",
                 "High-Low" = "#EF8A62", "Low-High" = "#67A9CF",
                 "Not significant" = "grey90")

  ggplot(rr_sf) +
    geom_sf(data = modzcta_rr_sf, fill = NA, color = "grey60", linewidth = 0.25) +
    geom_sf(aes(fill = cluster), color = NA) +
    scale_fill_manual(values = lisa_cols, name = "Cluster", drop = FALSE) +
    coord_sf(datum = NA) +
    labs(title = WAVE_LABS[w]) +
    THEME_MAP +
    theme(legend.position = if (w == tail(ny_waves, 1)) "right" else "none")
})

p_lisa <- wrap_plots(lisa_plots, nrow = 1) +
  plot_layout(guides = "collect") +
  plot_annotation(
    title = "LISA Hotspot Cluster Maps — NYC (FDR-adjusted p < 0.005)",
    subtitle = "High-High (red) = significant spatial concentration of elevated risk",
    theme = theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 12),
                  plot.subtitle = element_text(hjust = 0.5, color = "grey40", size = 9))
  ) &
  theme(legend.position = "right")
save_pub(p_lisa, "ny_lisa_cross_wave.png", w_in = 14, h_in = 5.5)


# ═══════════════════════════════════════════════════════════════════════════════
# 23. DISAGGREGATION CONSERVATION
# ═══════════════════════════════════════════════════════════════════════════════

md_path <- file.path(BASE, "NY_wave1", "models", "model_data_gt.rds")
if (file.exists(md_path)) {
  md <- readRDS(md_path)
  gt <- as.data.table(md$gt_dt)
  # NOTE: legacy schema stored `md$modzcta_cases`; current schema does not.
  # The disaggregation-conservation figure relies on that field, so skip if
  # absent rather than crash. Regenerate via the disaggregation pipeline if
  # this figure becomes needed again.
  if (!is.null(md$modzcta_cases) && "geo_id" %in% names(gt)) {
    allocated <- gt[geo_id != "unknown" & geo_id != "99999",
                    .(allocated = sum(Y, na.rm = TRUE)), by = .(geo_id, time_id)]

    orig <- as.data.table(md$modzcta_cases)
    cons <- merge(orig, allocated, by = c("geo_id", "time_id"))

    p_disagg <- ggplot(cons, aes(x = observed, y = allocated)) +
      geom_abline(intercept = 0, slope = 1, linetype = "dashed",
                  color = "red", linewidth = 0.6) +
      geom_point(alpha = 0.12, size = 1, color = "steelblue") +
      coord_equal() +
      labs(title = "Disaggregation Conservation Check — NYC Wave 1",
           subtitle = "MODZCTA-day: observed vs sum of grid allocations",
           x = "Observed MODZCTA-day Cases", y = "Sum of Grid Allocations") +
      THEME_PUB
    save_pub(p_disagg, "disagg_conservation.png", w_in = 6, h_in = 6)
    rm(allocated, cons)
  } else {
    message("[disagg] skipping — model_data_gt.rds lacks $modzcta_cases (legacy schema).")
  }
  rm(md, gt)
  gc(verbose = FALSE)
}


# ═══════════════════════════════════════════════════════════════════════════════
# 24–25. OD FLOWS — Copy existing (already correct)
# ═══════════════════════════════════════════════════════════════════════════════

for (fname in c("ny_od_flows_2022_01.png", "ut_od_flows_2022_01.png")) {
  src <- file.path("thesis_latex", "figures", fname)
  file.copy(src, file.path(OUT_DIR, fname), overwrite = TRUE)
}

