# _utils.R
# Shared helper functions for the spatiotemporal risk pipeline.
# Sourced by every stage after 00_setup.

# ---- project-root anchor (reproducibility) ----
# Anchor to THIS file's location: helpers/_utils.R lives at <root>/helpers/.
# The script's path is recovered from the source() call frame so cwd does not
# matter (works from project root, pipeline/, /tmp, RStudio, knitr).
.find_script_path <- function() {
  for (i in rev(seq_along(sys.frames()))) {
    f <- sys.frames()[[i]]
    if (!is.null(f$ofile)) return(normalizePath(f$ofile, mustWork = FALSE))
  }
  # Fallback: knitr / Rmd
  if (requireNamespace("knitr", quietly = TRUE)) {
    kf <- tryCatch(knitr::current_input(dir = TRUE), error = function(e) NULL)
    if (!is.null(kf)) return(normalizePath(kf, mustWork = FALSE))
  }
  NA_character_
}
.this_file <- .find_script_path()
if (is.na(.this_file))
  stop("Could not determine helpers/_utils.R location; source() it explicitly.")
.PROJECT_ROOT <- dirname(dirname(.this_file))   # <root>/helpers/_utils.R -> <root>
if (normalizePath(getwd(), mustWork = FALSE) != normalizePath(.PROJECT_ROOT, mustWork = FALSE))
  setwd(.PROJECT_ROOT)
# When knitting an Rmd, knitr resets cwd to the Rmd's directory before every
# chunk.  Pin its root.dir to the project so all chunks resolve relative paths
# from the project root (otherwise stages create stray outputs/ in pipeline/).
if (requireNamespace("knitr", quietly = TRUE))
  knitr::opts_knit$set(root.dir = .PROJECT_ROOT)
rm(.find_script_path, .this_file)

# ---- pipeline bootstrap ----
load_pipeline <- function(base = file.path("outputs", "spatiotemporal_risk")) {
  suppressPackageStartupMessages(library(data.table))
  cfg   <- readRDS(file.path(base, "CONFIG.rds"))
  wruns <- readRDS(file.path(base, "WAVE_RUNS.rds"))

  # Wave-native pipeline: WAVE_RUNS is the canonical batch driver.
  # ALL_RUNS is a derived view (wave x month) for stages that need to iterate
  # over the constituent calendar months (e.g. 02a year-specific rasters,
  # 02b NYT slicing, 03 Advan parquet files).  All output paths use wave dirs.
  runs <- rbindlist(lapply(seq_len(nrow(wruns)), function(i) {
    months <- wruns$wave_months[[i]]
    data.table(
      state    = wruns$state[i],
      month    = months,
      wave_id  = wruns$wave_id[i],
      acs_year = wruns$acs_year[i],
      run_id   = paste0(wruns$state[i], "_", wruns$wave_id[i]),
      out_dir  = file.path(base, paste0(wruns$state[i], "_", wruns$wave_id[i]))
    )
  }))

  # Assign to caller's environment
  assign("CONFIG",    cfg,   envir = parent.frame())
  assign("ALL_RUNS",  runs,  envir = parent.frame())
  assign("WAVE_RUNS", wruns, envir = parent.frame())

  # Wave-level labels only (month_labels removed from CONFIG)
  ml <- setNames(sapply(cfg$waves, `[[`, "label"),
                 sapply(cfg$waves, `[[`, "id"))
  assign("mlabels", ml, envir = parent.frame())

  invisible(list(CONFIG = cfg, ALL_RUNS = runs, WAVE_RUNS = wruns))
}

load_spatial <- function(states = unique(WAVE_RUNS$state)) {
  sd <- list()
  for (st in states) {
    cbg <- if (st == "UT") {
      raw <- sf::st_read(CONFIG$ut_counties_path, quiet = TRUE) |>
        sf::st_transform(get_utm_crs(st))
      raw$COUNTYFP10 <- sprintf("%03d", as.integer(raw$FIPS))
      raw$NAME10     <- raw$NAME
      raw$GEOID10    <- raw$FIPS_STR
      raw
    } else {
      sf::st_read(get_cbg_path(st), quiet = TRUE) |>
        sf::st_transform(get_utm_crs(st)) |>
        (\(x) x[substr(x$GEOID, 3, 5) %in% get_metro_fips(st, 3), ])()
    }
    grid <- sf::st_read(
      file.path(CONFIG$grids_dir, sprintf("grid_100m_%s.gpkg", st)),
      quiet = TRUE) |> sf::st_transform(get_utm_crs(st))
    sd[[st]] <- list(cbg_utm = cbg, grid_sf = grid)
  }
  sd
}

# Canonical batch driver for the wave-native pipeline.
# `fn` is called as fn(wave_id, wave_months); `state` is assigned in the
# caller's environment for each iteration.
run_waves_by_state <- function(fn, wave_runs = WAVE_RUNS) {
  results <- list()
  for (st in unique(wave_runs$state)) {
    st_waves <- wave_runs[wave_runs$state == st, ]
    assign("state", st, envir = parent.frame())
    for (i in seq_len(nrow(st_waves))) {
      res <- fn(st_waves$wave_id[i], st_waves$wave_months[[i]])
      results[[st_waves$run_id[i]]] <- res
    }
  }
  invisible(results)
}

# ---- formatting ----
fmt <- function(x) format(x, big.mark = ",", scientific = FALSE)

# ---- date helpers ----
get_date_range <- function(month_code) {
  parts <- strsplit(month_code, "-")[[1]]
  year <- as.integer(parts[1])
  month <- as.integer(parts[2])
  start_date <- as.Date(sprintf("%d-%02d-01", year, month))
  end_date <- lubridate::ceiling_date(start_date, "month") - 1
  list(start = start_date, end = end_date,
       n_days = as.integer(end_date - start_date) + 1, year = year, month = month)
}

get_wave_date_range <- function(wave_months) {
  ranges <- lapply(wave_months, get_date_range)
  start <- ranges[[1]]$start
  end   <- ranges[[length(ranges)]]$end
  list(start = start, end = end,
       n_days = as.integer(end - start) + 1L,
       month_ranges = ranges)
}

# ---- geography helpers ----
standardize_fips <- function(fips) {
  fips_clean <- trimws(as.character(fips))
  fips_int <- suppressWarnings(as.integer(fips_clean))
  result <- sprintf("%05d", fips_int)
  result[is.na(fips_int)] <- NA_character_
  result
}

standardize_cbg <- function(x) {
  # Force CBG GEOIDs into consistent 12-character zero-padded strings.
  # Handles: integer64, numeric (.0 suffix), scientific notation, whitespace.
  if (inherits(x, "integer64")) {
    x <- format(x, scientific = FALSE)
  } else {
    x <- as.character(x)
  }
  x <- trimws(x)
  # Handle scientific notation
  sci_idx <- grepl("[eE]", x) & !is.na(x)
  if (any(sci_idx)) {
    x[sci_idx] <- format(as.numeric(x[sci_idx]), scientific = FALSE)
    x[sci_idx] <- trimws(x[sci_idx])
  }
  # Strip .0 decimal artifacts from float conversion
  x <- sub("\\.0*$", "", x)
  # Zero-pad to 12 digits
  short_idx <- !is.na(x) & nchar(x) < 12 & grepl("^[0-9]+$", x)
  if (any(short_idx)) {
    x[short_idx] <- formatC(as.numeric(x[short_idx]), width = 12, format = "d", flag = "0")
  }
  x
}

get_utm_crs <- function(state) {
  if (state == "NY") CONFIG$ny_utm_crs else CONFIG$ut_utm_crs
}

get_metro_fips <- function(state, digits = 5) {
  if (digits == 5) {
    if (state == "NY") CONFIG$nyc_fips_5 else CONFIG$ut_fips_5
  } else {
    if (state == "NY") CONFIG$nyc_fips_3 else CONFIG$ut_fips_3
  }
}

get_cbg_path <- function(state) {
  if (state == "NY") CONFIG$ny_cbg_path else CONFIG$ut_cbg_path
}

get_landuse_path <- function(state, year = NULL) {
  if (state == "NY") return(CONFIG$ny_landuse_path)
  # UT: prefer year-matched LIR parcel GPKG (statewide, parcel polygons).
  # Fall back to county-level summary CSV when GPKG missing.
  if (!is.null(year)) {
    gpkg <- here::here(sprintf("data/landuse/utah_parcels_lir/Parcels_LIR_%s.gpkg",
                               as.character(year)))
    if (file.exists(gpkg)) return(gpkg)
  }
  CONFIG$ut_landuse_summary
}

get_buildings_path <- function(state) {
  if (state == "NY") CONFIG$ny_buildings_path else CONFIG$ut_buildings_path
}

# Per-POI footprint area (m^2) via spatial join to MS Building Footprints.
# Cached at outputs/.../_cache/poi_area_<state>.feather, keyed by persistent_id,
# so the heavy GeoJSON read happens once per state across all waves.
#
# Returns a data.table(persistent_id, poi_area_sqm); rows with NA lon/lat
# or no containing polygon get NA area.
build_poi_area_cache <- function(poi_dt, state, force = FALSE) {
  cache_dir  <- file.path(CONFIG$output_dir, "_cache")
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  cache_file <- file.path(cache_dir, sprintf("poi_area_%s.feather", state))

  if (file.exists(cache_file) && !isTRUE(force)) {
    cached <- as.data.table(arrow::read_feather(cache_file))
    if (all(unique(poi_dt$persistent_id) %in% cached$persistent_id))
      return(cached)
  }

  uniq <- unique(poi_dt[!is.na(poi_lon) & !is.na(poi_lat),
                        .(persistent_id, poi_lon, poi_lat)])
  utm_crs <- get_utm_crs(state)
  bld <- sf::st_read(get_buildings_path(state), quiet = TRUE)
  bld <- sf::st_transform(bld, utm_crs)
  bld$.bld_area_sqm <- as.numeric(sf::st_area(bld))
  bld <- bld[, ".bld_area_sqm"]

  pts <- sf::st_as_sf(uniq, coords = c("poi_lon", "poi_lat"), crs = 4326)
  pts <- sf::st_transform(pts, utm_crs)
  joined <- sf::st_join(pts, bld, join = sf::st_within, left = TRUE)

  out <- data.table::data.table(
    persistent_id = joined$persistent_id,
    poi_area_sqm  = joined$.bld_area_sqm)
  out <- unique(out, by = "persistent_id")
  arrow::write_feather(out, cache_file)
  out
}

get_state_prefix <- function(state) {
  if (state == "NY") "36" else "49"
}

# ---- distance-decay kernel fitting (Meyer & Held 2014; Gonzalez 2008) ----
# Compares three forms by AIC on Poisson-binned OD distances:
#   power-law:           f(d) = (1 + d/sigma)^(-alpha)
#   truncated power-law: f(d) = (1 + d/sigma)^(-alpha) * exp(-d/kappa)
#   Gaussian:            f(d) = exp(-d^2 / (2 sigma^2))
# Truncated form (Gonzalez et al. 2008, Brockmann 2006) typically wins for
# COVID-era OD: kappa captures wave-specific suppression of long trips.
fit_distance_kernel <- function(distance_m, weights, n_bins = 25) {
  ok <- is.finite(distance_m) & is.finite(weights) &
        distance_m > 0 & weights > 0
  d <- distance_m[ok]; w <- weights[ok]
  if (length(d) < 50) {
    return(list(form = "powerlaw", sigma = 1000, alpha = 2, kappa = NA_real_,
                aic_pl = NA_real_, aic_tpl = NA_real_, aic_gauss = NA_real_,
                n = length(d), note = "insufficient data, using defaults"))
  }
  brks <- exp(seq(log(max(min(d), 50)), log(max(d) + 1), length.out = n_bins + 1))
  bin <- findInterval(d, brks, all.inside = TRUE)
  agg <- tapply(w, bin, sum)
  mid <- (brks[-1] + brks[-length(brks)]) / 2
  y <- as.numeric(agg); x <- mid[as.integer(names(agg))]
  area <- diff(brks)[as.integer(names(agg))]
  # Poisson nLL of binned counts vs kernel (areas absorb 2*pi*r)
  nll_pl <- function(par) {
    sg <- exp(par[1]); al <- exp(par[2])
    mu <- (1 + x / sg)^(-al) * area
    mu <- mu * (sum(y) / sum(mu))
    -sum(stats::dpois(round(y), mu, log = TRUE))
  }
  nll_tpl <- function(par) {
    sg <- exp(par[1]); al <- exp(par[2]); kp <- exp(par[3])
    mu <- (1 + x / sg)^(-al) * exp(-x / kp) * area
    mu <- mu * (sum(y) / sum(mu))
    -sum(stats::dpois(round(y), mu, log = TRUE))
  }
  nll_gauss <- function(par) {
    sg <- exp(par[1])
    mu <- exp(-x^2 / (2 * sg^2)) * area
    mu <- mu * (sum(y) / sum(mu))
    -sum(stats::dpois(round(y), mu, log = TRUE))
  }
  fit_pl <- tryCatch(stats::optim(c(log(2000), log(2)), nll_pl,
                                  method = "Nelder-Mead"),
                     error = function(e) NULL)
  fit_tpl <- tryCatch(stats::optim(c(log(2000), log(2), log(20000)), nll_tpl,
                                   method = "Nelder-Mead"),
                      error = function(e) NULL)
  fit_g  <- tryCatch(stats::optim(log(3000), nll_gauss,
                                  method = "Brent",
                                  lower = log(50), upper = log(1e6)),
                     error = function(e) NULL)
  pl_ok  <- !is.null(fit_pl)  && is.finite(fit_pl$value)
  tpl_ok <- !is.null(fit_tpl) && is.finite(fit_tpl$value)
  g_ok   <- !is.null(fit_g)   && is.finite(fit_g$value)
  out <- list(
    sigma       = if (pl_ok)  exp(fit_pl$par[1])  else 1000,
    alpha       = if (pl_ok)  exp(fit_pl$par[2])  else 2,
    sigma_tpl   = if (tpl_ok) exp(fit_tpl$par[1]) else NA_real_,
    alpha_tpl   = if (tpl_ok) exp(fit_tpl$par[2]) else NA_real_,
    kappa       = if (tpl_ok) exp(fit_tpl$par[3]) else NA_real_,
    sigma_gauss = if (g_ok)   exp(fit_g$par)      else NA_real_,
    aic_pl    = if (pl_ok)  2 * 2 + 2 * fit_pl$value  else NA_real_,
    aic_tpl   = if (tpl_ok) 2 * 3 + 2 * fit_tpl$value else NA_real_,
    aic_gauss = if (g_ok)   2 * 1 + 2 * fit_g$value   else NA_real_,
    n = length(d)
  )
  aics <- c(powerlaw = out$aic_pl, tpowerlaw = out$aic_tpl, gauss = out$aic_gauss)
  aics <- aics[is.finite(aics)]
  out$form <- if (length(aics)) names(aics)[which.min(aics)] else "powerlaw"
  out
}

get_road_paths <- function(state) {
  if (state == "NY") {
    file.path("data", "transport", "newyork_roads2", "tl_2023_36_prisecroads.shp")
  } else {
    # Utah: per-county TIGER/Line shapefiles
    list.files(CONFIG$ut_roads_dir, pattern = "\\.shp$", full.names = TRUE)
  }
}

# Wave-level output dir (canonical; per-month staging dirs are deprecated)
get_wave_dir <- function(state, wave_id) {
  file.path(CONFIG$output_dir, paste0(state, "_", wave_id))
}

# Time-weighted population for a wave: averages year-specific population CSVs
# weighted by the fraction of wave months falling in each calendar year.
# Returns data.table(grid_id, population, E, log_E).
#
# `wave_months` is a character vector of "YYYY-MM" strings. For each unique
# year, reads `population_{state}_{year}.csv` (produced by 02a) and blends.
# If a needed yearly file is missing, falls back to the closest available year.
population_for_wave <- function(state, wave_months) {
  yrs <- substr(wave_months, 1, 4)
  w   <- prop.table(table(yrs))   # weights per year, sum to 1

  pop_path <- function(y) file.path(CONFIG$grids_dir,
                                    sprintf("population_%s_%s.csv", state, y))

  # Resolve missing years to nearest available
  available <- list.files(CONFIG$grids_dir,
                          pattern = sprintf("^population_%s_\\d{4}\\.csv$", state))
  avail_yrs <- as.integer(sub(".*_(\\d{4})\\.csv$", "\\1", available))
  if (length(avail_yrs) == 0L)
    stop(sprintf("No population_%s_*.csv files in %s", state, CONFIG$grids_dir))

  resolve_year <- function(y) {
    yi <- as.integer(y)
    if (yi %in% avail_yrs) return(as.character(yi))
    nearest <- avail_yrs[which.min(abs(avail_yrs - yi))]
    warning(sprintf("population_%s_%s.csv missing; using %s", state, y, nearest))
    as.character(nearest)
  }

  parts <- lapply(names(w), function(y) {
    dt <- data.table::fread(pop_path(resolve_year(y)))
    dt[, .(grid_id, population = population * w[[y]])]
  })

  blended <- data.table::rbindlist(parts)[
    , .(population = sum(population)), by = grid_id]
  blended[, E      := pmax(round(population, 2), 0.1)]
  blended[, log_E  := log(E)]
  blended[, population := round(population, 2)]
  blended[]
}

# ---- visualization (loaded by all stages that plot) ----

# Map theme: lon/lat axes, clean panels, bold titles
theme_map <- function(base_size = 11) {
  ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(
      axis.title       = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_line(color = "grey90", linewidth = 0.3),
      panel.grid.minor = ggplot2::element_blank(),
      plot.title       = ggplot2::element_text(face = "bold", hjust = 0.5, size = base_size + 1),
      plot.subtitle    = ggplot2::element_text(color = "grey40", hjust = 0.5),
      legend.position  = "right",
      strip.background = ggplot2::element_rect(fill = "gray90"),
      strip.text       = ggplot2::element_text(face = "bold")
    )
}
# Keep old name as alias so existing code doesn't break
theme_ts <- theme_map

# Scale bar (bottom-right) + north arrow (top-left)
map_extras <- function(scale_loc = "br", arrow_loc = "tl") {
  list(
    ggspatial::annotation_scale(location = scale_loc, width_hint = 0.15,
                                style = "ticks", line_width = 0.8),
    ggspatial::annotation_north_arrow(location = arrow_loc, which_north = "true",
                                      height = grid::unit(1, "cm"),
                                      width = grid::unit(0.8, "cm"),
                                      style = ggspatial::north_arrow_minimal())
  )
}

# Faded CBG base layer (call once, add to ggplot with +)
cbg_base_layer <- function(cbg_sf = NULL, state = "NY", fill = "grey95", color = "grey70",
                           linewidth = 0.2, alpha = 0.4) {
  if (is.null(cbg_sf)) {
    if (state == "UT") {
      cbg_sf <- sf::st_read(CONFIG$ut_counties_path, quiet = TRUE) |>
        sf::st_transform(get_utm_crs(state))
      linewidth <- 0.4
      color <- "grey50"
    } else {
      cbg_sf <- sf::st_read(get_cbg_path(state), quiet = TRUE) |>
        sf::st_transform(get_utm_crs(state))
      cbg_sf <- cbg_sf[substr(cbg_sf$GEOID, 3, 5) %in% get_metro_fips(state, 3), ]
    }
  }
  ggplot2::geom_sf(data = cbg_sf, fill = fill, color = color,
                   linewidth = linewidth, alpha = alpha, inherit.aes = FALSE)
}

# For UT maps: aggregate grid-level values to county polygons (291k cells are sub-pixel at state scale)
# grid_sf must have 'fips_5' column; value_col is the column to aggregate; fun is the aggregation function
ut_county_agg_sf <- function(grid_sf, value_col, fun = mean, na.rm = TRUE) {
  counties <- sf::st_read(CONFIG$ut_counties_path, quiet = TRUE) |>
    sf::st_transform(get_utm_crs("UT"))
  counties$GEOID10 <- counties$FIPS_STR
  gdt <- data.table::as.data.table(sf::st_drop_geometry(grid_sf))
  agg <- gdt[!is.na(get(value_col)), .(agg_val = fun(get(value_col), na.rm = na.rm)), by = fips_5]
  data.table::setnames(agg, "agg_val", value_col)
  merge(counties, agg, by.x = "GEOID10", by.y = "fips_5", all.x = TRUE)
}

# ---- shared constants ----
EXPOSURE_FLOOR <- 0.05
OFFSET_MIN     <- -4
OFFSET_MAX     <- 20      # generous safety cap; was 6, which silently clipped UT county-level offsets up to ~15
# State-specific tighter cap for NY: NY offsets stay in ~[-4, 6], so clipping at
# 6 trims pathological tail values that destabilise the likelihood (extreme
# E_inla = exp(offset)) without affecting realistic NY predictors. UT keeps 20.
N_RHO_BINS     <- 20
# May 13 2026: dropped 300 -> 40 after LinMob_ST_IV NY wave 2 stalled
# >30h with K=300 (Type IV block = 300 * 90 = 27,000 latent vars, on top
# of 96,964-dim BYM2). K=40 -> 3,600 latent vars; expected runtime
# ~4-6h for LinMob_ST_IV / NLMob_ST_IV on NY full-domain.
N_ST_CLUSTERS  <- 40
# Visualization-only flag: month-boundary despiking is a legacy patch
# from the monthly-Advan era. Weekly Advan inputs have no month-boundary
# step artifact, so the despike function would either no-op or, worse,
# flatten genuine signal. Set to TRUE only if reverting to monthly inputs.
USE_DESPIKE   <- FALSE
                        # Larger cluster count -> each cluster aggregates fewer
                        # observations -> healthier likelihood/prior balance and
                        # fewer NaN/INF rescues at the linear predictor.

hyper_bym2 <- list(prec = list(prior = "pc.prec", param = c(1, 0.01)),
                   phi  = list(prior = "pc",      param = c(0.5, 0.5)))

POTENTIAL_COVARIATES <- c(
  "ses_contrast_std", "pct_older_std", "pct_younger_std", "pct_employed_std",
  "pct_unemployed_std", "building_density_std", "landuse_commercial",
  "log_activity_density_std", "hourly_concentration_std",
  "daily_visit_cv_std", "weekday_dominance_std", "log_repeat_ratio_std",
  "log_contact_density_hourly_std", "log_daytime_pop_std",
  "is_weekend", "commercial_x_weekend", "activity_x_weekend",
  "avg_category_risk_std", "log_distance_from_home_std", "inf_pressure_std"
)

# ---- stage-10 shared helpers ----
# Common functions used by 12_visualization.Rmd (parameterized: NY / UT)

safe_save <- function(fig_dir, fname, p, w = 10, h = 8) {
  grDevices::png(file.path(fig_dir, fname),
                 width = w, height = h, units = "in", res = 300, bg = "white")
  on.exit(grDevices::dev.off(), add = TRUE)
  print(p)
  invisible(NULL)
}

main_annotation <- function(title_text, month_label)
  plot_annotation(title = paste(title_text, "—", month_label),
                  theme = theme(plot.title = element_text(
                    face = "bold", hjust = 0.5, size = 14)))

load_wave <- function(wave_id, state, grid_sf0) {
  out_dir <- file.path(CONFIG$output_dir, paste0(state, "_", wave_id))
  tbl_dir <- file.path(out_dir, "tables")
  res09b  <- readRDS(file.path(out_dir, "models", "model_results_unified.rds"))
  list(
    wave_id = wave_id, grid_sf = grid_sf0, out_dir = out_dir,
    fig_dir    = file.path(out_dir, "figures"),
    best_model = res09b$best_model,
    res09b     = res09b,
    label      = res09b$wave_label %||% wave_id,
    comp_09b   = data.table::fread(file.path(tbl_dir, "model_comparison_unified.csv")),
    fixed_09b  = data.table::fread(file.path(tbl_dir, "unified_fixed_effects.csv")),
    spatial_09b = data.table::fread(file.path(tbl_dir, "unified_spatial_effects.csv")),
    temporal_09b = data.table::fread(file.path(tbl_dir, "unified_temporal_effects.csv")),
    fitted_09b = data.table::fread(file.path(tbl_dir, "unified_fitted_values.csv"))
  )
}

# ---- batch-reporting redistribution (UT only, Nash et al. 2023) ----
# Detects zero-runs in daily case counts and redistributes the subsequent
# dump-day total back across the gap using growth-adjusted weights.
# Conservation: sum of cases per county is exactly preserved.

redistribute_batch_cases <- function(dt, id_col = "fips", cases_col = "new_cases",
                                     min_gap = 2, threshold_frac = 0.02) {
  setorderv(dt, c(id_col, "date"))
  n_redist <- 0L
  total_moved <- 0

  for (gid in unique(dt[[id_col]])) {
    idx <- which(dt[[id_col]] == gid)
    x   <- dt[[cases_col]][idx]
    n   <- length(x)
    if (n < 3 || sum(x) == 0) next

    med_x  <- median(x[x > 0], na.rm = TRUE)
    thresh <- max(1, threshold_frac * med_x)
    is_low <- x <= thresh

    # Walk through and find gap -> dump patterns
    i <- 1L
    while (i <= n) {
      if (is_low[i]) {
        # Start of a potential gap
        gap_start <- i
        while (i <= n && is_low[i]) i <- i + 1L
        gap_len <- i - gap_start
        if (gap_len >= min_gap && i <= n && x[i] > thresh) {
          # i is now the dump day — include it in the window
          window <- gap_start:i
          window_total <- sum(x[window])
          if (window_total > 0) {
            # Growth-adjusted weights: estimate local trend from flanking days
            pre_val  <- if (gap_start > 1) max(x[gap_start - 1], 1) else max(window_total / length(window), 1)
            post_val <- if (i < n) max(x[min(i + 1, n)], 1) else pre_val
            # Linear interpolation weights across the window
            w_len <- length(window)
            weights <- seq(pre_val, post_val, length.out = w_len)
            weights <- weights / sum(weights)
            # Distribute as integers, conserve total exactly
            alloc <- floor(weights * window_total)
            remainder <- as.integer(window_total - sum(alloc))
            # Assign remainder to positions with largest fractional parts
            if (remainder > 0) {
              frac <- (weights * window_total) - alloc
              top_idx <- order(frac, decreasing = TRUE)[seq_len(remainder)]
              alloc[top_idx] <- alloc[top_idx] + 1L
            }
            x[window] <- as.integer(alloc)
            n_redist <- n_redist + 1L
            total_moved <- total_moved + window_total - alloc[w_len]
          }
          i <- i + 1L  # move past dump day
        }
        # else gap too short; i already advanced past it
      } else {
        i <- i + 1L
      }
    }
    dt[[cases_col]][idx] <- x
  }
  dt
}

# ---- preprocessing helpers (shared by 09 and 09b) ----

add_spatial_orthogonalization <- function(df, grid_sf) {
  nb <- tryCatch(suppressWarnings(spdep::poly2nb(grid_sf, queen = TRUE, snap = 1)), error = function(e) NULL)
  lw <- if (!is.null(nb)) tryCatch(spdep::nb2listw(nb, style = "W", zero.policy = TRUE), error = function(e) NULL) else NULL
  if (is.null(lw)) return(df)
  dt <- as.data.table(df)
  ga <- dt[, .(avg = mean(log_rho, na.rm = TRUE)), by = idx_spatial][order(idx_spatial)]
  lv <- spdep::lag.listw(lw, ga$avg)
  ga[, resid := residuals(lm(avg ~ lv, na.action = na.exclude))]
  ga[, spatial_lag := lv]
  dt <- merge(dt, ga[, .(idx_spatial, avg, resid, spatial_lag)], by = "idx_spatial", all.x = TRUE)
  dt[, log_rho_resid := resid + (log_rho - avg)]
  dt[, mob_spatial_lag_std := scale(spatial_lag)[, 1]]
  dt[, c("avg", "resid", "spatial_lag") := NULL]
  as.data.frame(dt)
}

add_temporal_centering <- function(df) {
  rho_col <- if ("log_rho_resid" %in% names(df)) "log_rho_resid" else "log_rho"
  dt <- as.data.table(df)
  dt[, rho_daily_mean := mean(get(rho_col), na.rm = TRUE), by = idx_temporal]
  if (!"mob_day_mean" %in% names(dt)) dt[, mob_day_mean := rho_daily_mean]
  dt[, (rho_col) := get(rho_col) - rho_daily_mean]
  dt[, rho_daily_mean := NULL]
  as.data.frame(dt)
}

add_temporal_lags <- function(df, lags = c(3, 5, 7, 14, 21)) {
  dt <- as.data.table(df); setorder(dt, idx_spatial, idx_temporal)
  for (lag_k in lags) {
    lag_col <- paste0("mob_lag", lag_k)
    dt[, (lag_col) := shift(log_rho, lag_k, type = "lag"), by = idx_spatial]
    # LOCF: fill leading NAs with earliest available value (avoids jump artifacts)
    dt[is.na(get(lag_col)), (lag_col) := log_rho[1], by = idx_spatial]
  }
  # Negative-control: future mobility (lead by 7d). Should have ~zero
  # coefficient in NLMob; significant beta indicates over-fitting / random-walk
  # leakage. Filled at the right boundary by LOCF on the trailing value.
  dt[, mob_lag_neg := shift(log_rho, 7L, type = "lead"), by = idx_spatial]
  dt[is.na(mob_lag_neg), mob_lag_neg := log_rho[.N], by = idx_spatial]
  as.data.frame(dt)
}

add_rho_bins <- function(df, n_bins = N_RHO_BINS) {
  brks <- seq(min(df$log_rho, na.rm = TRUE) - 0.01,
              max(df$log_rho, na.rm = TRUE) + 0.01, length.out = n_bins + 1)
  df$log_rho_idx <- as.integer(cut(df$log_rho, breaks = brks, include.lowest = TRUE))
  df$log_rho_idx[is.na(df$log_rho_idx)] <- 1L
  df
}

# ---- Power-law adjacency-order baseline rho ----  REMOVED (Base_PL dropped)
# Was an internal scaffolding model never written up in the thesis. The
# six model variants are: Base, Temp, LinMob, NLMob, LinMob_ST_IV, NLMob_ST_IV.


add_st_clusters <- function(df, grid_sf, out_dir, n_clusters = N_ST_CLUSTERS) {
  n_grids <- max(df$idx_spatial)
  k <- min(n_clusters, n_grids)
  if (k >= n_grids) {
    df$idx_st_cluster <- df$idx_spatial
    # Each grid is its own cluster — build graph from spatial adjacency
    suppressWarnings(nb_cl <- spdep::poly2nb(grid_sf, queen = TRUE, snap = 1))
    card_nb <- sapply(nb_cl, function(x) sum(x > 0))
    if (any(card_nb == 0)) {
      coords_nb <- sf::st_coordinates(suppressWarnings(sf::st_centroid(grid_sf)))
      knn_nb <- spdep::knearneigh(coords_nb, k = min(8, nrow(grid_sf) - 1))
      nb_cl <- spdep::knn2nb(knn_nb, sym = TRUE)
    }
    cluster_graph <- file.path(out_dir, "models", "cluster_adj.graph")
    spdep::nb2INLA(cluster_graph, nb_cl)
    return(list(df = df, cluster_graph = cluster_graph))
  }
  coords <- sf::st_coordinates(suppressWarnings(sf::st_centroid(grid_sf)))
  set.seed(42); cl <- kmeans(coords, centers = k, nstart = 5)$cluster
  lookup <- data.table(grid_id = grid_sf$grid_id, idx_st_cluster = as.integer(cl))
  dt <- as.data.table(df)
  if ("idx_st_cluster" %in% names(dt)) dt[, idx_st_cluster := NULL]
  dt <- merge(dt, lookup, by = "grid_id", all.x = TRUE)
  grid_cl <- data.table(grid_id = grid_sf$grid_id, cl = as.integer(cl))
  suppressWarnings(nb_grid <- spdep::poly2nb(grid_sf, queen = TRUE, snap = 1))
  # Fix islands: if any cell has zero neighbors (e.g. sparse UT grid), use KNN fallback
  card_nb <- sapply(nb_grid, function(x) sum(x > 0))
  if (any(card_nb == 0)) {
    coords_nb <- sf::st_coordinates(suppressWarnings(sf::st_centroid(grid_sf)))
    knn_nb <- spdep::knearneigh(coords_nb, k = min(8, nrow(grid_sf) - 1))
    nb_grid <- spdep::knn2nb(knn_nb, sym = TRUE)
  }
  edges <- unique(rbindlist(lapply(seq_along(nb_grid), function(i) {
    ni <- nb_grid[[i]]; ni <- ni[ni > 0]
    if (length(ni) == 0) return(NULL)
    data.table(a = grid_cl$cl[i], b = grid_cl$cl[ni])
  })))
  edges <- edges[a != b]
  nb_cl <- lapply(seq_len(k), function(i) {
    nbrs <- sort(unique(edges[a == i, b]))
    if (length(nbrs) == 0) 0L else nbrs
  })
  class(nb_cl) <- "nb"; attr(nb_cl, "region.id") <- as.character(seq_len(k))
  # Fix isolated clusters: connect to nearest cluster by centroid distance
  cl_centers <- t(sapply(seq_len(k), function(i) colMeans(coords[cl == i, , drop = FALSE])))
  iso <- which(sapply(nb_cl, function(x) length(x) == 1 && x[1] == 0L))
  if (length(iso) > 0) {
    others <- setdiff(seq_len(k), iso)
    for (ci in iso) {
      if (length(others) == 0) next
      dists <- sqrt(rowSums((cl_centers[others, , drop = FALSE] - matrix(cl_centers[ci, ], nrow = length(others), ncol = 2, byrow = TRUE))^2))
      nearest <- others[which.min(dists)]
      nb_cl[[ci]] <- as.integer(nearest)
      nb_cl[[nearest]] <- sort(unique(c(nb_cl[[nearest]][nb_cl[[nearest]] > 0L], as.integer(ci))))
    }
    class(nb_cl) <- "nb"; attr(nb_cl, "region.id") <- as.character(seq_len(k))
  }
  cluster_graph <- file.path(out_dir, "models", "cluster_adj.graph")
  spdep::nb2INLA(cluster_graph, nb_cl)
  list(df = as.data.frame(dt), cluster_graph = cluster_graph)
}

stabilize_offsets <- function(df) {
  df$offset_total <- pmax(pmin(df$offset_total, OFFSET_MAX), OFFSET_MIN)
  df$E_inla <- exp(df$offset_total)
  efl <- max(quantile(df$E_inla, 0.005, na.rm = TRUE), EXPOSURE_FLOOR)
  df$E_inla <- pmax(df$E_inla, efl)
  # Option 2: log_rho is no longer in the offset, so all variants share the
  # same exposure offset; they differ only in which mobility covariate (none
  # vs log_rho_resid) enters the linear predictor. E_nomob is kept as an
  # alias for backward compatibility with get_E() and downstream extractors.
  df$E_nomob <- df$E_inla
  df
}

select_best_lag <- function(df) {
  lag_candidates <- paste0("mob_lag", c(3, 5, 7, 14, 21))
  lag_candidates <- lag_candidates[lag_candidates %in% names(df)]
  lag_candidates <- lag_candidates[sapply(lag_candidates, function(v) sd(df[[v]], na.rm = TRUE) > 0.01)]
  if (length(lag_candidates) == 0) return(NULL)
  lag_cors <- sapply(lag_candidates, function(v) {
    ok <- !is.na(df$Y) & !is.na(df[[v]]) & df$Y > 0
    if (sum(ok) > 100) abs(cor(df$Y[ok], df[[v]][ok])) else 0
  })
  lag_candidates[which.max(lag_cors)]
}

# Multi-lag selector (AOAS Slater 2025 §3.2): retains all lags whose marginal
# |cor(Y, mob_lag_k)| is at least 50% of the strongest lag, so a wave with a
# clear single peak still gets one term but a wave with a broad peak gets the
# full multi-lag block. Returns a character vector (possibly length 1).
select_mob_lags <- function(df, rel_threshold = 0.5, abs_floor = 0.02) {
  lag_candidates <- paste0("mob_lag", c(3, 5, 7, 14, 21))
  lag_candidates <- lag_candidates[lag_candidates %in% names(df)]
  lag_candidates <- lag_candidates[sapply(lag_candidates, function(v) sd(df[[v]], na.rm = TRUE) > 0.01)]
  if (length(lag_candidates) == 0) return(character(0))
  lag_cors <- sapply(lag_candidates, function(v) {
    ok <- !is.na(df$Y) & !is.na(df[[v]]) & df$Y > 0
    if (sum(ok) > 100) abs(cor(df$Y[ok], df[[v]][ok])) else 0
  })
  if (max(lag_cors) < abs_floor) return(character(0))
  lag_candidates[lag_cors >= rel_threshold * max(lag_cors)]
}

# ---- model helpers (shared by 09 and 09b) ----

rename_quants <- function(dt) {
  nms <- names(dt)
  nms[nms == "0.025quant"] <- "q025"; nms[nms == "0.5quant"] <- "q50"; nms[nms == "0.975quant"] <- "q975"
  setnames(dt, names(dt), nms)
}

get_E <- function(mn, df) {
  if (mn %in% c("Base", "Temp")) df$E_nomob
  else                           df$E_inla
}

# Model labels (self-descriptive names).  Order = increasing complexity;
# also defines the canonical row order in the comparison table.
# Both ST_III ablation variants (NLMob x ST_III and LinMob x ST_III) were
# dropped May 2026: they only existed as warmstart seeds for the Type IV
# fits and were never cited; their removal halves the comparison-table cost.
# Key NLMob_ST_IV (formerly ST_IV) renamed for self-documentation.
MODEL_LABELS <- c(
  Base          = "Baseline (BYM2+RW2, no mobility)",
  Temp          = "Temporal-only (RW2, no spatial RE)",
  LinMob        = "Linear mobility (BYM2+RW2+rho)",
  NLMob         = "Nonlinear mobility (BYM2+RW1+rw1(rho))",
  LinMob_ST_IV  = "LinMob + Knorr-Held Type IV (ICAR x RW1)",
  NLMob_ST_IV   = "NLMob + Knorr-Held Type IV (ICAR x RW1)"
)

# Models that should run with CPO enabled.
# NOTE: ST_IV CPO at the full-domain scale (~1.4k grids x ~120 days = ~170k
# observations per wave with augmented adjacency) repeatedly OOM-kills
# RStudio + macOS even at num.threads="1:1" (3 crashes, last after 9h of
# slow swap thrash). LCPO is a nice-to-have; WAIC + DIC are the canonical
# Bayesian model-comparison metrics for the thesis and are computed for
# every model. Disable CPO globally for the main fits — the CV path
# (11a/11b) computes out-of-sample predictive accuracy directly via
# fold-wise log-scores, which is the methodologically sounder substitute
# for LCPO anyway.
CPO_MODELS <- character(0)

# Mobility-aware models -- receive the selected mob lags + negative control.
MOB_MODELS <- c("LinMob", "NLMob", "LinMob_ST_IV", "NLMob_ST_IV")

# Display names: identity map (the on-disk keys are now self-documenting).
MODEL_DISPLAY <- c(
  Base         = "Base",
  Temp         = "Temp",
  LinMob       = "LinMob",
  NLMob        = "NLMob",
  LinMob_ST_IV = "LinMob_ST_IV",
  NLMob_ST_IV  = "NLMob_ST_IV"
)
display_model <- function(x) {
  out <- unname(MODEL_DISPLAY[as.character(x)])
  ifelse(is.na(out), as.character(x), out)
}

detect_flags <- function(df) {
  has <- function(x) x %in% names(df) && isTRUE(sd(df[[x]], na.rm = TRUE) > 0.01)
  list(rho = "log_rho" %in% names(df), contact = has("contact_intensity_std"),
       exposure_import = has("exposure_import_std"), spatial_lag = has("mob_spatial_lag_std"))
}

build_formulas <- function(covs, graph, flags, dispersion = 1.1, cluster_graph = NULL, n_grids = Inf,
                           state = NULL, wave_id = NULL) {
  fx <- paste(covs, collapse = " + ")
  small_domain <- is.finite(n_grids) && n_grids < 50

  # NY Wave 3 boundary-fit fix (May 15 2026):
  # The default PC priors allowed the BYM2 spatial and ST_IV cluster
  # precisions to drift to ~e^25-30 (data-driven near-zero variance under
  # Omicron homogeneity), making INLA's Newton step explore degenerate
  # corners that take >5 min per fn evaluation. Tighten the precision
  # PC priors to keep the optimizer in a well-conditioned region.
  tight <- isTRUE(state == "NY") && isTRUE(wave_id == "wave3")

  if (small_domain) {
    bym2 <- sprintf("f(idx_spatial,model='bym2',graph='%s',scale.model=TRUE,hyper=list(prec=list(prior='pc.prec',param=c(2,0.01)),phi=list(prior='pc',param=c(0.5,0.8))))", graph)
  } else if (tight) {
    # P(sigma_spatial > 0.3) = 0.01 ; strong shrinkage prior on phi toward structured.
    bym2 <- sprintf("f(idx_spatial,model='bym2',graph='%s',scale.model=TRUE,hyper=list(prec=list(prior='pc.prec',param=c(0.3,0.01)),phi=list(prior='pc',param=c(0.5,0.5))))", graph)
  } else {
    bym2 <- sprintf("f(idx_spatial,model='bym2',graph='%s',scale.model=TRUE,hyper=hyper_bym2)", graph)
  }
  rw_u <- if (dispersion < 1.5) 1 else if (dispersion < 3) 0.75 else 0.5
  if (tight) rw_u <- min(rw_u, 0.3)
  rho_u <- if (tight) 0.3 else 0.5
  rw2 <- sprintf("f(idx_temporal,model='rw2',constr=TRUE,scale.model=TRUE,hyper=list(prec=list(prior='pc.prec',param=c(%s,0.01))))", rw_u)
  rw1_t <- sprintf("f(idx_temporal,model='rw1',constr=TRUE,scale.model=TRUE,hyper=list(prec=list(prior='pc.prec',param=c(%s,0.01))))", rw_u)
  rw1_rho <- sprintf("f(log_rho_idx,model='rw1',constr=TRUE,scale.model=TRUE,hyper=list(prec=list(prior='pc.prec',param=c(%s,0.01))))", rho_u)
  rv <- if (flags$spatial_lag) "log_rho_resid" else "log_rho"

  fms <- list(Base = as.formula(paste(c("Y ~ 1", fx, bym2, rw2), collapse = " + ")))
  if (flags$rho) {
    lm <- c("Y ~ 1", fx, bym2, rw2, rv)
    if (flags$contact) lm <- c(lm, "contact_intensity_std")
    if (flags$spatial_lag) lm <- c(lm, "mob_spatial_lag_std")
    fms$LinMob <- as.formula(paste(lm, collapse = " + "))
    nl <- c("Y ~ 1", fx, bym2, rw1_t, rw1_rho)
    if (flags$contact) nl <- c(nl, "contact_intensity_std")
    if (flags$spatial_lag) nl <- c(nl, "mob_spatial_lag_std")
    fms$NLMob <- as.formula(paste(nl, collapse = " + "))
    if (flags$contact && !is.null(cluster_graph)) {
      st_u <- if (small_domain) 0.5
              else if (dispersion > 5) 0.05
              else if (dispersion > 3) 0.08
              else 0.15
      # NY Wave 3 ST_IV fix (May 16 2026): the dispersion-driven st_u=0.05
      # creates a prior-likelihood conflict (Pr(sigma_cluster>0.05)=0.01 →
      # risk-ratio <1.05) that Omicron's spatial heterogeneity fights hard,
      # producing a degenerate posterior ridge where Newton stalls
      # (maxld=+32M, fn=15 in 3.5h). Override to 0.3 so the cluster term
      # can express variance commensurate with the BYM2 prior.
      if (tight) st_u <- 0.3
      st_diag <- if (small_domain) 0.5 else 1e-1
      st_iv  <- sprintf("f(idx_st_cluster,model='besag',graph='%s',scale.model=TRUE,group=idx_temporal,control.group=list(model='rw1',hyper=list(prec=list(prior='pc.prec',param=c(%s,0.01)))),hyper=list(prec=list(prior='pc.prec',param=c(%s,0.01))),diagonal=%s)", cluster_graph, st_u, st_u, st_diag)

      # LinMob x Type IV ST interaction (LinMob_ST_IV).
      lm_st4 <- c("Y ~ 1", fx, bym2, rw2, rv, st_iv)
      if (flags$contact)      lm_st4 <- c(lm_st4, "contact_intensity_std")
      if (flags$spatial_lag)  lm_st4 <- c(lm_st4, "mob_spatial_lag_std")
      fms$LinMob_ST_IV <- as.formula(paste(lm_st4, collapse = " + "))

      # NLMob x Type IV ST interaction (NLMob_ST_IV; was ST_IV pre May 2026).
      st4_terms <- c("Y ~ 1", fx, bym2, rw1_t, rw1_rho, "contact_intensity_std", st_iv)
      if (flags$spatial_lag) st4_terms <- c(st4_terms, "mob_spatial_lag_std")
      fms$NLMob_ST_IV <- as.formula(paste(st4_terms, collapse = " + "))
    }
    fms$Temp <- as.formula(paste(c("Y ~ 1", fx, rw2, rv), collapse = " + "))
  }

  # Reorder to canonical increasing-complexity order.
  ord <- intersect(names(MODEL_LABELS), names(fms))
  fms[ord]
}

nb_hyper_from <- function(dispersion) {
  if (dispersion < 1.5) c(1, 0.01)
  else if (dispersion < 3) c(0.7, 0.01)
  else if (dispersion < 5) c(0.5, 0.01)
  else if (dispersion < 100) c(0.3, 0.05)
  else c(0.1, 0.1)
}

# ---- pre-fit diagnostics: catch known INLA failure modes early ----
diagnose_prefit <- function(df, cluster_graph = NULL) {
  diag <- list(warnings = character(), severe = character())
  dt <- data.table::as.data.table(df)

  # 1. Degenerate clusters: zero/near-zero variance in Y or log_rho
  has_rho <- "log_rho" %in% names(dt)
  if ("idx_st_cluster" %in% names(dt)) {
    cl_stats <- dt[, .(
      var_Y = var(Y, na.rm = TRUE),
      var_rho = if (has_rho) var(log_rho, na.rm = TRUE) else NA_real_,
      n_obs = .N,
      n_nonzero = sum(Y > 0, na.rm = TRUE),
      pct_zero = mean(Y == 0, na.rm = TRUE)
    ), by = idx_st_cluster]
    degen_y <- cl_stats[!is.finite(var_Y) | var_Y < 1e-8]
    degen_rho <- cl_stats[!is.na(var_rho) & (!is.finite(var_rho) | var_rho < 1e-8)]
    all_zero <- cl_stats[n_nonzero == 0]
    if (nrow(degen_y) > 0)
      diag$severe <- c(diag$severe, sprintf("%d clusters have zero/degenerate Y variance", nrow(degen_y)))
    if (nrow(all_zero) > 0)
      diag$severe <- c(diag$severe, sprintf("%d clusters have ALL-zero Y (singular Hessian risk)", nrow(all_zero)))
    if (nrow(degen_rho) > 0)
      diag$warnings <- c(diag$warnings, sprintf("%d clusters have zero rho variance", nrow(degen_rho)))
    diag$cluster_stats <- cl_stats
  }

  # 2. Exposure extremes
  if ("E_inla" %in% names(dt)) {
    e_q <- quantile(dt$E_inla, c(0.01, 0.99), na.rm = TRUE)
    n_tiny <- sum(dt$E_inla < 0.1, na.rm = TRUE)
    n_huge <- sum(dt$E_inla > 1000, na.rm = TRUE)
    if (n_tiny > 0.05 * nrow(dt))
      diag$warnings <- c(diag$warnings, sprintf("%d obs (%.1f%%) have E_inla < 0.1 (log-likelihood explosion risk)",
                                                 n_tiny, 100 * n_tiny / nrow(dt)))
    if (n_huge > 0)
      diag$warnings <- c(diag$warnings, sprintf("%d obs have E_inla > 1000", n_huge))
    diag$exposure_range <- e_q
  }

  # 3. Data sparsity
  zf <- mean(dt$Y == 0, na.rm = TRUE)
  if (zf > 0.98)
    diag$severe <- c(diag$severe, sprintf("Extreme sparsity: %.1f%% zeros (likelihood instability)", 100 * zf))
  else if (zf > 0.95)
    diag$warnings <- c(diag$warnings, sprintf("High sparsity: %.1f%% zeros", 100 * zf))

  # 4. Covariate extremes (possible predictor explosion)
  for (v in c("log_rho", "log_rho_resid", "contact_intensity_std")) {
    if (v %in% names(dt)) {
      rng <- range(dt[[v]], na.rm = TRUE)
      if (any(!is.finite(rng))) {
        diag$severe <- c(diag$severe, sprintf("%s contains Inf/NaN", v))
      } else if (diff(rng) > 20) {
        diag$warnings <- c(diag$warnings, sprintf("%s range is extreme: [%.1f, %.1f]", v, rng[1], rng[2]))
      }
    }
  }

  # 5. Dimensionality check
  n_grids <- max(dt$idx_spatial, na.rm = TRUE)
  n_days <- max(dt$idx_temporal, na.rm = TRUE)
  latent_dim <- n_grids * n_days
  if (latent_dim > 500000)
    diag$warnings <- c(diag$warnings, sprintf("Very large latent field: %d grids x %d days = %s",
                                               n_grids, n_days, format(latent_dim, big.mark = ",")))

  diag
}

# ---- merge_degenerate_clusters: fix zero-variance ST clusters ----
merge_degenerate_clusters <- function(df, cluster_graph_file) {
  dt <- data.table::as.data.table(df)

  cl_stats <- dt[, .(
    var_Y = var(Y, na.rm = TRUE),
    n_nonzero = sum(Y > 0, na.rm = TRUE)
  ), by = idx_st_cluster]

  # Identify degenerate clusters: zero variance OR all zeros
  bad_cl <- cl_stats[!is.finite(var_Y) | var_Y < 1e-8 | n_nonzero == 0, idx_st_cluster]
  if (length(bad_cl) == 0) return(list(df = as.data.frame(dt), changed = FALSE))

  # Read cluster adjacency to find valid neighbors
  if (file.exists(cluster_graph_file)) {
    g_lines <- readLines(cluster_graph_file)
    n_cl <- as.integer(g_lines[1])
    nb_cl <- vector("list", n_cl)
    for (i in seq_len(n_cl)) {
      parts <- as.integer(strsplit(trimws(g_lines[i + 1]), "\\s+")[[1]])
      # INLA graph format per row: <node_id> <k_i> <nbr_1> ... <nbr_{k_i}>
      # so neighbours start at parts[3]. The previous parts[-1] kept k_i as
      # a fake neighbour (e.g. node-with-5-neighbours had cluster `5` added)
      # which corrupted the merge-target selection below.
      nb_cl[[i]] <- if (length(parts) >= 3) parts[3:length(parts)] else integer(0)
    }
  } else {
    return(list(df = as.data.frame(dt), changed = FALSE))
  }

  good_cl <- setdiff(cl_stats$idx_st_cluster, bad_cl)
  if (length(good_cl) == 0) return(list(df = as.data.frame(dt), changed = FALSE))

  # Build merge_map: bad cluster -> target good cluster (a neighbor when
  # possible, else any good cluster). Used both to relabel df and to remap
  # the original adjacency below (preserving sparse topology).
  merge_map <- setNames(seq_along(nb_cl), seq_along(nb_cl))
  for (bc in bad_cl) {
    if (bc <= length(nb_cl)) {
      neighbors <- nb_cl[[bc]]
      valid_nb <- intersect(neighbors, good_cl)
      target <- if (length(valid_nb) > 0) valid_nb[1] else good_cl[1]
    } else {
      target <- good_cl[1]
    }
    merge_map[as.character(bc)] <- target
    dt[idx_st_cluster == bc, idx_st_cluster := target]
  }

  # Re-number clusters to be contiguous
  old_ids <- sort(unique(dt$idx_st_cluster))
  new_map <- setNames(seq_along(old_ids), old_ids)
  dt[, idx_st_cluster := as.integer(new_map[as.character(idx_st_cluster)])]

  # Rebuild cluster adjacency by remapping the ORIGINAL sparse adjacency
  # through merge_map -> new_map. CRITICAL: never use a complete graph K_n
  # here, because besag on K_n grouped by time is computationally
  # catastrophic (dense precision; ST_III/ST_IV hang for hours).
  n_new <- length(old_ids)
  nb_new <- vector("list", n_new)
  for (i in seq_along(nb_cl)) {
    src_old <- i
    src_after_merge <- as.integer(merge_map[as.character(src_old)])
    if (!(src_after_merge %in% old_ids)) next
    src_new <- as.integer(new_map[as.character(src_after_merge)])
    nbrs_old <- nb_cl[[i]]
    nbrs_old <- nbrs_old[nbrs_old > 0]
    if (length(nbrs_old) == 0) next
    nbrs_after_merge <- as.integer(merge_map[as.character(nbrs_old)])
    nbrs_after_merge <- nbrs_after_merge[nbrs_after_merge %in% old_ids]
    if (length(nbrs_after_merge) == 0) next
    nbrs_new <- as.integer(new_map[as.character(nbrs_after_merge)])
    nbrs_new <- setdiff(unique(nbrs_new), src_new)
    nb_new[[src_new]] <- sort(unique(c(nb_new[[src_new]], nbrs_new)))
  }
  # Connect any isolated cluster to its nearest non-empty neighbor list to
  # keep the graph connected (besag requires connected support).
  iso <- which(sapply(nb_new, length) == 0)
  for (ci in iso) {
    others <- setdiff(seq_len(n_new), ci)
    if (length(others) == 0) next
    target <- others[1]
    nb_new[[ci]] <- as.integer(target)
    nb_new[[target]] <- sort(unique(c(nb_new[[target]], as.integer(ci))))
  }
  # Coerce to integer entries (spdep::nb2INLA expects integer vectors; empty
  # neighbour sets are encoded as 0L per the nb convention).
  nb_new <- lapply(nb_new, function(x) if (length(x) == 0) 0L else as.integer(x))
  class(nb_new) <- "nb"
  attr(nb_new, "region.id") <- as.character(seq_len(n_new))
  spdep::nb2INLA(cluster_graph_file, nb_new)

  list(df = as.data.frame(dt), changed = TRUE, n_clusters = n_new)
}

fit_model <- function(formula, df, E_vec, model_name = NULL,
                      zero_frac = 0, dispersion = 1.1,
                      cpo = FALSE, config = TRUE, cv_mode = FALSE,
                      n_grids = Inf, warmstart_theta = NULL) {
  t0 <- Sys.time()
  is_st <- !is.null(model_name) && grepl("ST", model_name)
  small_domain <- is.finite(n_grids) && n_grids < 50
  # CV-mode: force gaussian latent strategy. CV scoring uses posterior
  # MEAN ± SD only (no full marginals are read back), so simplified.laplace's
  # extra cost buys us nothing in CV. ~30% wall-time saving on non-ST models.
  strat <- if (cv_mode) "gaussian"
           else if (small_domain || zero_frac > 0.90 || is_st) "gaussian"
           else "simplified.laplace"
  diag  <- if (small_domain) 0.5
           else if (is_st) 1e-1 else if (zero_frac > 0.90) 1e-2 else 1e-3
  nb_hyper <- nb_hyper_from(dispersion)
  # Auto-switch to Poisson when overdispersion is negligible. With
  # dispersion < 1.1 the NB precision posterior degenerates to a delta
  # at +Inf (NY wave 2 Base: mean 41,401 sd 25 — boundary spike). The
  # pc.prec prior fights the data, the log-likelihood surface w.r.t.
  # theta_NB is flat near the boundary, and Newton iteration wanders
  # for hours. Poisson is the correct family in this regime AND
  # eliminates one hyperparameter (4 -> 3 for Base/LinMob; 5 -> 4 for
  # NLMob), shrinking the joint optimization. WAIC is still comparable
  # within a wave (which is the only ranking we use). May 11 2026.
  use_poisson <- isTRUE(dispersion < 1.1)
  fam <- if (use_poisson) "poisson" else "nbinomial"
  # May 18 2026 (NY-rescue): forcing int.strategy="grid"+diff.logdens=4 when
  # cpo=TRUE was the dominant cost on 4M-row NY ST_IV folds (origin 1 took
  # >1 h on PARDISO and didn't finish). INLA still computes CPO/PIT/failure
  # arrays under int.strategy="eb"+diff.logdens=2 -- it's an empirical-Bayes
  # plug-in CPO instead of grid-integrated CPO, but the LCPO and the PIT
  # uniformity diagnostic are unchanged to the precision we care about
  # for thesis-level cross-validation. Keep CV-fast settings even when cpo
  # is on. Set env CV_CPO_GRID=1 to restore the slow integrated path.
  cpo_grid <- isTRUE(as.logical(Sys.getenv("CV_CPO_GRID", "FALSE")))
  int_strat   <- if (cv_mode && cpo && cpo_grid) "grid" else "eb"
  diff_logdens <- if (cv_mode && cpo && cpo_grid) 4 else if (cv_mode) 2 else 3
  cv_h_val     <- 0.015
  # Restrict summary quantiles in CV to the 2.5/97.5 endpoints used by scoring.
  cv_quantiles <- c(0.025, 0.975)

  # Per-model thread budget. INLA's "4:1" splits work across 4 outer
  # threads; each replicates a chunk of the working memory. The danger
  # case was ST_IV with cpo=TRUE, where LOO predictive densities for ~5M
  # observations made the per-thread footprint blow past 24 GB physical
  # RAM. With CPO globally disabled (CPO_MODELS = character(0)) the
  # ST_IV peak drops by ~6-8x. However on a 24 GB Mac with other apps
  # running, "4:1" still OOMs ST_IV (May 5 2026: 2x SIGKILL). Honour the
  # session-level inla.setOption(num.threads=...) so callers can drop to
  # "2:1" or "1:1" without editing this file. Default to "4:1".
  primary_threads <- INLA::inla.getOption("num.threads") %||% "4:1"

  # Warmstart hook: if a hyperparameter mode is supplied (e.g. from a
  # previously-fitted nested model such as ST_III -> ST_IV), pass it via
  # control.mode. INLA pads the supplied theta vector with zeros for any
  # extra hyperparameters in the larger model, so a partial warmstart is
  # safe and typically halves the inner-Newton iterations needed to find
  # the joint posterior mode.
  ctrl_mode <- if (!is.null(warmstart_theta) && length(warmstart_theta))
                 list(theta = warmstart_theta, restart = TRUE)
               else list()

  # CV diagnostic capture (May 17 2026 evening). When in cv_mode we
  # surface the R-level error message AND any inla stderr (safe-retry
  # banners, lambda_lim crashes, etc.) so the side-log records why a
  # rung failed instead of swallowing into an opaque NULL.
  last_inla_msg <- character(0)
  run_inla <- function(safe_mode = FALSE, diag_override = NULL, strat_override = NULL) {
    d <- if (!is.null(diag_override)) diag_override
         else if (safe_mode) max(diag, 5e-2) else diag
    s <- if (!is.null(strat_override)) strat_override else strat
    # cluster max-quality overrides (May 13 2026). Reads three options set
    # by the SLURM wrapper via env vars (see pipeline/09_fit_model.Rmd
    # cluster override block). Locally these stay NULL so behaviour is
    # unchanged.
    .fit_int   <- getOption("fit_int_strategy", NULL)
    .fit_h     <- getOption("fit_h_val",        NULL)
    .fit_dld   <- getOption("fit_diff_logdens", NULL)
    if (!is.null(.fit_int)) int_strat   <- .fit_int
    if (!is.null(.fit_dld)) diff_logdens <- .fit_dld
    # Finite-difference step for the Hessian of log p(theta | y) at the
    # EB mode. Smaller h -> more precise hyperpar credible intervals but
    # ~k^2 extra full sparse solves on a 4.36 M latent vector. We only
    # consume point summaries (means + quantiles) downstream.
    # May 11 2026 (evening): bumped primary path from h=0.01 -> h=0.005
    # after Poisson Base wave 2 crashed Newton-Raphson 3x with
    # "lambda < 1.0 / lambda_lim" before rescue at h=0.005. With NB
    # gone, the hyperpar likelihood surface is flatter ($\phi_b$ stuck
    # at boundary, only $\tau_b,\tau_c$ identified) and h=0.01 was
    # overshooting. The ~k^2 cost increase is small at k=3-5 hyperpars
    # (Poisson) and the wall-time saving from avoiding rescue retries
    # dominates.
    h_val <- if (cv_mode) cv_h_val else 0.005
    if (!is.null(.fit_h)) h_val <- .fit_h
    # return.marginals.predictor stores a FULL marginal density per linear
    # predictor cell. On NY wave 2 that is 4.36 M densities; combined with
    # CCD integration this drives INLA past 35 GB virtual on a 24 GB box
    # and causes swap thrash (May 10–11 2026: LinMob wave 2 hung 9 h).
    # No downstream pipeline reads marginals.predictor — only
    # summary.fitted.values (means + quantiles) and summary.random are
    # consumed. Force it OFF unconditionally; keep config = caller's
    # choice so inla.posterior.sample is still possible if needed.
    cc <- list(dic = !cv_mode, waic = !cv_mode,
               cpo = cpo, config = config,
               return.marginals.predictor = FALSE)
    if (cv_mode) cc$quantiles <- cv_quantiles
    tryCatch({
      msgs <- capture.output({
        .res <- INLA::inla(formula, family = fam, data = df, E = E_vec,
             control.inla = list(strategy = s, int.strategy = int_strat,
                                 diff.logdens = diff_logdens,
                                 diagonal = d,
                                 h = h_val, cmin = 0,
                                 tolerance.step = 1e-2,
                                 max.iter = getOption("fit_max_iter", 25),
                                 control.vb = list(enable = FALSE)),
             control.compute = cc,
             control.family = if (use_poisson) list() else
               list(hyper = list(theta = list(prior = "pc.prec", param = nb_hyper))),
             control.fixed = list(mean = 0, prec = 0.3),
             control.predictor = list(compute = TRUE, link = 1),
             control.mode = ctrl_mode,
             safe = safe_mode,
             num.threads = primary_threads, verbose = FALSE)
      }, type = "message")
      if (cv_mode && length(msgs) > 0) {
        keep <- grep("(failed|error|crash|lambda|abort|cannot|invalid|NaN|Inf|singular)",
                     msgs, ignore.case = TRUE, value = TRUE)
        if (length(keep)) last_inla_msg <<- c(last_inla_msg,
                                              paste(head(keep, 4), collapse = " | "))
      }
      if (is.null(.res$summary.fitted.values)) {
        if (cv_mode) last_inla_msg <<- c(last_inla_msg, "summary.fitted.values=NULL")
        return(NULL)
      }
      # May 15 2026: accept degenerate-theta fits where WAIC integration
      # fails. NY-w3 Base converged at sigma_spat ~ 0.003 / sigma_temp ~ 1e-8;
      # WAIC became NaN at this boundary but the fit (means, quantiles,
      # DIC, CPO) is otherwise valid. Previously this NULL triggered the
      # full 4-step fallback ladder which dropped CPO and cost ~3 hours.
      .res
    }, error = function(e) {
      if (cv_mode) last_inla_msg <<- c(last_inla_msg,
                                       paste0("R-error: ", conditionMessage(e)))
      NULL
    })
  }
  result <- run_inla(safe_mode = FALSE)
  # May 17 2026 (late): CV-mode warmstart-mismatch rescue. UT_wave1/wave3
  # 10a failed 2/2 origins with "Your model has 5 hyperparameter(s) which
  # is different from the 6 hyperparameter(s) given in 'control.mode'".
  # The warmstart_theta is computed once against the main-fit model, but
  # the CV refit's effective hyperparameter count can differ when the
  # holdout truncation collapses a temporal/spline component. Dropping
  # control.mode lets INLA find its own mode from scratch — slower per
  # origin but unblocks every fold that hits this. Cheap to attempt
  # unconditionally as the first fallback; if warmstart was fine the
  # initial fit already succeeded and this code is not reached.
  if (is.null(result) && length(ctrl_mode)) {
    gc(full = TRUE)
    if (cv_mode) last_inla_msg <- c(last_inla_msg, "<retry: drop warmstart>")
    # FIX 2026-05-17 PM: an empty list() for control.mode let INLA fall
    # back to its internally stored (mismatched) theta, so the
    # "X hyperparameter(s) different from Y" error kept firing across
    # the whole rescue ladder. UT_wave1 10a job 39880573 had 3 origins
    # FAILED with the same dim-mismatch repeated 9x in stderr. Use
    # the documented "start fresh, no warm theta" sentinel instead;
    # warmstart_theta is also nulled so the deeper rungs below cannot
    # reintroduce it.
    ctrl_mode <- list(restart = TRUE)
    warmstart_theta <- NULL
    result <- run_inla(safe_mode = FALSE)
  }
  # May 17 2026: CV-mode rescue ladder. Previously *all* fallbacks were
  # gated by `!cv_mode`, so if the first attempt produced non-finite
  # output (sparse cell-day counts at wave shoulders), the fold was lost
  # with no retry. Run on UT_wave1/UT_wave3 had 15/15 origins FAILED.
  # Allow two cheap safe-mode retries in CV; deeper fallbacks remain
  # full-fit only (they request finer integration grids that don't make
  # sense at CV's coarsened settings).
  if (is.null(result)) {
    gc(full = TRUE)
    result <- run_inla(safe_mode = TRUE, diag_override = if (small_domain) 0.8 else max(diag, 5e-2))
  }
  if (is.null(result)) {
    gc(full = TRUE)
    result <- run_inla(safe_mode = TRUE, diag_override = if (small_domain) 1.0 else 0.1, strat_override = "gaussian")
  }
  # CV-mode Poisson fallback (May 17 2026 evening). When NB safe-retry
  # ladder exhausts and we're in cv_mode, try a Poisson refit. This
  # eliminates the size hyperparameter that is the most common
  # convergence-failure source under temporal masking (dispersion at
  # boundary).
  if (is.null(result) && cv_mode) {
    gc(full = TRUE)
    .saved_fam <- fam
    .saved_use_p <- use_poisson
    fam <- "poisson"
    use_poisson <- TRUE
    result <- run_inla(safe_mode = TRUE,
                       diag_override = if (small_domain) 1.0 else 0.1,
                       strat_override = "gaussian")
    if (!is.null(result)) {
      last_inla_msg <- c(last_inla_msg, "<rescued by Poisson fallback>")
    } else {
      fam <- .saved_fam
      use_poisson <- .saved_use_p
    }
  }
  if (cv_mode && is.null(result) && length(last_inla_msg)) {
    message("[fit_model CV-mode FAIL] ", paste(last_inla_msg, collapse = " ;; "))
  }
  if (is.null(result) && !cv_mode) {
    gc(full = TRUE)
    r3_diag <- if (small_domain) 1.0 else 0.1
    tryCatch({
      invisible(capture.output({
        .res <- INLA::inla(formula, family = fam, data = df, E = E_vec,
             control.inla = list(strategy = "gaussian", int.strategy = "eb",
                                 diagonal = r3_diag, h = 0.005, cmin = 0,
                                 tolerance.step = 1e-2, max.iter = 25,
                                 control.vb = list(enable = FALSE)),
             # May 15 2026: honour caller's cpo flag in deep fallbacks.
             # Previously we dropped CPO here which made LCPO unavailable
             # for any pair that needed the 4th/5th retry.
             control.compute = list(dic = TRUE, waic = TRUE, cpo = cpo, config = config),
             control.family = if (use_poisson) list() else
               list(hyper = list(theta = list(prior = "pc.prec", param = nb_hyper))),
             control.fixed = list(mean = 0, prec = 0.3),
             control.predictor = list(compute = TRUE, link = 1),
             safe = TRUE,
             num.threads = "2:1", verbose = FALSE)
      }, type = "message"))
      if (!is.null(.res$summary.fitted.values)) result <- .res
    }, error = function(e) NULL)
  }
  if (is.null(result) && !cv_mode) {
    gc(full = TRUE)
    r4_diag <- if (small_domain) 1.0 else 0.2
    r4_strat <- if (small_domain) "gaussian" else "simplified.laplace"
    tryCatch({
      invisible(capture.output({
        .res <- INLA::inla(formula, family = fam, data = df, E = E_vec,
             control.inla = list(strategy = r4_strat, int.strategy = "eb",
                                 diagonal = r4_diag, h = 0.005, cmin = 0,
                                 tolerance.step = 1e-2, max.iter = 25,
                                 control.vb = list(enable = FALSE)),
             control.compute = list(dic = TRUE, waic = TRUE, cpo = cpo, config = config),
             control.family = if (use_poisson) list() else
               list(hyper = list(theta = list(prior = "pc.prec", param = nb_hyper))),
             control.fixed = list(mean = 0, prec = 0.3),
             control.predictor = list(compute = TRUE, link = 1),
             control.mode = list(restart = TRUE),
             safe = TRUE,
             num.threads = "1:1", verbose = FALSE)
      }, type = "message"))
      if (!is.null(.res$summary.fitted.values)) result <- .res
    }, error = function(e) NULL)
  }

  runtime <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
  if (is.null(result)) return(list(result = NULL, runtime = runtime, status = "FAILED",
                                   convergence = NULL))

  if (cv_mode) return(list(result = result, runtime = runtime))

  waic_val <- result$waic$waic
  pD <- result$dic$p.eff
  max_fit <- max(result$summary.fitted.values$mean, na.rm = TRUE)
  # May 15 2026: NA WAIC at degenerate-theta boundaries is no longer fatal.
  # The fit is OK if fitted means are finite and DIC's p.eff is sensible.
  # WAIC NA gets propagated downstream; comparison_table uses DIC as a
  # backup ordering when WAIC is missing.
  waic_bad <- !is.null(waic_val) && !is.na(waic_val) &&
              (is.infinite(waic_val) || waic_val > 1e7)
  pD_bad   <- !is.null(pD) && !is.na(pD) && pD < 0
  status <- if (waic_bad || pD_bad) "UNSTABLE"
            else if (max_fit > 1e6) "IMPLAUSIBLE" else "OK"
  # Salvage CPO when caller asked for it and main path returned no CPO
  # (e.g. degenerate-theta fits that converged in gaussian mode).
  if (status == "OK" && cpo && is.null(result$cpo)) {
    result <- tryCatch(INLA::inla.cpo(result, verbose = FALSE),
                       error = function(e) result)
  }
  # Recompute failed CPO values for reliable LCPO (cap to avoid hours-long LOO refitting)
  if (status == "OK" && cpo && !is.null(result$cpo$failure)) {
    n_fail <- sum(result$cpo$failure > 0, na.rm = TRUE)
    if (n_fail > 0 && n_fail <= 5000) {
      result <- tryCatch(INLA::inla.cpo(result, verbose = FALSE),
                         error = function(e) result)
    }
  }

  list(result = if (status == "OK") result else NULL, runtime = runtime, status = status,
       convergence = list(
         mlik_int    = result$mlik["log marginal-likelihood (integration)", 1],
         mlik_gauss  = result$mlik["log marginal-likelihood (Gaussian)", 1],
         mode_status = result$mode$mode.status,
         n_hyperpar  = nrow(result$summary.hyperpar),
         cpu_total_s = sum(result$cpu.used, na.rm = TRUE),
         max_fitted  = max_fit, waic = waic_val, pD = pD
       ))
}

# ---- extract_model_outputs: pull RR, fixed/spatial/temporal/fitted from INLA result ----
extract_model_outputs <- function(result, df, grid_index, model_name, E_vec = df$E_inla) {
  if (is.null(result)) return(NULL)
  out <- list()

  if (!is.null(result$summary.fixed)) {
    fx <- data.table::as.data.table(result$summary.fixed, keep.rownames = "covariate")
    rename_quants(fx)
    fx[, `:=`(RR = exp(mean), RR_lower = exp(q025), RR_upper = exp(q975),
              significant = (q025 > 0) | (q975 < 0))]
    out$fixed_effects <- fx
  }

  if ("idx_spatial" %in% names(result$summary.random)) {
    sp <- data.table::as.data.table(result$summary.random$idx_spatial)
    names(sp)[1] <- "idx_spatial"; rename_quants(sp)
    sp <- merge(sp, grid_index, by = "idx_spatial", all.x = TRUE)
    sp[, RR_spatial := exp(mean)]
    out$spatial_effects <- sp
  }

  if ("idx_temporal" %in% names(result$summary.random)) {
    tp <- data.table::as.data.table(result$summary.random$idx_temporal)
    names(tp)[1] <- "idx_temporal"; rename_quants(tp)
    tp[, RR_temporal := exp(mean)]
    out$temporal_effects <- tp[order(idx_temporal)]
  }

  if (!is.null(result$summary.fitted.values)) {
    ft <- data.table::as.data.table(result$summary.fitted.values)
    data.table::setnames(ft, names(ft), paste0("fitted_", names(ft)))
    ft[, `:=`(Y = df$Y, E = E_vec, grid_id = df$grid_id, time_id = df$time_id,
              pearson_resid = (df$Y - fitted_mean) / sqrt(pmax(fitted_mean, 1e-6)))]
    if (!is.null(result$summary.linear.predictor)) {
      lp <- data.table::as.data.table(result$summary.linear.predictor)
      ft[, `:=`(RR = exp(lp$mean), RR_lower = exp(lp[["0.025quant"]]),
                RR_upper = exp(lp[["0.975quant"]]),
                RR_exceed1   = 1 - pnorm(0,          mean = lp$mean, sd = lp$sd),
                RR_exceed1.5 = 1 - pnorm(log(1.5),   mean = lp$mean, sd = lp$sd),
                RR_exceed2   = 1 - pnorm(log(2),     mean = lp$mean, sd = lp$sd),
                CI_width = exp(lp[["0.975quant"]]) - exp(lp[["0.025quant"]]))]
    }
    out$fitted_values <- ft
  }

  if (!is.null(result$summary.hyperpar))
    out$hyperparameters <- data.table::as.data.table(result$summary.hyperpar, keep.rownames = "parameter")
  if ("log_rho_idx" %in% names(result$summary.random))
    out$rw1_curve <- data.table::as.data.table(result$summary.random$log_rho_idx)

  # Spatiotemporal interaction effects (Knorr-Held ST_I–IV models)
  if ("idx_st_cluster" %in% names(result$summary.random)) {
    st <- data.table::as.data.table(result$summary.random$idx_st_cluster)
    names(st)[1] <- "idx_st"; rename_quants(st)
    st[, RR_st := exp(mean)]
    out$st_interaction <- st
  }

  out
}

# ---- select_covariates: validate, filter, and VIF-screen covariates ----
select_covariates <- function(df, potential, vif_threshold = 5) {
  avail <- c()
  for (v in potential) {
    if (!v %in% names(df)) next
    ok <- if (v %in% c("landuse_commercial")) sum(df[[v]], na.rm = TRUE) > 10
          else if (v == "is_weekend" || grepl("_x_", v)) TRUE
          else sd(df[[v]], na.rm = TRUE) > 0.01
    if (ok) avail <- c(avail, v)
  }
  # -- VIF screening: drop covariates with high collinearity --
  numeric_avail <- avail[sapply(avail, function(v) is.numeric(df[[v]]) && !v %in% c("is_weekend"))]
  if (length(numeric_avail) >= 2) {
    X <- as.data.frame(df[, numeric_avail, drop = FALSE])
    X <- X[complete.cases(X), , drop = FALSE]
    if (nrow(X) >= 2 * length(numeric_avail)) {
      repeat {
        vifs <- sapply(seq_along(numeric_avail), function(j) {
          r2 <- summary(lm(X[[j]] ~ ., data = X[, -j, drop = FALSE]))$r.squared
          1 / (1 - r2)
        })
        names(vifs) <- numeric_avail
        worst <- which.max(vifs)
        if (vifs[worst] > vif_threshold && length(numeric_avail) > 1) {
          avail <- setdiff(avail, numeric_avail[worst])
          numeric_avail <- numeric_avail[-worst]
          X <- X[, -worst, drop = FALSE]
        } else break
      }
    }
  }
  avail
}

# ── Temporal effect despiking ─────────────────────────────────────
# Remove month-boundary composition spikes from temporal effects.
# Detects single-day jumps > threshold at known month boundaries
# and replaces with linear interpolation from neighbors.
despike_temporal <- function(tp, month_breaks, col = "mean", threshold = 0.15) {
  x <- tp[[col]]; n <- length(x)
  for (mb in month_breaks) {
    idx <- which.min(abs(as.numeric(tp$date - mb)))
    if (idx < 3 || idx > n - 2) next
    d_in  <- x[idx] - x[idx - 1]
    d_out <- x[idx + 1] - x[idx]
    if (abs(d_in) > threshold && sign(d_in) != sign(d_out) &&
        abs(d_out) > threshold * 0.5)
      x[idx] <- (x[idx - 1] + x[idx + 1]) / 2
  }
  x
}
