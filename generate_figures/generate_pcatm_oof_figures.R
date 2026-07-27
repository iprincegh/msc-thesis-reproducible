#!/usr/bin/env Rscript
# generate_pcatm_oof_figures.R
# Out-of-fold PCAtM summary and per-fold diagnostics (winsorised variants).
#+ Reads per-wave `validation/pcatm_oof_winsorised.csv` and writes PNGs.

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(tidyr)
  library(readr)
  library(stringr)
})

# Locate project root robustly: resolve script path from --file= arg if present,
# else fall back to cwd.
.script_path <- (function() {
  args <- commandArgs(trailingOnly = FALSE)
  m <- regmatches(args, regexpr("(?<=^--file=).+", args, perl = TRUE))
  if (length(m)) normalizePath(m[1], mustWork = FALSE) else NA_character_
})()
PROJECT_ROOT <- if (!is.na(.script_path)) {
  normalizePath(file.path(dirname(.script_path), ".."),
                winslash = "/", mustWork = TRUE)
} else getwd()

ST_DIR       <- file.path(PROJECT_ROOT, "outputs", "spatiotemporal_risk")
CROSS_DIR    <- file.path(ST_DIR, "_cross_wave")
FIG_DIR      <- file.path(PROJECT_ROOT, "thesis_latex", "figures")
SUMMARY_CSV  <- file.path(CROSS_DIR, "pcatm_oof_winsorised_summary.csv")

stopifnot(file.exists(SUMMARY_CSV))
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

wave_levels <- c("NY_wave1", "NY_wave2", "NY_wave3",
                 "UT_wave1", "UT_wave2", "UT_wave3")
wave_labels <- c("NYC W1", "NYC W2", "NYC W3",
                 "Utah W1", "Utah W2", "Utah W3")

state_palette <- c(NYC = "#2C7FB8", Utah = "#D95F0E")

# ---------------------------------------------------------------------------
# 1. Cross-wave summary figure
# ---------------------------------------------------------------------------
summ <- read_csv(SUMMARY_CSV, show_col_types = FALSE) |>
  mutate(
    wave   = factor(wave, levels = wave_levels, labels = wave_labels),
    state  = ifelse(str_starts(as.character(wave), "NYC"), "NYC", "Utah"),
    point  = ifelse(prefer_clip, pcatm_oof_clip2_robust, pcatm_oof_robust),
    lo     = pcatm_draws_clip_lo,
    hi     = pcatm_draws_clip_hi,
    label  = sprintf("%.1f%%", 100 * point),
    degen_note = ifelse(n_folds_degenerate > 0,
                        sprintf("%d/%d folds dropped",
                                n_folds_degenerate, n_folds_total),
                        sprintf("%d folds", n_folds_robust))
  )

p_summary <- ggplot(summ, aes(x = wave, y = point, fill = state)) +
  geom_hline(yintercept = 0, colour = "grey60", linewidth = 0.3) +
  geom_col(width = 0.65, alpha = 0.85) +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.18, linewidth = 0.45) +
  geom_text(aes(label = label,
                y = ifelse(point >= 0,
                           pmax(hi, point) + 0.04,
                           pmin(lo, point) - 0.04)),
            size = 3.1, fontface = "bold") +
  geom_text(aes(label = degen_note, y = -0.08),
            size = 2.4, colour = "grey35", vjust = 1) +
  scale_fill_manual(values = state_palette, guide = "none") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                     breaks = scales::pretty_breaks(6)) +
  labs(x = NULL, y = "PCATM (OOF, winsorised)",
       title = "Out-of-fold PCATM by wave",
       subtitle = "Bars: robust point estimate (clip-preferred when extreme \u03b7 present); whiskers: 95% credible interval over posterior draws") +
  theme_minimal(base_size = 11) +
  theme(plot.title    = element_text(face = "bold"),
        plot.subtitle = element_text(colour = "grey30", size = 9),
        panel.grid.major.x = element_blank())

ggsave(file.path(FIG_DIR, "pcatm_oof_winsorised_summary.png"),
       p_summary, width = 9, height = 4.6, dpi = 300, bg = "white")
message("wrote: pcatm_oof_winsorised_summary.png")

# ---------------------------------------------------------------------------
# 2. Per-fold dotplot
# ---------------------------------------------------------------------------
per_fold <- lapply(wave_levels, function(w) {
  f <- file.path(ST_DIR, w, "validation", "pcatm_oof_winsorised.csv")
  if (!file.exists(f)) return(NULL)
  read_csv(f, show_col_types = FALSE) |>
    mutate(wave = w)
}) |>
  bind_rows() |>
  mutate(
    wave_lab = factor(wave, levels = wave_levels, labels = wave_labels),
    state    = ifelse(str_starts(wave, "NY"), "NYC", "Utah"),
    point    = pcatm_oof_clip2,
    is_degenerate = as.logical(is_degenerate)
  )

y_lo <- -1.0; y_hi <- 1.2
kept_only <- per_fold |> filter(!is_degenerate)
drop_tab  <- per_fold |> group_by(wave_lab) |>
             summarise(n_drop = sum(is_degenerate), .groups = "drop") |>
             filter(n_drop > 0)

p_fold <- ggplot(kept_only,
                 aes(x = wave_lab, y = point, colour = state)) +
  geom_hline(yintercept = 0, colour = "grey60", linewidth = 0.3) +
  geom_jitter(width = 0.18, height = 0, size = 2.0, alpha = 0.8, shape = 16) +
  stat_summary(fun = median, geom = "crossbar",
               width = 0.45, linewidth = 0.35, colour = "black", fatten = 1.4) +
  geom_text(data = drop_tab,
            aes(x = wave_lab, y = y_lo + 0.05,
                label = sprintf("%d degenerate folds dropped", n_drop)),
            inherit.aes = FALSE, size = 2.8, colour = "grey25", fontface = "italic") +
  scale_colour_manual(values = state_palette, guide = "none") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                     breaks = scales::pretty_breaks(6),
                     limits = c(y_lo, y_hi), oob = scales::squish) +
  labs(x = NULL, y = "Per-fold PCATM (clip\u03b7\u00b12)",
       title = "Per-fold OOF PCATM variability (kept folds)",
       subtitle = "Black bar: median across kept folds. Degenerate folds (>50% mob-cov columns near-constant) excluded.") +
  theme_minimal(base_size = 11) +
  theme(plot.title    = element_text(face = "bold"),
        plot.subtitle = element_text(colour = "grey30", size = 9),
        panel.grid.major.x = element_blank())

ggsave(file.path(FIG_DIR, "pcatm_oof_per_fold.png"),
       p_fold, width = 9, height = 4.8, dpi = 300, bg = "white")
message("wrote: pcatm_oof_per_fold.png")

invisible(NULL)
