## generate_cpo_pit_plots.R
## Cross-wave CPO/PIT diagnostics (LCPO summary + PIT histograms).
## Inputs: per-wave `validation/temporal_cv_results.rds` (fold_diag, holdout_long).
## Output: `outputs/spatiotemporal_risk/_cross_wave/cpo_pit/` (PDFs + CSV).

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

waves <- c("NY_wave1", "NY_wave2", "NY_wave3",
           "UT_wave1", "UT_wave2", "UT_wave3")
base  <- "outputs/spatiotemporal_risk"
out   <- file.path(base, "_cross_wave", "cpo_pit")
dir.create(out, showWarnings = FALSE, recursive = TRUE)

set.seed(20260521L)  # reproducible randomised PIT

## ---- randomised PIT for count predictives ---------------------------------
## y ~ Poisson(mu) or NegBin(mu, size). Randomised PIT smooths the
## discreteness so a U(0,1) histogram is the calibration target.
rand_pit <- function(y, mu, family, size = NA_real_) {
  y  <- pmax(0L, as.integer(round(y)))
  mu <- pmax(1e-9, mu)
  if (family == "poisson") {
    Fy   <- ppois(y, mu)
    Fym1 <- ppois(y - 1L, mu)
  } else {  # nbinomial
    sz <- ifelse(is.finite(size) & size > 0, size, 1e8)
    Fy   <- pnbinom(y, size = sz, mu = mu)
    Fym1 <- pnbinom(y - 1L, size = sz, mu = mu)
  }
  u <- runif(length(y))
  Fym1 + u * (Fy - Fym1)
}

## ---- gather per-wave numerics + cell PIT ----------------------------------
pit_long    <- vector("list", length(waves))
pit_w1_long <- vector("list", length(waves))
wave_stat   <- vector("list", length(waves))

for (w in waves) {
  rds_fp <- file.path(base, w, "validation", "temporal_cv_results.rds")
  if (!file.exists(rds_fp)) { warning("missing rds: ", w); next }
  r  <- readRDS(rds_fp)
  fd <- as.data.table(r$fold_diag)
  hl <- as.data.table(r$holdout_long)
  ws <- as.data.table(r$window_scores)

  ## attach family / nb_size per origin
  fam_tab <- fd[, .(origin, family, nb_size)]
  hl <- merge(hl, fam_tab, by = "origin", all.x = TRUE)
  hl <- hl[is.finite(mu) & is.finite(y_true)]

  ## randomised PIT (vectorise per family group for speed)
  hl <- hl[!is.na(family)]
  hl[, pit := NA_real_]
  for (fam in unique(hl$family)) {
    idx <- which(hl$family == fam)
    if (!length(idx)) next
    sz  <- if (fam == "nbinomial") hl$nb_size[idx] else NA_real_
    hl$pit[idx] <- rand_pit(hl$y_true[idx], hl$mu[idx], fam, sz)
  }

  pit_long[[w]] <- data.table(wave = w, origin = hl$origin, pit = hl$pit)

  ## per-origin Wasserstein-1 to U(0,1), computed from the per-cell PIT we just made
  ## (avoids dependency on possibly-stale ws$pit_W1)
  w1_per_origin <- hl[!is.na(pit), .(
    pit_W1 = mean(abs(sort(pit) - (seq_len(.N) - 0.5) / .N))
  ), by = origin]

  pit_w1_long[[w]] <- w1_per_origin[, .(wave = w, origin, pit_W1)]

  wave_stat[[w]] <- data.table(
    wave             = w,
    family           = paste(sort(unique(fd$family)), collapse = "/"),
    n_origins        = nrow(fd),
    mean_LCPO        = mean(fd$LCPO,     na.rm = TRUE),
    sd_LCPO          = sd(fd$LCPO,       na.rm = TRUE),
    cpo_failures_tot = sum(fd$cpo_failures, na.rm = TRUE),
    mean_pit_W1      = mean(w1_per_origin$pit_W1, na.rm = TRUE),
    n_cells_pit      = nrow(hl)
  )
}

pit_dt   <- rbindlist(pit_long,  fill = TRUE)
stat_dt  <- rbindlist(wave_stat, fill = TRUE)
fwrite(stat_dt, file.path(out, "cpo_pit_per_wave.csv"))
print(stat_dt)

## ---- per-origin LCPO long table (for the bar chart) -----------------------
origins_long <- rbindlist(lapply(waves, function(w) {
  rds_fp <- file.path(base, w, "validation", "temporal_cv_results.rds")
  if (!file.exists(rds_fp)) return(NULL)
  fd <- as.data.table(readRDS(rds_fp)$fold_diag)
  fd[, .(wave = w, origin, LCPO, family)]
}), fill = TRUE)

origins_long[, wave := factor(wave, levels = waves)]
stat_dt[, wave    := factor(wave, levels = waves)]
pit_dt [, wave    := factor(wave, levels = waves)]

## ---- plot 1: PIT histograms (6-panel) -------------------------------------
p_pit <- ggplot(pit_dt, aes(x = pit)) +
  geom_histogram(breaks = seq(0, 1, by = 0.05),
                 fill = "#4C72B0", colour = "white", linewidth = 0.2) +
  geom_hline(data = pit_dt[, .(yref = .N / 20), by = wave],
             aes(yintercept = yref),
             linetype = "dashed", colour = "grey30", linewidth = 0.4) +
  facet_wrap(~ wave, ncol = 3, scales = "free_y") +
  labs(x = "Randomised PIT", y = "Count",
       title = "PIT histograms by wave (pooled across origins, h = 3)",
       subtitle = "Dashed line = uniform reference") +
  theme_bw(base_size = 10) +
  theme(strip.background = element_rect(fill = "grey92", colour = NA),
        panel.grid.minor = element_blank())

ggsave(file.path(out, "pit_histograms_by_wave.pdf"),
       p_pit, width = 9, height = 5.5)
ggsave("thesis_latex/figures/cv_pit_histograms_oof.png",
       p_pit, width = 9, height = 5.5, dpi = 300, bg = "white")

## ---- plot 2: LCPO bar + per-origin dots ----------------------------------
p_lcpo <- ggplot(stat_dt, aes(x = wave, y = mean_LCPO, fill = family)) +
  geom_col(width = 0.65, alpha = 0.85) +
  geom_errorbar(aes(ymin = mean_LCPO - sd_LCPO,
                    ymax = mean_LCPO + sd_LCPO),
                width = 0.18, colour = "grey25") +
  geom_jitter(data = origins_long, aes(x = wave, y = LCPO),
              inherit.aes = FALSE, width = 0.12, height = 0,
              shape = 21, fill = "white", colour = "grey15",
              size = 1.6, stroke = 0.3) +
  scale_fill_manual(values = c(poisson = "#4C72B0", nbinomial = "#C44E52")) +
  labs(x = NULL, y = "Mean LCPO (per origin)",
       fill = "INLA family",
       title = "Cross-wave log-CPO (h = 3 day forecast)",
       subtitle = "Bars: mean across origins ± 1 sd; dots: individual origins") +
  theme_bw(base_size = 10) +
  theme(panel.grid.minor = element_blank(),
        legend.position = "bottom")

ggsave(file.path(out, "lcpo_by_wave.pdf"),
       p_lcpo, width = 8, height = 4.5)
ggsave("thesis_latex/figures/cv_lcpo_by_wave.png",
       p_lcpo, width = 8, height = 4.5, dpi = 300, bg = "white")

## ---- plot 3: pit_W1 per origin -------------------------------------------
ws_long <- rbindlist(pit_w1_long, fill = TRUE)
ws_long[, wave := factor(wave, levels = waves)]

p_w1 <- ggplot(ws_long, aes(x = factor(origin), y = pit_W1)) +
  geom_col(fill = "#55A868", alpha = 0.85, width = 0.7) +
  facet_wrap(~ wave, ncol = 3, scales = "free_x") +
  labs(x = "Origin (training endpoint, day index)",
       y = expression(W[1]*"(PIT, U(0,1))"),
       title = "PIT Wasserstein-1 distance to uniform, per origin",
       subtitle = "Lower = better calibration; 0 = exact uniform PIT") +
  theme_bw(base_size = 10) +
  theme(strip.background = element_rect(fill = "grey92", colour = NA),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 0, size = 7))

ggsave(file.path(out, "pit_w1_by_wave_origin.pdf"),
       p_w1, width = 9, height = 5)
ggsave("thesis_latex/figures/cv_pit_w1_by_origin.png",
       p_w1, width = 9, height = 5, dpi = 300, bg = "white")

message("Wrote plots to ", out)
