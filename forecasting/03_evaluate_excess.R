# 03_evaluate_excess.R
# Mexico measles — World Cup 2026 forecast benchmark
#
# Reads the bootstrap predictive curves saved by run_forecasts.m and computes,
# for every analysis x state x model (GLM, GRM, RICH, ENS):
#   tables/weekly_forecast_eval.csv   per forecast week: observed, median,
#                                     50/95% PI, exceedance, ratio, excess, WIS,
#                                     predictive percentile of the observation
#   tables/cumulative_excess.csv      primary (15 Jun-20 Jul) and full
#                                     (15 Jun-3 Aug) windows: observed vs
#                                     cumulative forecast distribution
#   tables/forecast_summary.csv       window-level MAE, WIS, 95% coverage,
#                                     number of weeks above the 95% PI
#
# All quantities come from the same predictive draws (observation noise
# included), so models and the ensemble are treated identically.

suppressPackageStartupMessages({ library(dplyr); library(tidyr) })

PROJ <- if (interactive() && requireNamespace("rstudioapi", quietly = TRUE) &&
            rstudioapi::isAvailable()) {
  dirname(rstudioapi::getSourceEditorContext()$path)
} else {
  args <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", args[grepl("^--file=", args)])
  if (length(f)) dirname(normalizePath(f)) else getwd()
}
setwd(PROJ)

MODELS <- c("GLM", "GRM", "RICH", "ENS")
WIS_ALPHAS <- c(0.02, 0.05, seq(0.1, 0.9, by = 0.1))   # K = 11, as in the paper

# ── Weighted interval score from predictive draws ────────────────────────────
wis_from_draws <- function(draws, y, alphas = WIS_ALPHAS) {
  m  <- median(draws)
  is <- vapply(alphas, function(a) {
    l <- quantile(draws, a / 2, names = FALSE, type = 7)
    u <- quantile(draws, 1 - a / 2, names = FALSE, type = 7)
    (u - l) + (2 / a) * (l - y) * (y < l) + (2 / a) * (y - u) * (y > u)
  }, numeric(1))
  (0.5 * abs(y - m) + sum((alphas / 2) * is)) / (length(alphas) + 0.5)
}

q <- function(x, p) quantile(x, p, names = FALSE, type = 7)

manifest <- read.csv("manifest.csv", stringsAsFactors = FALSE)
groups   <- distinct(manifest, analysis, state, input_file)

weekly <- list(); cumul <- list()

for (g in seq_len(nrow(groups))) {
  an <- groups$analysis[g]; st <- groups$state[g]
  dates <- read.csv(sub("\\.txt$", "_dates.csv", groups$input_file[g]),
                    stringsAsFactors = FALSE) %>%
    mutate(week_start = as.Date(week_start))

  for (mdl in MODELS) {
    run_dir <- file.path("runs", an, st, mdl)
    f <- file.path(run_dir, "curves_predictive.csv")
    if (!file.exists(file.path(run_dir, "DONE")) || !file.exists(f)) {
      message(sprintf("  %s/%s/%s: not run yet, skipped", an, st, mdl)); next
    }
    draws <- as.matrix(read.csv(f, header = FALSE))
    if (nrow(draws) != nrow(dates)) {
      stop(sprintf("%s: %d curve rows but %d weeks in dates file",
                   run_dir, nrow(draws), nrow(dates)))
    }

    fc  <- which(dates$role != "calibration")
    # evaluation weeks are the same calendar weeks in every analysis
    ev  <- which(dates$role %in% c("forecast_primary", "forecast_extended"))
    for (i in fc) {
      d <- draws[i, ]; y <- dates$cases[i]
      weekly[[length(weekly) + 1]] <- tibble(
        analysis = an, state = st, model = mdl,
        week_start = dates$week_start[i], horizon = i - min(fc) + 1,  # weeks after calibration end
        window = dates$role[i],
        observed = y, nowcast_mature = dates$nowcast_mature[i],
        median = q(d, .5), lower50 = q(d, .25), upper50 = q(d, .75),
        lower95 = q(d, .025), upper95 = q(d, .975),
        above95 = y > q(d, .975), below95 = y < q(d, .025),
        ratio = ifelse(q(d, .5) > 0, y / q(d, .5), NA_real_),
        excess = y - q(d, .5),
        pred_percentile = mean(d < y) + 0.5 * mean(d == y),
        wis = wis_from_draws(d, y))
    }

    for (w in c("primary", "full")) {
      rows <- if (w == "primary") which(dates$role == "forecast_primary") else ev
      cd <- colSums(draws[rows, , drop = FALSE]); y <- sum(dates$cases[rows])
      cumul[[length(cumul) + 1]] <- tibble(
        analysis = an, state = st, model = mdl, window = w,
        first_week = min(dates$week_start[rows]),
        last_week  = max(dates$week_start[rows]),
        n_weeks = length(rows), observed = y,
        median = q(cd, .5), lower95 = q(cd, .025), upper95 = q(cd, .975),
        excess = y - q(cd, .5),
        excess_lower95 = y - q(cd, .975), excess_upper95 = y - q(cd, .025),
        ratio = ifelse(q(cd, .5) > 0, y / q(cd, .5), NA_real_),
        # one-sided upper-tail probability: P(forecast >= observed)
        p_upper = mean(cd >= y),
        above95 = y > q(cd, .975))
    }
  }
}

if (!length(weekly)) stop("No completed runs found under runs/ - run run_forecasts.m first")

weekly <- bind_rows(weekly); cumul <- bind_rows(cumul)

summary_tbl <- bind_rows(
  weekly %>% filter(window == "forecast_primary") %>% mutate(window = "primary"),
  weekly %>% filter(window != "forecast_lead") %>% mutate(window = "full")) %>%
  group_by(analysis, state, model, window) %>%
  summarise(n_weeks = n(),
            weeks_above95 = sum(above95), weeks_below95 = sum(below95),
            coverage95 = mean(!above95 & !below95),
            mae = mean(abs(observed - median)),
            mean_wis = mean(wis), .groups = "drop")

dir.create("tables", showWarnings = FALSE)
write.csv(weekly,      "tables/weekly_forecast_eval.csv", row.names = FALSE)
write.csv(cumul,       "tables/cumulative_excess.csv",    row.names = FALSE)
write.csv(summary_tbl, "tables/forecast_summary.csv",     row.names = FALSE)

cat("\nCumulative excess (primary window):\n")
print(as.data.frame(cumul %>% filter(window == "primary") %>%
  transmute(analysis, state, model, observed, median = round(median, 1),
            PI95 = sprintf("%.0f-%.0f", lower95, upper95),
            excess = round(excess, 1), p_upper = round(p_upper, 3))),
  row.names = FALSE)
cat("\nWritten: tables/weekly_forecast_eval.csv, cumulative_excess.csv, forecast_summary.csv\n")
