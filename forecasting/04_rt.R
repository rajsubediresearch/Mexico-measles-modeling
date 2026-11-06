# 04_rt.R
# Mexico measles — effective reproduction number
#
# Two estimates, using the same generation interval
# (gamma, mean 11.7 days, SD 2 days, discretised to weeks):
#
#   1. Model-based Rt: renewal equation applied to each bootstrap trajectory
#      of the fitted/forecast incidence (curves_process.csv), as in v8:
#        R_t = I_t / sum_{j>=1} I_{t-j} g_j
#      Values in forecast weeks are projections of the pre-tournament fit.
#      First two weeks of each series are not reported (unstable).
#
#   2. Data-based Rt: Cori et al. (2013) estimator applied to the observed
#      mature weekly counts, over a sliding window of TAU weeks, so that
#      tournament-period weeks are informed by tournament-period data.
#
# Outputs: tables/rt_model.csv, tables/rt_data.csv, tables/rt_summary.csv

suppressPackageStartupMessages({ library(readxl); library(dplyr); library(tidyr) })

PROJ <- if (interactive() && requireNamespace("rstudioapi", quietly = TRUE) &&
            rstudioapi::isAvailable()) {
  dirname(rstudioapi::getSourceEditorContext()$path)
} else {
  args <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", args[grepl("^--file=", args)])
  if (length(f)) dirname(normalizePath(f)) else getwd()
}
setwd(PROJ)

# ── Settings ─────────────────────────────────────────────────────────────────
GI_MEAN_DAYS <- 11.7
GI_SD_DAYS   <- 2
TAU          <- 2          # Cori sliding window (weeks)
PRIOR_SHAPE  <- 1          # Cori gamma prior on R: mean 5, SD 5
PRIOR_SCALE  <- 5
MATURE_RELEASE <- "2026-09-21"
LAST_WEEK    <- as.Date("2026-08-03")
STATES <- c(CDMX = "CIUDAD DE MEXICO", JAL = "JALISCO")
MODELS <- c("GLM", "GRM", "RICH", "ENS")

# ── Weekly generation-interval distribution ──────────────────────────────────
gi_shape <- (GI_MEAN_DAYS / GI_SD_DAYS)^2
gi_scale <- GI_SD_DAYS^2 / GI_MEAN_DAYS
g <- diff(pgamma(7 * (0:4), shape = gi_shape, scale = gi_scale))
g <- g / sum(g)                       # g[j] = P(generation interval in week j)
cat("Weekly generation interval g_1..g_4:", round(g, 4), "\n")

# total infectiousness: Lambda_t = sum_j I_{t-j} g_j   (matrix-friendly)
lambda_mat <- function(I) {           # I: T x N
  I <- as.matrix(I); Tn <- nrow(I)
  L <- matrix(0, Tn, ncol(I))
  for (j in seq_along(g)) {
    if (Tn > j) L[(j + 1):Tn, ] <- L[(j + 1):Tn, ] + g[j] * I[1:(Tn - j), , drop = FALSE]
  }
  L
}

q <- function(x, p) quantile(x, p, names = FALSE, na.rm = TRUE)

# ══════════════════════════════════════════════════════════════════════════════
# 1. Model-based Rt from bootstrap trajectories
# ══════════════════════════════════════════════════════════════════════════════
manifest <- read.csv("manifest.csv", stringsAsFactors = FALSE)
groups   <- distinct(manifest, analysis, state, input_file)
rt_model <- list()

for (k in seq_len(nrow(groups))) {
  an <- groups$analysis[k]; st <- groups$state[k]
  dates <- read.csv(sub("\\.txt$", "_dates.csv", groups$input_file[k])) %>%
    mutate(week_start = as.Date(week_start))
  for (mdl in MODELS) {
    f <- file.path("runs", an, st, mdl, "curves_process.csv")
    if (!file.exists(f)) next
    I <- as.matrix(read.csv(f, header = FALSE))
    stopifnot(nrow(I) == nrow(dates))
    L <- lambda_mat(I)
    R <- I / L
    R[L <= 0] <- NA
    R[1:2, ] <- NA                    # unstable first steps, not reported
    rt_model[[length(rt_model) + 1]] <- tibble(
      analysis = an, state = st, model = mdl,
      week_start = dates$week_start, role = dates$role,
      median = apply(R, 1, q, .5),
      lower95 = apply(R, 1, q, .025),
      upper95 = apply(R, 1, q, .975)) %>%
      filter(!is.na(median))
  }
}
rt_model <- bind_rows(rt_model)

# ══════════════════════════════════════════════════════════════════════════════
# 2. Data-based Rt (Cori) from mature observed counts
# ══════════════════════════════════════════════════════════════════════════════
mature <- read_excel(file.path("data", "releases",
                               paste0("ensemble_confirmed_", MATURE_RELEASE, ".xlsx")),
                     sheet = "all_geographies") %>%
  mutate(week_start = as.Date(week_start))

rt_data <- list()
for (st in names(STATES)) {
  d <- mature %>% filter(geography == STATES[[st]], week_start <= LAST_WEEK) %>%
    arrange(week_start)
  I <- d$observed
  L <- as.vector(lambda_mat(matrix(I, ncol = 1)))
  for (t in seq_along(I)) {
    if (t < TAU + 2) next
    win <- (t - TAU + 1):t
    a <- PRIOR_SHAPE + sum(I[win])
    b <- 1 / (1 / PRIOR_SCALE + sum(L[win]))
    if (sum(L[win]) <= 0) next
    rt_data[[length(rt_data) + 1]] <- tibble(
      state = st, week_start = d$week_start[t],
      window_weeks = TAU, cases_in_window = sum(I[win]),
      mean = a * b, median = qgamma(.5, a, scale = b),
      lower95 = qgamma(.025, a, scale = b), upper95 = qgamma(.975, a, scale = b),
      p_above1 = 1 - pgamma(1, a, scale = b))
  }
}
rt_data <- bind_rows(rt_data)

# ══════════════════════════════════════════════════════════════════════════════
# Summary
# ══════════════════════════════════════════════════════════════════════════════
first_week <- function(df, cond) {
  w <- df$week_start[cond]; if (length(w)) min(w) else as.Date(NA)
}
summ_model <- rt_model %>%
  filter(role == "calibration" | TRUE) %>%
  group_by(analysis, state, model) %>%
  summarise(
    first_week_median_below1 = first_week(cur_data(), median < 1 & role == "calibration"),
    first_week_upper_below1  = first_week(cur_data(), upper95 < 1 & role == "calibration"),
    rt_calib_end = sprintf("%.2f (%.2f-%.2f)",
                           median[role == "calibration"][sum(role == "calibration")],
                           lower95[role == "calibration"][sum(role == "calibration")],
                           upper95[role == "calibration"][sum(role == "calibration")]),
    .groups = "drop")

summ_data <- rt_data %>%
  mutate(period = case_when(
    week_start >= as.Date("2026-06-15") & week_start <= as.Date("2026-07-20") ~ "primary",
    week_start >= as.Date("2026-07-27") ~ "extended",
    week_start >= as.Date("2026-05-11") & week_start <= as.Date("2026-06-08") ~ "pre_tournament_5wk",
    TRUE ~ NA_character_)) %>%
  filter(!is.na(period)) %>%
  group_by(state, period) %>%
  summarise(weeks = n(), median_of_weekly_medians = median(median),
            max_upper95 = max(upper95), weeks_p_above1_gt0.5 = sum(p_above1 > 0.5),
            .groups = "drop")

dir.create("tables", showWarnings = FALSE)
write.csv(rt_model,   "tables/rt_model.csv",   row.names = FALSE)
write.csv(rt_data,    "tables/rt_data.csv",    row.names = FALSE)
write.csv(summ_model, "tables/rt_summary_model.csv", row.names = FALSE)
write.csv(summ_data,  "tables/rt_summary_data.csv",  row.names = FALSE)

cat("\nModel-based Rt (main analysis):\n")
print(as.data.frame(filter(summ_model, analysis == "main")), row.names = FALSE)
cat("\nData-based Rt (Cori, ", TAU, "-week window), by period:\n", sep = "")
print(as.data.frame(summ_data), row.names = FALSE)
cat("\nWritten: tables/rt_model.csv, rt_data.csv, rt_summary_model.csv, rt_summary_data.csv\n")
