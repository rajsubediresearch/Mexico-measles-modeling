# 01_prepare_inputs.R
# Mexico measles — World Cup 2026 forecast benchmark (QuantDiffForecast)
#
# Builds everything the MATLAB driver needs, from the weekly release files in
# data/releases/ (ensemble_confirmed_<release>.xlsx, produced by the nowcasting
# pipeline):
#   inputs/<analysis>/<STATE>_<analysis>.txt     toolbox input (t, cases)
#   inputs/<analysis>/<STATE>_<analysis>_dates.csv  t -> week, role, raw/nowcast
#   options/<analysis>/options_<analysis>_<STATE>_<MODEL>.m
#   manifest.csv                                  one row per MATLAB run
#   tables/data_revisions.csv                     hold-out counts by release
#
# Run from the forecasting/ folder (or set PROJ below). Re-running overwrites
# the generated files; nothing in data/ is modified.

suppressPackageStartupMessages({
  library(readxl); library(dplyr); library(tidyr); library(lubridate)
})

PROJ <- if (interactive() && requireNamespace("rstudioapi", quietly = TRUE) &&
            rstudioapi::isAvailable()) {
  dirname(rstudioapi::getSourceEditorContext()$path)
} else {
  args <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", args[grepl("^--file=", args)])
  if (length(f)) dirname(normalizePath(f)) else getwd()
}
setwd(PROJ)

# ══════════════════════════════════════════════════════════════════════════════
# CONFIGURATION
# ══════════════════════════════════════════════════════════════════════════════

MATURE_RELEASE <- "2026-09-21"   # frozen analysis release (mature counts)

STATES <- tribble(
  ~code,  ~geography,          ~label,        ~calib_start,
  "CDMX", "CIUDAD DE MEXICO",  "Mexico City", "2025-11-24",
  "JAL",  "JALISCO",           "Jalisco",     "2025-09-08"
) %>% mutate(calib_start = as.Date(calib_start))

# Forecast window: onset weeks 15 Jun - 3 Aug 2026 (8 weeks)
#   primary  = weeks beginning 15 Jun - 20 Jul (direct tournament exposures)
#   extended = weeks beginning 27 Jul - 3 Aug (one further generation)
FIRST_EVAL_WEEK    <- as.Date("2026-06-15")   # common to all analyses
LAST_FORECAST_WEEK <- as.Date("2026-08-03")
PRIMARY_LAST_WEEK  <- as.Date("2026-07-20")

# k_lb_mult: lower bound of K as a multiple of cumulative calibration cases
ANALYSES <- tribble(
  ~analysis,     ~calib_end,   ~calib_release, ~calib_value, ~k_lb_mult, ~description,
  "main",        "2026-06-08", MATURE_RELEASE, "observed",   1,
    "Mature counts; calibration ends with onset week 8-14 Jun",
  "SA_realtime", "2026-06-08", "2026-07-13",   "v8",         1,
    "Calibration on 13 Jul release as used in v8 (nowcast-adjusted); evaluated on mature counts",
  "SB_end0601",  "2026-06-01", MATURE_RELEASE, "observed",   1,
    "Calibration ends one week earlier (onset week 1-7 Jun)",
  "SC_Klow",     "2026-06-08", MATURE_RELEASE, "observed",   0.5,
    "As main, with the lower bound of K relaxed to half the cumulative calibration cases"
) %>% mutate(calib_end = as.Date(calib_end))

MODELS <- c("GLM", "GRM", "RICH")

# Estimation settings (shared by all runs)
METHOD1        <- 3     # MLE, negative binomial VAR = mean + alpha*mean (NB1)
NUMSTARTPOINTS <- 100
B              <- 300
K_UB_MULT      <- 10    # K upper bound = 10 x cumulative calibration cases

# ══════════════════════════════════════════════════════════════════════════════
# HELPERS
# ══════════════════════════════════════════════════════════════════════════════

read_release <- function(release) {
  f <- file.path("data", "releases", paste0("ensemble_confirmed_", release, ".xlsx"))
  if (!file.exists(f)) stop("Missing release file: ", f)
  read_excel(f, sheet = "all_geographies") %>%
    mutate(week_start = as.Date(week_start), release = release)
}

# value used by the v8 / nowcasting pipeline: nowcast where it exceeds observed
v8_value <- function(observed, nowcast) {
  ifelse(!is.na(nowcast) & nowcast > observed, nowcast, observed)
}

state_series <- function(rel_df, geo) {
  rel_df %>% filter(geography == geo) %>%
    select(week_start, observed, nowcast_estimate) %>%
    arrange(week_start)
}

fmt <- function(x) format(signif(x, 6), scientific = FALSE, trim = TRUE)
vec <- function(x) paste0("[", paste(sapply(x, fmt), collapse = " "), "]")
cellstr <- function(x) paste0("{", paste0("'", x, "'", collapse = ","), "}")

model_spec <- function(model, cum_cases, k_lb_mult = 1) {
  K_lb <- k_lb_mult * cum_cases; K_ub <- K_UB_MULT * cum_cases; K_init <- 1.2 * cum_cases
  switch(model,
    GLM  = list(fc = "GLM",  name = "GLM model",
                label = c("r","p","K"),
                LB = c(0.001, 0.01, K_lb), UB = c(20, 1, K_ub),
                init = c(0.5, 0.8, K_init)),
    GRM  = list(fc = "GRM",  name = "GRM model",
                label = c("r","p","a","K"),
                LB = c(0.001, 0.01, 0.01, K_lb), UB = c(20, 1, 10, K_ub),
                init = c(0.5, 0.8, 1, K_init)),
    RICH = list(fc = "RICH", name = "Richards model",
                label = c("r","a","K"),
                LB = c(0.001, 0.01, K_lb), UB = c(20, 10, K_ub),
                init = c(0.5, 1, K_init)),
    stop("Unknown model ", model))
}

write_options <- function(path, fn_name, cadfilename, spec, windowsize, horizon,
                          analysis, state_label, cal_start, cal_end) {
  txt <- c(
    sprintf("function [cadfilename1,caddisease,datatype,dist1,numstartpoints,B,model,params,vars,getperformance,forecastingperiod,windowsize1,tstart1,tend1,printscreen1] = %s", fn_name),
    "% Generated by 01_prepare_inputs.R - edit the R config, not this file.",
    sprintf("%% Analysis: %s | %s | %s", analysis, state_label, spec$name),
    sprintf("%% Calibration: onset weeks %s to %s (%d weeks); forecast horizon %d weeks",
            cal_start, cal_end, windowsize, horizon),
    "",
    "global method1",
    "",
    sprintf("cadfilename1 = '%s';", cadfilename),
    "caddisease   = 'measles';",
    "datatype     = 'cases';",
    "",
    sprintf("method1 = %d;   %% MLE negative binomial, VAR = mean + alpha*mean (NB1)", METHOD1),
    sprintf("dist1   = %d;", METHOD1),
    "",
    sprintf("numstartpoints = %d;", NUMSTARTPOINTS),
    sprintf("B              = %d;", B),
    "",
    sprintf("model.fc   = @%s;", spec$fc),
    sprintf("model.name = '%s';", spec$name),
    "",
    sprintf("params.num     = %d;", length(spec$label)),
    sprintf("params.label   = %s;", cellstr(spec$label)),
    sprintf("params.LB      = %s;", vec(spec$LB)),
    sprintf("params.UB      = %s;", vec(spec$UB)),
    sprintf("params.initial = %s;", vec(spec$init)),
    sprintf("params.fixed   = %s;", vec(rep(0, length(spec$label)))),
    "params.fixI0   = 1;   % C(0) fixed to the first observation",
    "params.composite      = [];",
    "params.composite_name = '';",
    "params.extra0         = [];",
    "",
    "vars.num       = 1;",
    "vars.label     = {'C'};",
    "vars.initial   = 1;   % overridden because params.fixI0 = 1",
    "vars.fit_index = 1;",
    "vars.fit_diff  = 1;   % fit incidence (derivative of C)",
    "",
    "getperformance    = 1;",
    sprintf("forecastingperiod = %d;", horizon),
    sprintf("windowsize1       = %d;", windowsize),
    "tstart1           = 1;",
    "tend1             = 1;",
    "printscreen1      = 0;   % figures are made by 05_figures.R",
    "end"
  )
  writeLines(txt, path)
}

# ══════════════════════════════════════════════════════════════════════════════
# BUILD
# ══════════════════════════════════════════════════════════════════════════════

releases <- sub("^ensemble_confirmed_(.*)\\.xlsx$", "\\1",
                list.files("data/releases", pattern = "^ensemble_confirmed_.*\\.xlsx$"))
rel_all <- bind_rows(lapply(releases, read_release))
mature  <- rel_all %>% filter(release == MATURE_RELEASE)

dir.create("inputs",  showWarnings = FALSE)
dir.create("options", showWarnings = FALSE)
dir.create("tables",  showWarnings = FALSE)

manifest <- list()
summary_rows <- list()

for (a in seq_len(nrow(ANALYSES))) {
  an <- ANALYSES[a, ]
  calib_rel <- rel_all %>% filter(release == an$calib_release)
  if (nrow(calib_rel) == 0) stop("Release not found: ", an$calib_release)
  in_dir  <- file.path("inputs",  an$analysis); dir.create(in_dir,  showWarnings = FALSE)
  opt_dir <- file.path("options", an$analysis); dir.create(opt_dir, showWarnings = FALSE)

  for (s in seq_len(nrow(STATES))) {
    st <- STATES[s, ]
    weeks <- seq(st$calib_start, LAST_FORECAST_WEEK, by = "week")

    mat <- state_series(mature,    st$geography) %>% rename(obs_mature = observed, nowcast_mature = nowcast_estimate)
    cal <- state_series(calib_rel, st$geography) %>% rename(obs_calrel = observed, nowcast_calrel = nowcast_estimate)

    d <- tibble(week_start = weeks) %>%
      left_join(mat, by = "week_start") %>%
      left_join(cal, by = "week_start") %>%
      mutate(t    = row_number() - 1L,
             role = case_when(week_start <= an$calib_end ~ "calibration",
                              # forecast weeks before the common evaluation
                              # window (only when calibration ends earlier)
                              week_start <  FIRST_EVAL_WEEK ~ "forecast_lead",
                              week_start <= PRIMARY_LAST_WEEK ~ "forecast_primary",
                              TRUE ~ "forecast_extended"),
             calib_value = if (an$calib_value == "v8")
                             v8_value(obs_calrel, nowcast_calrel) else obs_calrel,
             # calibration weeks: value from the calibration release;
             # forecast weeks: mature observed counts (evaluation data)
             cases = ifelse(role == "calibration", calib_value, obs_mature))

    if (any(is.na(d$cases))) {
      stop(sprintf("%s/%s: missing counts for weeks %s", an$analysis, st$code,
                   paste(d$week_start[is.na(d$cases)], collapse = ", ")))
    }
    if (d$cases[1] == 0) stop(st$code, ": first calibration week is zero")

    windowsize <- sum(d$role == "calibration")
    horizon    <- sum(d$role != "calibration")
    cum_cal    <- sum(d$cases[d$role == "calibration"])

    cad <- sprintf("%s_%s", st$code, an$analysis)
    write.table(d %>% transmute(t, cases = round(cases, 4)),
                file.path(in_dir, paste0(cad, ".txt")),
                sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
    write.csv(d %>% select(t, week_start, role, cases, obs_mature, nowcast_mature,
                           obs_calrel, nowcast_calrel),
              file.path(in_dir, paste0(cad, "_dates.csv")), row.names = FALSE)

    summary_rows[[length(summary_rows) + 1]] <- tibble(
      analysis = an$analysis, state = st$code,
      calib_start = st$calib_start, calib_end = an$calib_end,
      windowsize = windowsize, horizon = horizon, cum_calibration = cum_cal,
      non_integer_calib = sum(d$cases != round(d$cases)))

    for (m in seq_along(MODELS)) {
      spec <- model_spec(MODELS[m], cum_cal, an$k_lb_mult)
      fn   <- sprintf("options_%s_%s_%s", an$analysis, st$code, MODELS[m])
      write_options(file.path(opt_dir, paste0(fn, ".m")), fn, cad, spec,
                    windowsize, horizon, an$analysis, st$label,
                    st$calib_start, an$calib_end)
      manifest[[length(manifest) + 1]] <- tibble(
        run_id = sprintf("%s__%s__%s", an$analysis, st$code, MODELS[m]),
        analysis = an$analysis, state = st$code, model = MODELS[m],
        input_file = file.path("inputs", an$analysis, paste0(cad, ".txt")),
        options_dir = opt_dir, options_fn = fn,
        windowsize = windowsize, horizon = horizon,
        seed = 1000L * a + 10L * s + m,
        out_dir = file.path("runs", an$analysis, st$code, MODELS[m]))
    }
  }
}

manifest <- bind_rows(manifest)
write.csv(manifest, "manifest.csv", row.names = FALSE)

# ── Data revisions: hold-out weeks across releases (comment 1 / S-A) ─────────
revisions <- rel_all %>%
  inner_join(STATES %>% select(code, geography), by = "geography") %>%
  filter(week_start >= as.Date("2026-05-25"), week_start <= LAST_FORECAST_WEEK) %>%
  transmute(state = code, week_start, release, observed, nowcast_estimate) %>%
  arrange(state, week_start, release)
write.csv(revisions, "tables/data_revisions.csv", row.names = FALSE)

# ── Report ───────────────────────────────────────────────────────────────────
smry <- bind_rows(summary_rows)
cat("\nInputs prepared (mature release:", MATURE_RELEASE, ")\n")
print(as.data.frame(smry), row.names = FALSE)
if (any(smry$non_integer_calib > 0)) {
  cat("\nNote: non-integer calibration counts present (nowcast values) in:",
      paste(smry$analysis[smry$non_integer_calib > 0],
            smry$state[smry$non_integer_calib > 0], collapse = "; "), "\n")
}
cat(sprintf("\n%d runs written to manifest.csv\n", nrow(manifest)))
