# 05_figures.R
# Mexico measles — publication figures (all labels set here; no manual editing)
#
#   figures/ (with titles) and figures/manuscript/ (no in-image titles, for the
#   journal; captions carry the text):
#   fig1_states_incidence        weekly incidence, all states (mature release)
#   figures/fig2_forecast_CDMX           fit + forecast, Mexico City (A: ensemble, full
#   figures/fig3_forecast_JAL            series; B-E: each model, zoom on forecast window)
#   figures/fig4_cumulative_excess       observed vs cumulative forecast, all analyses
#   figures/figS_rt                      model-based and data-based Rt
#   figures/figS_revisions               hold-out week counts across data releases
#
# Each figure is saved as PNG (300 dpi) and PDF. Run 03_evaluate_excess.R and
# 04_rt.R first.

suppressPackageStartupMessages({
  library(readxl); library(dplyr); library(tidyr); library(ggplot2)
  library(patchwork); library(scales)
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
dir.create("figures", showWarnings = FALSE)

# ══════════════════════════════════════════════════════════════════════════════
# Labels, colours, periods (edit here)
# ══════════════════════════════════════════════════════════════════════════════
MATURE_RELEASE <- "2026-09-21"
STATE_LABEL <- c(CDMX = "Mexico City", JAL = "Jalisco")
STATE_GEO   <- c(CDMX = "CIUDAD DE MEXICO", JAL = "JALISCO")
MODEL_LABEL <- c(GLM  = "Generalized logistic (GLM)",
                 GRM  = "Generalized Richards (GRM)",
                 RICH = "Richards (RICH)",
                 ENS  = "Ensemble (equal-weight linear pool)")
ANALYSIS_LABEL <- c(main        = "Main (mature counts)",
                    SA_realtime = "Calibrated on 13 Jul release",
                    SB_end0601  = "Calibration ends 1 Jun",
                    SC_Klow     = "Relaxed lower bound on K")
# Match periods at the host stadiums (FIFA schedule)
MATCHES <- tibble(state = c("CDMX", "JAL"),
                  start = as.Date(c("2026-06-11", "2026-06-11")),
                  end   = as.Date(c("2026-07-05", "2026-06-26")))
CALIB_END      <- as.Date("2026-06-08")
PRIMARY_END    <- as.Date("2026-07-20")
ZOOM_FROM      <- as.Date("2026-05-04")

COL_FIT   <- "#1f4e79"   # fit / forecast
COL_OBS   <- "black"
COL_EVAL  <- "#c0392b"   # observed counts in forecast weeks
COL_MATCH <- "#f2c94c"   # match period shading
COL_DATA_RT <- "#c0392b"

theme_pub <- function(base = 11) {
  theme_bw(base_size = base) +
    theme(panel.grid.minor = element_blank(),
          strip.background = element_rect(fill = "grey95", colour = "grey70"),
          strip.text = element_text(face = "bold"),
          plot.title = element_text(face = "bold"),
          plot.tag = element_text(face = "bold"),
          legend.position = "bottom")
}

save_fig <- function(p, name, w, h) {
  ggsave(file.path(OUT_DIR, paste0(name, ".png")), p, width = w, height = h, dpi = 300)
  ggsave(file.path(OUT_DIR, paste0(name, ".pdf")), p, width = w, height = h)
  cat("  saved ", OUT_DIR, "/", name, " (.png, .pdf)\n", sep = "")
}
# Titles/subtitles inside the image are dropped in manuscript mode (the
# caption carries them); panel titles and tags are kept.
ttl <- function(x) if (MS) NULL else x

q <- function(x, p) quantile(x, p, names = FALSE)

# Curves -> weekly bands (predictive PI from noisy draws; line = process median)
read_bands <- function(analysis, state, model) {
  dates <- read.csv(file.path("inputs", analysis,
                              sprintf("%s_%s_dates.csv", state, analysis))) %>%
    mutate(week_start = as.Date(week_start))
  base <- file.path("runs", analysis, state, model)
  P  <- as.matrix(read.csv(file.path(base, "curves_predictive.csv"), header = FALSE))
  Pr <- as.matrix(read.csv(file.path(base, "curves_process.csv"),    header = FALSE))
  dates %>% mutate(
    analysis = analysis, state = state, model = model,
    fit     = apply(Pr, 1, q, .5),
    lower95 = apply(P, 1, q, .025), upper95 = apply(P, 1, q, .975),
    lower50 = apply(P, 1, q, .25),  upper50 = apply(P, 1, q, .75),
    obs_type = ifelse(role == "calibration", "Calibration data",
                      "Observed (evaluation)"))
}

for (MS in c(FALSE, TRUE)) {
OUT_DIR <- if (MS) file.path("figures", "manuscript") else "figures"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
cat(if (MS) "\n-- Manuscript versions (no in-image titles) --\n" else "-- Annotated versions --\n")

# ══════════════════════════════════════════════════════════════════════════════
# Figure 1: weekly incidence by state
# ══════════════════════════════════════════════════════════════════════════════
cat("Figure 1\n")
# Preferred source: the cleaned line list of the mature release (all states).
# It is not in git (patient-level data); without it the release file is used,
# which lacks states where the nowcast step failed (fixed in run_pipeline_loop.R
# for later releases).
LL_FILE <- file.path("..", "outputs", paste0("02_cleaned_data_", MATURE_RELEASE),
                     "measles_confirmed.rds")
if (file.exists(LL_FILE)) {
  mature <- readRDS(LL_FILE) %>%
    filter(!is.na(FEC_INI_EXANT)) %>%
    mutate(week_start = lubridate::floor_date(as.Date(FEC_INI_EXANT), "week",
                                              week_start = 1)) %>%
    count(geography = state_notif, week_start, name = "observed")
  cat("  source: cleaned line list (", LL_FILE, ")\n", sep = "")
} else {
  warning("Cleaned line list not found (", LL_FILE, "); Figure 1 uses the release ",
          "file, which may omit states without a nowcast.")
  mature <- read_excel(file.path("data", "releases",
                                 paste0("ensemble_confirmed_", MATURE_RELEASE, ".xlsx")),
                       sheet = "all_geographies") %>%
    mutate(week_start = as.Date(week_start)) %>%
    filter(geography != "National")
}

tot <- mature %>% group_by(geography) %>% summarise(total = sum(observed)) %>%
  arrange(desc(total))
f1 <- mature %>%
  mutate(label = tools::toTitleCase(tolower(geography)),
         label = recode(label, "Ciudad De Mexico" = "Mexico City"),
         host  = geography %in% c(STATE_GEO, "NUEVO LEON")) %>%
  left_join(tot, by = "geography") %>%
  mutate(label = factor(label, levels = unique(label[order(-total)])))

p1 <- ggplot(f1, aes(week_start, observed, fill = host)) +
  annotate("rect", xmin = as.Date("2026-06-11"), xmax = as.Date("2026-07-05"),
           ymin = -Inf, ymax = Inf, fill = COL_MATCH, alpha = .35) +
  geom_col(width = 6) +
  facet_wrap(~ label, scales = "free_y", ncol = 5) +
  scale_fill_manual(values = c(`TRUE` = COL_EVAL, `FALSE` = "grey45"),
                    labels = c(`TRUE` = "World Cup host jurisdiction",
                               `FALSE` = "Other state"), name = NULL) +
  scale_x_date(date_breaks = "4 months", date_labels = "%b\n%Y") +
  labs(title = ttl("Weekly confirmed measles cases by state of notification, Mexico, 2025-2026"),
       subtitle = ttl(sprintf("Week of rash onset; DGE release of %s. Shading: World Cup matches in Mexico (11 Jun to 5 Jul 2026). Free y-axes.",
                          format(as.Date(MATURE_RELEASE), "%d %b %Y"))),
       x = "Week of rash onset", y = "Weekly confirmed cases") +
  theme_pub(10)
save_fig(p1, "fig1_states_incidence", if (MS) 10 else 15, if (MS) 11 else 15)

# ══════════════════════════════════════════════════════════════════════════════
# Figures 2-3 — fit and forecast per jurisdiction
# ══════════════════════════════════════════════════════════════════════════════
forecast_panel <- function(b, st, zoom = TRUE, title = NULL) {
  m <- MATCHES %>% filter(state == st)
  if (zoom) b <- b %>% filter(week_start >= ZOOM_FROM)
  ggplot(b, aes(week_start)) +
    annotate("rect", xmin = m$start, xmax = m$end, ymin = -Inf, ymax = Inf,
             fill = COL_MATCH, alpha = .35) +
    geom_vline(xintercept = CALIB_END + 6.5, linetype = "dashed", colour = "grey30") +
    geom_vline(xintercept = PRIMARY_END + 6.5, linetype = "dotted", colour = "grey30") +
    geom_ribbon(aes(ymin = lower95, ymax = upper95), fill = COL_FIT, alpha = .15) +
    geom_ribbon(aes(ymin = lower50, ymax = upper50), fill = COL_FIT, alpha = .25) +
    geom_line(aes(y = fit), colour = COL_FIT, linewidth = .8) +
    geom_point(aes(y = cases, shape = obs_type, colour = obs_type), size = 1.9) +
    scale_shape_manual(values = c("Calibration data" = 16, "Observed (evaluation)" = 17), name = NULL) +
    scale_colour_manual(values = c("Calibration data" = COL_OBS, "Observed (evaluation)" = COL_EVAL), name = NULL) +
    scale_x_date(date_breaks = if (zoom) "1 month" else "2 months",
                 date_labels = "%d %b\n%Y") +
    scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, .05))) +
    labs(title = title, x = "Week of rash onset", y = "Weekly confirmed cases") +
    theme_pub(11)
}

for (st in names(STATE_LABEL)) {
  fig_no <- c(CDMX = 2, JAL = 3)[[st]]
  cat("Figure", fig_no, "\n")
  bands <- lapply(names(MODEL_LABEL), function(mdl) read_bands("main", st, mdl))
  names(bands) <- names(MODEL_LABEL)
  pA <- forecast_panel(bands$ENS, st, zoom = FALSE,
                       title = paste0(MODEL_LABEL[["ENS"]], ": full calibration window"))
  pB <- lapply(names(MODEL_LABEL), function(mdl)
    forecast_panel(bands[[mdl]], st, zoom = TRUE,
                   title = paste0(MODEL_LABEL[[mdl]], ": from ", format(ZOOM_FROM, "%d %b"))))
  p <- (pA / wrap_plots(pB, ncol = 2)) +
    plot_layout(heights = c(1, 2), guides = "collect") +
    plot_annotation(
      tag_levels = "A",
      title = ttl(sprintf("%s: forecast benchmark from pre-tournament data and observed measles incidence", STATE_LABEL[[st]])),
      subtitle = ttl(paste0("Line: model fit and forecast (median of bootstrap trajectories); bands: 50% and 95% prediction intervals (300 bootstraps, with observation noise).\n",
                        "Dashed line: end of calibration (onset week of 8 Jun). Dotted line: end of primary evaluation window (onset week of 20 Jul). Shading: matches hosted in ",
                        STATE_LABEL[[st]], ".")),
      theme = theme(plot.title = element_text(face = "bold", size = 14),
                    legend.position = "bottom"))
  save_fig(p, sprintf("fig%d_forecast_%s", fig_no, st), 13, 12)
}

# ══════════════════════════════════════════════════════════════════════════════
# Figure 4 — cumulative excess, all analyses
# ══════════════════════════════════════════════════════════════════════════════
cat("Figure 4\n")
cx <- read.csv("tables/cumulative_excess.csv") %>%
  mutate(state_l = factor(STATE_LABEL[state], levels = STATE_LABEL),
         window_l = factor(ifelse(window == "primary",
                                  "Primary window (onsets 15 Jun to 26 Jul)",
                                  "Full window (onsets 15 Jun to 9 Aug)"),
                           levels = c("Primary window (onsets 15 Jun to 26 Jul)",
                                      "Full window (onsets 15 Jun to 9 Aug)")),
         analysis_l = factor(ANALYSIS_LABEL[analysis], levels = rev(ANALYSIS_LABEL)),
         model = factor(model, levels = rev(names(MODEL_LABEL))))
obs_lines <- cx %>% distinct(state_l, window_l, observed)

p4 <- ggplot(cx, aes(y = analysis_l, colour = model)) +
  geom_vline(data = obs_lines, aes(xintercept = observed), colour = COL_EVAL, linewidth = .8) +
  geom_linerange(aes(xmin = lower95, xmax = upper95),
                 position = position_dodge(width = .7), linewidth = .7) +
  geom_point(aes(x = median), position = position_dodge(width = .7), size = 2) +
  facet_wrap(state_l ~ window_l, scales = "free_x", ncol = 2) +
  scale_colour_manual(values = c(GLM = "#7f8c8d", GRM = "#16a085", RICH = "#8e44ad", ENS = COL_FIT),
                      labels = MODEL_LABEL, breaks = names(MODEL_LABEL), name = NULL) +
  labs(title = ttl("Observed cumulative cases versus the forecast benchmark"),
       subtitle = ttl("Points and bars: median and 95% prediction interval of cumulative forecast cases. Red line: observed cumulative cases (mature counts)."),
       x = "Cumulative confirmed cases", y = NULL) +
  theme_pub(11) + guides(colour = guide_legend(nrow = 2))
save_fig(p4, "fig4_cumulative_excess", 12, 7.5)

# ══════════════════════════════════════════════════════════════════════════════
# Figure S — reproduction number
# ══════════════════════════════════════════════════════════════════════════════
cat("Figure S (Rt)\n")
rtm <- read.csv("tables/rt_model.csv") %>% mutate(week_start = as.Date(week_start)) %>%
  filter(analysis == "main", model == "ENS")
rtd <- read.csv("tables/rt_data.csv") %>% mutate(week_start = as.Date(week_start))
RT_FROM <- as.Date("2026-01-05")

rt_panel <- function(st) {
  m <- MATCHES %>% filter(state == st)
  ggplot() +
    annotate("rect", xmin = m$start, xmax = m$end, ymin = -Inf, ymax = Inf,
             fill = COL_MATCH, alpha = .35) +
    geom_hline(yintercept = 1, linetype = "dashed") +
    geom_vline(xintercept = CALIB_END + 6.5, linetype = "dashed", colour = "grey30") +
    geom_ribbon(data = filter(rtm, state == st, week_start >= RT_FROM),
                aes(week_start, ymin = lower95, ymax = upper95), fill = COL_FIT, alpha = .2) +
    geom_line(data = filter(rtm, state == st, week_start >= RT_FROM),
              aes(week_start, median, linetype = "Model-based (ensemble fit; projection after calibration)"),
              colour = COL_FIT, linewidth = .8) +
    geom_pointrange(data = filter(rtd, state == st, week_start >= RT_FROM),
                    aes(week_start, median, ymin = lower95, ymax = upper95,
                        shape = "Data-based (Cori, 2-week window)"),
                    colour = COL_DATA_RT, size = .3) +
    scale_linetype_manual(values = 1, name = NULL) +
    scale_shape_manual(values = 16, name = NULL) +
    scale_x_date(date_breaks = "1 month", date_labels = "%b\n%Y") +
    coord_cartesian(ylim = c(0, 3.5)) +
    labs(title = STATE_LABEL[[st]], x = "Week of rash onset (end of estimation window)",
         y = expression(R[t])) +
    theme_pub(11)
}
pS <- (rt_panel("CDMX") / rt_panel("JAL")) + plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "A",
    title = ttl("Effective reproduction number before and during the World Cup"),
    subtitle = ttl("Generation interval: gamma, mean 11.7 days, SD 2 days. Bands/bars: 95% intervals. Shading: matches hosted in each jurisdiction."),
    theme = theme(plot.title = element_text(face = "bold", size = 14), legend.position = "bottom"))
save_fig(pS, "figS_rt", 11, 9)

# ══════════════════════════════════════════════════════════════════════════════
# Figure S — data revisions across releases
# ══════════════════════════════════════════════════════════════════════════════
cat("Figure S (revisions)\n")
rv <- read.csv("tables/data_revisions.csv") %>%
  mutate(week_start = as.Date(week_start), release = as.Date(release),
         state_l = factor(STATE_LABEL[state], levels = STATE_LABEL)) %>%
  filter(week_start >= as.Date("2026-06-01"), week_start <= as.Date("2026-07-20"))
pR <- ggplot(rv, aes(release, observed, colour = format(week_start, "%d %b"))) +
  geom_line() + geom_point(size = 1.6) +
  facet_wrap(~ state_l, scales = "free_y") +
  scale_x_date(breaks = sort(unique(rv$release)), date_labels = "%d %b") +
  labs(title = ttl("Reported weekly counts by data release"),
       subtitle = ttl("Each line is one onset week; counts rise as delayed reports are added."),
       x = "Data release", y = "Confirmed cases in onset week",
       colour = "Onset week beginning") +
  theme_pub(11) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
save_fig(pR, "figS_revisions", 12, 5.5)

}

cat("\nAll figures written to figures/ and figures/manuscript/\n")
