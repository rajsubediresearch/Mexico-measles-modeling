# 06_submission_figures.R
# Mexico measles — main-text figures sized for Journal of Travel Medicine submission
#
#   Figure 1  weekly confirmed cases by state of notification (mature release)
#   Figure 2  fit + forecast, Mexico City (A: ensemble, full calibration window;
#             B–E: each model and the ensemble, from May 2026)
#   Figure 3  same for Jalisco
#
# JTM asks for figures as separate files (not embedded in Word), line art at
# 600 dpi and 13 cm wide. Each figure is written as
#   FigureN.tif  600 dpi, LZW-compressed (upload file)
#   FigureN.pdf  vector copy (for checking; the .tif files are the uploads)
# Content matches 05_figures.R (same data, models and colours); only size, font
# size and label wording (British spelling, shorter panel titles) differ.
# Supplementary figures stay embedded in the supplement and are not exported here.
#
# Run 03_evaluate_excess.R and 04_rt.R first (as for 05_figures.R).
# Output folder: set SUBMIT_DIR below, or pass it as an environment variable:
#   SUBMIT_DIR="/path/to/final figures" Rscript 06_submission_figures.R

suppressPackageStartupMessages({
  library(readxl); library(dplyr); library(tidyr); library(ggplot2)
  library(patchwork); library(scales)
})

# Folder of this script: sourced file > RStudio/Positron editor > Rscript --file > working directory
PROJ <- local({
  of <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
  if (!is.null(of) && nzchar(of)) return(dirname(normalizePath(of)))
  if (interactive() && requireNamespace("rstudioapi", quietly = TRUE) &&
      rstudioapi::isAvailable()) {
    ed <- tryCatch(rstudioapi::getSourceEditorContext()$path, error = function(e) NULL)
    if (length(ed) == 1 && nzchar(ed) && basename(ed) == "06_submission_figures.R")
      return(dirname(ed))
  }
  args <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", args[grepl("^--file=", args)])
  if (length(f)) dirname(normalizePath(f)) else getwd()
})
setwd(PROJ)

SUBMIT_DIR <- Sys.getenv(
  "SUBMIT_DIR",
  "/Users/rajsubedi/Library/CloudStorage/Dropbox-GSUDropbox/Raj Subedi/Mexico_Measles_manuscript/Final Manuscript from Dr Chowell/final figures")
dir.create(SUBMIT_DIR, showWarnings = FALSE, recursive = TRUE)

# ══════════════════════════════════════════════════════════════════════════════
# Settings (edit here)
# ══════════════════════════════════════════════════════════════════════════════
WIDTH_CM <- 13       # JTM: 13 cm wide
DPI      <- 600      # JTM: line art at 600 dpi
BASE_PT  <- 7        # base font size at final print size

MATURE_RELEASE <- "2026-09-21"
STATE_LABEL <- c(CDMX = "Mexico City", JAL = "Jalisco")
STATE_GEO   <- c(CDMX = "CIUDAD DE MEXICO", JAL = "JALISCO")
MODEL_LABEL <- c(GLM  = "Generalised logistic (GLM)",
                 GRM  = "Generalised Richards (GRM)",
                 RICH = "Richards (RICH)",
                 ENS  = "Ensemble")
MATCHES <- tibble(state = c("CDMX", "JAL"),
                  start = as.Date(c("2026-06-11", "2026-06-11")),
                  end   = as.Date(c("2026-07-05", "2026-06-26")))
CALIB_END   <- as.Date("2026-06-08")
PRIMARY_END <- as.Date("2026-07-20")
ZOOM_FROM   <- as.Date("2026-05-04")

COL_FIT   <- "#1f4e79"
COL_OBS   <- "black"
COL_EVAL  <- "#c0392b"
COL_MATCH <- "#f2c94c"

theme_sub <- function(base = BASE_PT) {
  theme_bw(base_size = base) +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major = element_line(linewidth = .2),
          panel.border = element_rect(linewidth = .3),
          axis.ticks = element_line(linewidth = .25),
          strip.background = element_rect(fill = "grey95", colour = "grey70", linewidth = .3),
          strip.text = element_text(face = "bold", size = base - 1,
                                    margin = margin(1.5, 1, 1.5, 1)),
          plot.title = element_text(face = "bold", size = base),
          plot.tag = element_text(face = "bold", size = base + 1),
          legend.position = "bottom",
          legend.key.size = unit(3, "mm"),
          legend.margin = margin(0, 0, 0, 0),
          plot.margin = margin(3, 4, 3, 3))
}

save_sub <- function(p, name, height_cm) {
  w <- WIDTH_CM / 2.54; h <- height_cm / 2.54
  tif <- file.path(SUBMIT_DIR, paste0(name, ".tif"))
  pdf <- file.path(SUBMIT_DIR, paste0(name, ".pdf"))
  if (requireNamespace("ragg", quietly = TRUE)) {
    ggsave(tif, p, width = w, height = h, dpi = DPI, device = ragg::agg_tiff,
           compression = "lzw", bg = "white")
  } else {
    ggsave(tif, p, width = w, height = h, dpi = DPI, device = "tiff",
           compression = "lzw", bg = "white")
  }
  # base PDF device: works on macOS without XQuartz (cairo_pdf needs it)
  ggsave(pdf, p, width = w, height = h, device = grDevices::pdf, useDingbats = FALSE)
  cat(sprintf("  %s: %.0f x %.0f px (%g x %g cm at %d dpi) + PDF\n",
              name, w * DPI, h * DPI, WIDTH_CM, height_cm, DPI))
}

q <- function(x, p) quantile(x, p, names = FALSE)

read_bands <- function(analysis, state, model) {
  dates <- read.csv(file.path("inputs", analysis,
                              sprintf("%s_%s_dates.csv", state, analysis))) %>%
    mutate(week_start = as.Date(week_start))
  base <- file.path("runs", analysis, state, model)
  P  <- as.matrix(read.csv(file.path(base, "curves_predictive.csv"), header = FALSE))
  Pr <- as.matrix(read.csv(file.path(base, "curves_process.csv"),    header = FALSE))
  dates %>% mutate(
    fit     = apply(Pr, 1, q, .5),
    lower95 = apply(P, 1, q, .025), upper95 = apply(P, 1, q, .975),
    lower50 = apply(P, 1, q, .25),  upper50 = apply(P, 1, q, .75),
    obs_type = ifelse(role == "calibration", "Calibration data",
                      "Observed (forecast weeks)"))
}

cat("Writing submission figures to:\n  ", SUBMIT_DIR, "\n", sep = "")

# ══════════════════════════════════════════════════════════════════════════════
# Figure 1 — weekly incidence by state
# ══════════════════════════════════════════════════════════════════════════════
LL_FILE <- file.path("..", "outputs", paste0("02_cleaned_data_", MATURE_RELEASE),
                     "measles_confirmed.rds")
if (file.exists(LL_FILE)) {
  mature <- readRDS(LL_FILE) %>%
    filter(!is.na(FEC_INI_EXANT)) %>%
    mutate(week_start = lubridate::floor_date(as.Date(FEC_INI_EXANT), "week",
                                              week_start = 1)) %>%
    count(geography = state_notif, week_start, name = "observed")
} else {
  warning("Cleaned line list not found (", LL_FILE, "); Figure 1 uses the release ",
          "file, which may omit states without a nowcast. Use the line list for submission.")
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
           ymin = -Inf, ymax = Inf, fill = COL_MATCH, alpha = .45) +
  geom_col(width = 6) +
  facet_wrap(~ label, scales = "free_y", ncol = 5) +
  scale_fill_manual(values = c(`TRUE` = COL_EVAL, `FALSE` = "grey45"),
                    breaks = c("TRUE", "FALSE"),
                    labels = c(`TRUE` = "World Cup host jurisdiction",
                               `FALSE` = "Other state"), name = NULL) +
  scale_x_date(date_breaks = "6 months", date_labels = "%b\n%Y") +
  scale_y_continuous(breaks = breaks_pretty(n = 3)) +
  labs(x = "Week of rash onset", y = "Weekly confirmed cases") +
  theme_sub() +
  theme(axis.text = element_text(size = BASE_PT - 1.5),
        strip.text = element_text(size = BASE_PT - 1.5, face = "bold",
                                  margin = margin(1, 1, 1, 1)),
        panel.spacing = unit(1.2, "mm"))
n_rows <- ceiling(nlevels(f1$label) / 5)
save_sub(p1, "Figure1", min(22, 2.4 * n_rows + 2))   # ~19 cm for 32 states

# ══════════════════════════════════════════════════════════════════════════════
# Figures 2–3 — fit and forecast per jurisdiction
# ══════════════════════════════════════════════════════════════════════════════
forecast_panel <- function(b, st, zoom = TRUE, title = NULL) {
  m <- MATCHES %>% filter(state == st)
  if (zoom) b <- b %>% filter(week_start >= ZOOM_FROM)
  ggplot(b, aes(week_start)) +
    annotate("rect", xmin = m$start, xmax = m$end, ymin = -Inf, ymax = Inf,
             fill = COL_MATCH, alpha = .45) +
    geom_vline(xintercept = CALIB_END + 6.5, linetype = "dashed",
               colour = "grey30", linewidth = .3) +
    geom_vline(xintercept = PRIMARY_END + 6.5, linetype = "dotted",
               colour = "grey30", linewidth = .4) +
    geom_ribbon(aes(ymin = lower95, ymax = upper95), fill = COL_FIT, alpha = .15) +
    geom_ribbon(aes(ymin = lower50, ymax = upper50), fill = COL_FIT, alpha = .25) +
    geom_line(aes(y = fit), colour = COL_FIT, linewidth = .5) +
    geom_point(aes(y = cases, shape = obs_type, colour = obs_type), size = .9) +
    scale_shape_manual(values = c("Calibration data" = 16, "Observed (forecast weeks)" = 17),
                       name = NULL) +
    scale_colour_manual(values = c("Calibration data" = COL_OBS,
                                   "Observed (forecast weeks)" = COL_EVAL), name = NULL) +
    scale_x_date(date_breaks = if (zoom) "1 month" else "2 months",
                 date_labels = if (zoom) "%d %b" else "%b\n%Y") +
    scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, .05))) +
    labs(title = title, x = "Week of rash onset", y = "Weekly confirmed cases") +
    theme_sub()
}

for (st in names(STATE_LABEL)) {
  fig_no <- c(CDMX = 2, JAL = 3)[[st]]
  bands <- lapply(names(MODEL_LABEL), function(mdl) read_bands("main", st, mdl))
  names(bands) <- names(MODEL_LABEL)
  pA <- forecast_panel(bands$ENS, st, zoom = FALSE,
                       title = "Ensemble: full calibration window")
  pB <- lapply(names(MODEL_LABEL), function(mdl)
    forecast_panel(bands[[mdl]], st, zoom = TRUE, title = MODEL_LABEL[[mdl]]))
  p <- (pA / wrap_plots(pB, ncol = 2)) +
    plot_layout(heights = c(1, 2), guides = "collect") +
    plot_annotation(tag_levels = "A") &
    theme(legend.position = "bottom")
  save_sub(p, sprintf("Figure%d", fig_no), 15)
}

cat("Done.\n")
