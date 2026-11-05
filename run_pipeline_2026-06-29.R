# run_pipeline_2026-06-29.R
# Mexico DGE Measles — Full Weekly Pipeline
# Author: Raj (GRA — Dr. Chowell Lab, Georgia State University)
#
# INSTRUCTIONS FOR NEXT WEEK:
#   1. Save this file with the new date tag (e.g. run_pipeline_2026-07-06.R)
#   2. Update DATA_FILE to the new Excel filename
#   3. Update DATE_TAG to the new date
#   4. Run the whole script — all outputs saved to dated subfolders automatically
#
# Output structure:
#   outputs/
#     02_cleaned_data_2026-06-29/
#     05_nowcast_ensemble_2026-06-29/
#     06_ensemble_input_2026-06-29/
#     07_toolbox_input_2026-06-29/
#     08_panel_plots_2026-06-29/

# ══════════════════════════════════════════════════════════════════════════════
# UPDATE THESE TWO LINES EACH WEEK
DATA_FILE <- "Base_Lunes_06_07_2026_ANONIMIZADA_GC.xlsx"
DATE_TAG  <- "2026-07-06"
# ══════════════════════════════════════════════════════════════════════════════

setwd("D:/PhD coursework/GRA/Spring 2026/Mexico_Measles/nowcasting")

library(readxl)
library(dplyr)
library(lubridate)
library(tidyr)
library(NobBS)
library(ggplot2)
library(writexl)

cat("══════════════════════════════════════════\n")
cat("Mexico Measles Pipeline —", DATE_TAG, "\n")
cat("Data file:", DATA_FILE, "\n")
cat("══════════════════════════════════════════\n\n")

# ── PARAMETERS ────────────────────────────────────────────────────────────────

MAX_D  <- 30
WIN    <- 60
MIN_N  <- 20
NADAPT <- 3000

series_defs <- list(
  confirmed = "CONFIRMADO",
  suspected = "EN ESTUDIO",
  total     = c("CONFIRMADO", "EN ESTUDIO")
)

# ── CREATE OUTPUT FOLDERS ─────────────────────────────────────────────────────

dirs <- list(
  cleaned  = paste0("outputs/02_cleaned_data_",       DATE_TAG),
  nowcast  = paste0("outputs/05_nowcast_ensemble_",   DATE_TAG),
  ensemble = paste0("outputs/06_ensemble_input_",     DATE_TAG),
  toolbox  = paste0("outputs/07_toolbox_input_",      DATE_TAG),
  plots    = paste0("outputs/08_panel_plots_",        DATE_TAG)
)
for (d in c(unlist(dirs),
            file.path(dirs$nowcast, "nowcast_by_state"),
            file.path(dirs$toolbox, c("confirmed","suspected","total")))) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)
}

# ══════════════════════════════════════════════════════════════════════════════
# STEP 1 — LOAD & CLEAN
# ══════════════════════════════════════════════════════════════════════════════

cat("[1/5] Loading and cleaning data...\n")

parse_serial <- function(x) as.Date(as.numeric(x), origin = "1899-12-30")

df_raw <- read_excel(DATA_FILE, sheet = 1, col_types = "text")

df <- df_raw %>%
  mutate(
    FEC_INI_EXANT   = parse_serial(FEC_INI_EXANT),
    FEC_PRI_CONSULT = parse_serial(FEC_PRI_CONSULT),
    delay_days      = as.numeric(FEC_PRI_CONSULT - FEC_INI_EXANT),
    IDE_EDA_ANO     = as.numeric(IDE_EDA_ANO),
    AÑO             = as.integer(AÑO),
    SEMANA          = as.integer(SEMANA),
    age_group       = cut(IDE_EDA_ANO,
                          breaks = c(0,5,10,15,20,30,40,50,Inf),
                          labels = c("0-4","5-9","10-14","15-19",
                                     "20-29","30-39","40-49","50+"),
                          right  = FALSE),
    state           = trimws(DES_EDO),
    state_notif     = gsub("\\s+", " ", trimws(DES_EDO_NOTIFICANTE))
  )

df_all       <- df
df_confirmed <- df %>% filter(CLAS_FINAL_SARAMPION == "CONFIRMADO")

saveRDS(df_all,       file.path(dirs$cleaned, "measles_all.rds"))
saveRDS(df_confirmed, file.path(dirs$cleaned, "measles_confirmed.rds"))

cat(sprintf("  Rows: %s | Confirmed: %s | EN ESTUDIO: %s | Discarded: %s\n",
            nrow(df),
            sum(df$CLAS_FINAL_SARAMPION == "CONFIRMADO"),
            sum(df$CLAS_FINAL_SARAMPION == "EN ESTUDIO"),
            sum(df$CLAS_FINAL_SARAMPION == "DESCARTADO")))
cat(sprintf("  Onset range: %s to %s\n",
            min(df$FEC_INI_EXANT, na.rm=TRUE),
            max(df$FEC_INI_EXANT, na.rm=TRUE)))

# ══════════════════════════════════════════════════════════════════════════════
# STEP 2 — NOWCASTING (NobBS)
# ══════════════════════════════════════════════════════════════════════════════

cat("\n[2/5] Running nowcasts...\n")

df_nodiscard <- df %>% filter(CLAS_FINAL_SARAMPION != "DESCARTADO")
GLOBAL_NOW   <- max(df_nodiscard$FEC_INI_EXANT, na.rm = TRUE)
cat(sprintf("  GLOBAL_NOW: %s\n", GLOBAL_NOW))

states <- c("National", sort(unique(df_nodiscard$state_notif)))

prep_nobbs <- function(data, now) {
  data %>%
    filter(!is.na(FEC_INI_EXANT), !is.na(FEC_PRI_CONSULT)) %>%
    transmute(onset_date  = as.Date(FEC_INI_EXANT),
              report_date = as.Date(FEC_PRI_CONSULT)) %>%
    filter(report_date <= as.Date(now)) %>%
    as.data.frame() %>%
    arrange(onset_date)
}

run_nobbs <- function(linelist, now, label) {
  cat(sprintf("  %-25s n=%-6s now=%s\n", label, nrow(linelist), now))
  result <- NobBS(
    data          = linelist,
    now           = as.Date(now),
    units         = "1 day",
    onset_date    = "onset_date",
    report_date   = "report_date",
    moving_window = WIN,
    max_D         = MAX_D,
    specs         = list(dist="NB", nBurnin=2000, nSamp=10000,
                         nThin=1, nAdapt=NADAPT, conf=0.95)
  )
  result$estimates %>%
    mutate(geography = label, now_used = as.Date(now))
}

aggregate_weekly <- function(daily_est) {
  daily_est %>%
    mutate(week_start = floor_date(onset_date, "week", week_start=1)) %>%
    group_by(geography, week_start, now_used) %>%
    summarise(estimate=sum(estimate,na.rm=TRUE),
              lower=sum(lower,na.rm=TRUE),
              upper=sum(upper,na.rm=TRUE),
              .groups="drop") %>%
    arrange(geography, week_start)
}

nowcast_results <- list(confirmed=list(), suspected=list(), total=list())
skipped <- list()

for (series_name in names(series_defs)) {
  cat(sprintf("\n  — %s —\n", series_name))
  cf <- series_defs[[series_name]]

  for (geo in states) {
    df_g <- if (geo=="National") df_nodiscard else
      df_nodiscard %>% filter(state_notif==geo)
    df_g <- df_g %>% filter(CLAS_FINAL_SARAMPION %in% cf)
    ll   <- prep_nobbs(df_g, GLOBAL_NOW)

    if (nrow(ll) < MIN_N) {
      skipped[[length(skipped)+1]] <- tibble(
        series=series_name, geography=geo, n=nrow(ll), reason=paste0("n < ",MIN_N))
      next
    }

    tryCatch({
      est  <- run_nobbs(ll, GLOBAL_NOW, geo)
      week <- aggregate_weekly(est)
      nowcast_results[[series_name]][[geo]] <- week
    }, error = function(e) {
      cat(sprintf("  ERROR %s: %s\n", geo, e$message))
      skipped[[length(skipped)+1]] <<- tibble(
        series=series_name, geography=geo, n=nrow(ll), reason=e$message)
    })
  }
}

# Save nowcast Excel
skipped_df <- if (length(skipped)>0) bind_rows(skipped) else
  tibble(series=character(),geography=character(),n=integer(),reason=character())

write_xlsx(
  c(lapply(names(nowcast_results), function(s) bind_rows(nowcast_results[[s]])) %>%
      setNames(names(nowcast_results)),
    list(skipped=skipped_df)),
  file.path(dirs$nowcast, paste0("nowcast_ensemble_", DATE_TAG, ".xlsx"))
)
cat(sprintf("\n  Saved nowcast Excel. Skipped: %s\n", nrow(skipped_df)))

# ══════════════════════════════════════════════════════════════════════════════
# STEP 3 — MERGE OBSERVED + NOWCASTED
# ══════════════════════════════════════════════════════════════════════════════

cat("\n[3/5] Merging observed + nowcasted...\n")

make_observed <- function(data, class_filter, geo) {
  df_g <- if (geo=="National") data else data %>% filter(state_notif==geo)
  df_g <- df_g %>% filter(CLAS_FINAL_SARAMPION %in% class_filter)
  if (nrow(df_g)==0) return(NULL)
  start    <- floor_date(min(df_g$FEC_INI_EXANT,na.rm=TRUE),"week",week_start=1)
  end      <- floor_date(max(df_g$FEC_INI_EXANT,na.rm=TRUE),"week",week_start=1)
  week_seq <- seq(start, end, by="week")
  df_g %>%
    filter(!is.na(FEC_INI_EXANT)) %>%
    mutate(week_start=floor_date(FEC_INI_EXANT,"week",week_start=1)) %>%
    count(week_start, name="observed") %>%
    right_join(tibble(week_start=week_seq), by="week_start") %>%
    mutate(observed=replace_na(observed,0L), geography=geo) %>%
    arrange(week_start) %>%
    select(geography, week_start, observed)
}

merge_series <- function(observed_df, nowcast_df) {
  observed_df %>%
    left_join(nowcast_df %>%
                select(geography, week_start,
                       nowcast_estimate=estimate,
                       nowcast_lower=lower,
                       nowcast_upper=upper),
              by=c("geography","week_start")) %>%
    mutate(
      nowcast_estimate = ifelse(!is.na(nowcast_estimate) &
                                  nowcast_estimate>observed,
                                nowcast_estimate, NA_real_),
      nowcast_lower    = ifelse(!is.na(nowcast_estimate), nowcast_lower, NA_real_),
      nowcast_upper    = ifelse(!is.na(nowcast_estimate), nowcast_upper, NA_real_)
    )
}

ensemble_sheets <- list()
for (series_name in names(series_defs)) {
  now_df   <- bind_rows(nowcast_results[[series_name]])
  cf       <- series_defs[[series_name]]
  geos     <- sort(unique(now_df$geography))
  all_data <- list()
  for (geo in geos) {
    obs <- make_observed(df_nodiscard, cf, geo)
    if (is.null(obs)) next
    all_data[[geo]] <- merge_series(obs, now_df %>% filter(geography==geo))
  }
  combined <- bind_rows(all_data)
  sheet_list <- c(list(all_geographies=combined), all_data)
  write_xlsx(sheet_list,
             file.path(dirs$ensemble,
                       paste0("ensemble_",series_name,"_",DATE_TAG,".xlsx")))
  ensemble_sheets[[series_name]] <- combined
  cat(sprintf("  Saved: ensemble_%s_%s.xlsx (%s geographies)\n",
              series_name, DATE_TAG, length(all_data)))
}

# ══════════════════════════════════════════════════════════════════════════════
# STEP 4 — TOOLBOX TXT FILES
# ══════════════════════════════════════════════════════════════════════════════

cat("\n[4/5] Exporting toolbox txt files...\n")

for (series_name in names(ensemble_sheets)) {
  df_s       <- ensemble_sheets[[series_name]]
  out_folder <- file.path(dirs$toolbox, series_name)
  geos       <- sort(unique(df_s$geography))

  for (geo in geos) {
    df_g <- df_s %>% filter(geography==geo) %>% arrange(week_start)
    start_dt <- as.character(min(df_g$week_start))
    end_dt   <- as.character(max(df_g$week_start))
    df_g <- df_g %>%
      mutate(cases = ifelse(!is.na(nowcast_estimate), nowcast_estimate, observed),
             t     = seq(0, nrow(.)-1))
    geo_clean <- gsub(" ","_", geo)
    fname     <- paste0(geo_clean,"_",start_dt,"_",end_dt,".txt")
    write.table(df_g %>% select(t, cases),
                file.path(out_folder, fname),
                sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)
  }
  cat(sprintf("  %s: %s files saved\n", series_name, length(geos)))
}

# ══════════════════════════════════════════════════════════════════════════════
# STEP 5 — PANEL PLOTS
# ══════════════════════════════════════════════════════════════════════════════

cat("\n[5/5] Generating panel plots...\n")

load_series <- function(series_name) {
  folder <- file.path(dirs$toolbox, series_name)
  files  <- list.files(folder, pattern="\\.txt$", full.names=FALSE)
  lapply(files, function(fn) {
    dates    <- regmatches(fn, gregexpr("\\d{4}-\\d{2}-\\d{2}", fn))[[1]]
    start_dt <- as.Date(dates[1])
    geo      <- gsub("_\\d{4}-\\d{2}-\\d{2}_\\d{4}-\\d{2}-\\d{2}\\.txt$","",fn)
    d <- read.table(file.path(folder,fn), sep="\t", header=FALSE,
                    col.names=c("t","cases"))
    d$geography <- geo
    d$week      <- start_dt + d$t*7
    d
  }) %>% bind_rows()
}

make_panel <- function(df, series_name) {
  geo_order <- df %>%
    group_by(geography) %>%
    summarise(total=sum(cases,na.rm=TRUE)) %>%
    arrange(desc(total)) %>% pull(geography)
  df$geography <- factor(df$geography, levels=geo_order)

  ggplot(df, aes(x=week, y=cases)) +
    geom_col(fill="#2c7bb6", width=5) +
    facet_wrap(~geography, scales="free_y", ncol=4) +
    scale_x_date(date_labels="%b\n%Y", date_breaks="3 months") +
    labs(title    = paste0("Mexico Measles — ",toupper(series_name)," cases (nowcasted)"),
         subtitle = paste0("Weekly incidence by state of notification | Data: ",DATE_TAG),
         x=NULL, y="Cases",
         caption="Last week nowcasted where estimate > observed | NobBS v1.1.1") +
    theme_minimal(base_size=15) +
    theme(strip.text       = element_text(size=20, face="bold"),
          strip.background = element_rect(fill="grey93", color=NA),
          axis.text.x      = element_text(size=15, hjust=0.5),
          axis.text.y      = element_text(size=15),
          axis.title.y     = element_text(size=15),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_line(color="grey90"),
          panel.spacing    = unit(1,"lines"),
          plot.title       = element_text(size=20, face="bold"),
          plot.subtitle    = element_text(size=18, color="grey40"),
          plot.caption     = element_text(size=16, color="grey50"),
          plot.margin      = margin(10,15,10,10))
}

for (s in names(series_defs)) {
  df_p      <- load_series(s) %>% filter(geography != "National")
  n_states  <- length(unique(df_p$geography))
  p         <- make_panel(df_p, s)
  out_file  <- file.path(dirs$plots, paste0("panel_",s,".png"))
  ggsave(out_file, p, width=16, height=ceiling(n_states/4)*3.2+2, dpi=200)
  cat(sprintf("  Saved: panel_%s.png (%s states)\n", s, n_states))
}

# ══════════════════════════════════════════════════════════════════════════════
# DONE
# ══════════════════════════════════════════════════════════════════════════════

cat("\n══════════════════════════════════════════\n")
cat("Pipeline complete:", DATE_TAG, "\n")
cat("Outputs:\n")
for (nm in names(dirs)) cat(sprintf("  %s\n", normalizePath(dirs[[nm]])))
cat("══════════════════════════════════════════\n")
