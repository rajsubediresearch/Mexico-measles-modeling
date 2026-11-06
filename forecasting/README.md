# Forecast benchmark: measles during the FIFA World Cup 2026 (Mexico City, Jalisco)

Phenomenological growth models (GLM, GRM, Richards) fitted with the
[QuantDiffForecast toolbox](https://github.com/rajsubediresearch/QuantDiffForecast-Toolbox)
to weekly confirmed measles cases up to the onset week containing the opening
match. Their forecasts give the incidence expected had the pre-tournament
trajectory continued, and observed counts are compared against them.

## Design

| | |
|---|---|
| Analysis release | DGE line list of 21 Sep 2026 (all evaluation weeks mature) |
| Series | Confirmed cases by state of notification and week of rash onset (Monday weeks) |
| Calibration | Mexico City 24 Nov 2025 – 8 Jun 2026 (29 wk); Jalisco 8 Sep 2025 – 8 Jun 2026 (40 wk) |
| Forecast window | Onset weeks 15 Jun – 3 Aug 2026 (8 wk) |
| Primary window | 15 Jun – 20 Jul: exposures during matches (11 Jun – 5 Jul) + 7–21 day exposure-to-rash interval |
| Extended window | 27 Jul – 3 Aug: one further generation (mean generation interval ≈ 12 days) |
| Models | GLM, GRM, RICH; equal-weight linear-pool ensemble (ENS) |
| Estimation | MLE, negative binomial NB1 (var = μ + αμ), 100 MultiStart points, 300 parametric bootstraps |
| Bounds on K | cumulative calibration cases to 10× that (`SC_Klow`: from 0.5×) |
| Outcomes | weekly exceedance of the 95% PI; cumulative excess over primary and full windows; WIS, MAE, coverage (secondary) |

### Analyses

| ID | Purpose |
|---|---|
| `main` | Calibration and evaluation on mature counts (21 Sep release) |
| `SA_realtime` | Calibration on the 13 Jul release exactly as in manuscript v8 (nowcast-adjusted); evaluated on mature counts. Shows how real-time reporting incompleteness changed the benchmark |
| `SB_end0601` | Calibration ends one week earlier (onset week 1–7 Jun), to remove onsets that could reflect pre-opening travel |
| `SC_Klow` | As `main`, with the lower bound of K relaxed to half the cumulative calibration cases (K hit its lower bound in the bootstrap fits) |

All settings live in the configuration block of `01_prepare_inputs.R`.
Option files are generated from it; edit the R config, not the `.m` files.

## Workflow

Run everything from this `forecasting/` folder.

1. **Nowcasting / weekly series** (project root): `run_pipeline_loop.R`
   produces `outputs/06_ensemble_input_<release>/ensemble_confirmed_<release>.xlsx`.
   Copy the release used here into `data/releases/`.
2. **`01_prepare_inputs.R`** (R) writes `inputs/`, `options/`, `manifest.csv`
   and `tables/data_revisions.csv`.
3. **`run_forecasts.m`** (MATLAB) runs every manifest row into
   `runs/<analysis>/<state>/<model>/`, then the ensembles into `.../ENS/`.
   Completed runs are skipped on re-run (delete `DONE` or set `FORCE = true`).
4. **`03_evaluate_excess.R`** (R) writes the evaluation tables to `tables/`.
5. **`04_rt.R`** (R): model-based Rt from the bootstrap trajectories and
   data-based Rt (Cori, 2-week window) from mature counts -> `tables/rt_*.csv`.
6. **`05_figures.R`** (R): all figures to `figures/` (PNG 300 dpi + PDF). All
   titles, labels, colours and periods are set at the top of the script.
   Figure 1 reads the cleaned line list from
   `../outputs/02_cleaned_data_<release>/` (not in git) so every state appears.
7. *Optional:* **`plot_toolbox_figures.m`** (MATLAB) redraws QuantDiffForecast's
   own diagnostic plots (fit and forecast, performance by horizon, parameter
   histograms) for every run from the saved results, without refitting, into
   `runs/<analysis>/<state>/<model>/toolbox_figures/` as PNG (300 dpi) and `.fig`.
   It needs the run's `.mat` files, which are kept locally but not in git.

## Toolbox

The toolbox is a git submodule pinned to a commit of the fork, so every
result can be traced to the exact code version (recorded in each
`run_info.txt`). After cloning this repository:

```bash
git submodule update --init
```

`run_forecasts.m` copies the toolbox code into `work/qdf/` (git-ignored) and
runs it from there, because the toolbox always reads `input/` and writes
`output/` next to its own files and most output file names omit the dataset.
The submodule itself is never modified.

## Folder map

```
data/releases/      weekly series per release (aggregated counts, all states)
data/v8_reference/  toolbox inputs used for manuscript v8 (GrowthPredict)
inputs/             generated toolbox inputs + *_dates.csv (t -> week, role, raw/nowcast)
options/            generated option files, one per analysis x state x model
manifest.csv        list of MATLAB runs
runs/               toolbox outputs, bootstrap curves, run_info.txt per run
tables/             evaluation tables
figures/            publication figures
toolbox/            QuantDiffForecast-Toolbox (submodule)
```
