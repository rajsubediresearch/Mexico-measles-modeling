% plot_toolbox_figures.m  (optional, after run_forecasts.m)
% Mexico measles — the toolbox's own diagnostic figures for every run
%
% Re-draws QuantDiffForecast's standard plots (fit + forecast with bootstrap
% curves, forecast performance by horizon, parameter histograms, ...) from the
% saved results - no refitting - and saves each figure as
%   runs/<analysis>/<state>/<model>/toolbox_figures/<run>_figNNN_<type>.png  (300 dpi)
%   runs/<analysis>/<state>/<model>/toolbox_figures/<run>_figNNN_<type>.fig  (editable)
%
% Publication figures are made separately by 05_figures.R.
%
% Usage: open in MATLAB from the forecasting/ folder and press Run.

clear; close all; clc
global method1 %#ok<GVMIS> used by the toolbox

% ══════════════════════════════════════════════════════════════════════════
% SETTINGS
FORCE         = false;      % true = redraw runs that already have figures
ONLY_ANALYSES = {};        % e.g. {'main'}; empty = all
ONLY_STATES   = {};        % e.g. {'CDMX'}
SAVE_FIG      = true;      % also save MATLAB .fig files
FIG_POSITION  = [50 50 1600 1000];   % window size in pixels while drawing
FONT_SIZE     = 11;        % fixed font size (pt) for saved figures; the toolbox
                           % scales fonts with screen DPI, so they are reset here
TITLE_SIZE    = 10;        % axis titles (histogram titles carry long CI text)
% ══════════════════════════════════════════════════════════════════════════

PROJ = fileparts(mfilename('fullpath'));
if isempty(PROJ), PROJ = pwd; end
QDF_CODE = fullfile(PROJ, 'toolbox', 'QuantDiffForecast-Toolbox', 'forecasting_odemodels code');
WORK     = fullfile(PROJ, 'work', 'qdf');
assert(isfolder(QDF_CODE), 'Toolbox not found: run  git submodule update --init');

% fresh working copy of the toolbox (same approach as run_forecasts.m)
warning('off', 'MATLAB:rmpath:DirNotFound');
rmpath(QDF_CODE);
if isfolder(WORK), rmpath(WORK); rmdir(WORK, 's'); end
mkdir(WORK); mkdir(fullfile(WORK, 'input')); mkdir(fullfile(WORK, 'output'));
copyfile(fullfile(QDF_CODE, '*.m'), WORK);
copyfile(fullfile(PROJ, 'toolbox_extras', '*.m'), WORK);   % GRM, RICH, linear-pool ensemble (not in upstream)
addpath(WORK, '-begin');

M = readtable(fullfile(PROJ, 'manifest.csv'), 'TextType', 'string', 'Delimiter', ',');
keep = true(height(M), 1);
if ~isempty(ONLY_ANALYSES), keep = keep & ismember(M.analysis, ONLY_ANALYSES); end
if ~isempty(ONLY_STATES),   keep = keep & ismember(M.state,    ONLY_STATES);   end
M = M(keep, :);

oldPos = get(groot, 'defaultFigurePosition');
set(groot, 'defaultFigurePosition', FIG_POSITION);
cleanupObj = onCleanup(@() set(groot, 'defaultFigurePosition', oldPos));

for k = 1:height(M)
    r      = M(k, :);
    runDir = fullfile(PROJ, r.out_dir);
    figDir = fullfile(runDir, 'toolbox_figures');

    if ~isfile(fullfile(runDir, 'DONE'))
        fprintf('[%2d/%d] %-26s not run yet, skipped\n', k, height(M), r.run_id); continue
    end
    if isfolder(figDir) && ~isempty(dir(fullfile(figDir, '*.png'))) && ~FORCE
        fprintf('[%2d/%d] %-26s figures exist, skipped\n', k, height(M), r.run_id); continue
    end
    if isempty(dir(fullfile(runDir, '*.mat')))
        fprintf('[%2d/%d] %-26s no .mat results (needed to redraw), skipped\n', k, height(M), r.run_id); continue
    end
    fprintf('[%2d/%d] %-26s ... ', k, height(M), r.run_id);

    % put this run's inputs and saved results where the toolbox expects them
    delete(fullfile(WORK, 'input',  '*'));
    delete(fullfile(WORK, 'output', '*'));
    copyfile(fullfile(PROJ, r.input_file), fullfile(WORK, 'input'));
    copyfile(fullfile(runDir, '*.mat'), fullfile(WORK, 'output'));
    copyfile(fullfile(runDir, '*.csv'), fullfile(WORK, 'output'));

    optDir = fullfile(PROJ, r.options_dir);
    addpath(optDir);
    try
        close all
        plotForecast_ODEModel(@() options_with_plots(char(r.options_fn)), ...
                              1, 1, r.windowsize, r.horizon);
        if ~isfolder(figDir), mkdir(figDir); end
        figs = findobj('Type', 'figure');
        [~, order] = sort([figs.Number]); figs = figs(order);
        for f = figs'
            normalise_figure(f, fig_type(f.Number), FONT_SIZE, TITLE_SIZE);
            base = sprintf('%s_fig%03d_%s', r.run_id, f.Number, fig_type(f.Number));
            exportgraphics(f, fullfile(figDir, [base '.png']), 'Resolution', 300);
            if SAVE_FIG, savefig(f, fullfile(figDir, [base '.fig'])); end
        end
        fprintf('%d figures saved\n', numel(figs));
    catch ME
        fprintf('FAILED: %s\n', ME.message);
    end
    close all
    rmpath(optDir);
end
fprintf('\nDone. Figures are in runs/<analysis>/<state>/<model>/toolbox_figures/\n');

% ── local functions ───────────────────────────────────────────────────────
function varargout = options_with_plots(fn)
% Calls a generated options file and switches on the toolbox's plotting
% (printscreen1, the 15th output), leaving the options file unchanged.
[varargout{1:15}] = feval(fn);
varargout{15} = 1;
end

function normalise_figure(f, type, fs, ts)
% Fixed physical size and font sizes, independent of the screen the figure
% was drawn on. Width scales with the number of panels side by side.
ax = findobj(f, 'Type', 'axes');
nPanels = max(1, numel(ax));
switch type
    case 'fit_forecast',          wh = [10 5];
    case 'forecast_performance',  wh = [10 7];
    case {'parameter_histograms', 'diagnostics'}
                                  wh = [3.6 * nPanels, 3.6];
    case 'state_variables',       wh = [8 5];
    otherwise,                    wh = [10 6];
end
set(f, 'Units', 'inches');
f.Position = [1 1 wh];
set(findall(f, '-property', 'FontSize'), 'FontSize', fs);
for a = ax'
    a.Title.FontSize = ts;
    a.Title.FontWeight = 'normal';
end
drawnow
end

function t = fig_type(n)
% Figure numbers used by plotForecast_ODEModel (100+, 200+, ... per window)
if n == 400,      t = 'state_variables';
elseif n < 200,   t = 'fit_forecast';
elseif n < 300,   t = 'forecast_performance';
elseif n < 400,   t = 'parameter_histograms';
elseif n < 500,   t = 'diagnostics';
else,             t = 'other';
end
end
