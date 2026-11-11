% run_forecasts.m  (step 2)
% Mexico measles — World Cup 2026 forecast benchmark (QuantDiffForecast)
%
% Runs every row of manifest.csv (written by 01_prepare_inputs.R) and stores
% each run's outputs in its own folder, runs/<analysis>/<state>/<model>/, so
% runs never overwrite each other. Then builds an equal-weight linear-pool
% ensemble of GLM, GRM and RICH for each analysis x state (runs/.../ENS/).
%
% The toolbox is used unmodified: its forecasting code (plus the project's
% toolbox_extras/: GRM, RICH, linear-pool ensemble) is copied to work/qdf/
% (git-ignored) because it always reads input/ and writes output/ next to
% itself. The copy is refreshed at the start of every session.
%
% Per run folder:
%   toolbox outputs (Forecast-*.csv, performance-*.csv, parameters-*.csv, ...)
%   curves_process.csv     T x B      bootstrap model trajectories (no noise)
%   curves_predictive.csv  T x 20B    trajectories with observation noise
%   run_result.mat         AICc and performance matrices
%   run_info.txt           seed, timings, toolbox commit, MATLAB version
%   DONE                   marker; delete it (or set FORCE) to rerun
%
% Usage: open this file in MATLAB and press Run, or type  run_forecasts
% in the Command Window from the forecasting/ folder.

clear; close all; clc
global method1 %#ok<GVMIS> used by the toolbox

% ══════════════════════════════════════════════════════════════════════════
% SETTINGS
FORCE          = false;   % true = rerun runs that already have a DONE marker
ONLY_ANALYSES  = {};      % e.g. {'main'}; empty = all analyses in manifest
ONLY_STATES    = {};      % e.g. {'CDMX'}; empty = all
BUILD_ENSEMBLE = true;
% ══════════════════════════════════════════════════════════════════════════

PROJ = fileparts(mfilename('fullpath'));
if isempty(PROJ), PROJ = pwd; end

QDF_ROOT = fullfile(PROJ, 'toolbox', 'QuantDiffForecast-Toolbox');
QDF_CODE = fullfile(QDF_ROOT, 'forecasting_odemodels code');
WORK     = fullfile(PROJ, 'work', 'qdf');

assert(isfolder(QDF_CODE), ['Toolbox not found at %s\n' ...
    'Run: git submodule update --init  (see README.md)'], QDF_CODE);
assert(isfile(fullfile(PROJ, 'manifest.csv')), ...
    'manifest.csv not found - run 01_prepare_inputs.R first');

% ── Fresh working copy of the toolbox code ────────────────────────────────
warning('off', 'MATLAB:rmpath:DirNotFound');
rmpath(QDF_CODE);                       % never use the submodule copy directly
if isfolder(WORK), rmpath(WORK); rmdir(WORK, 's'); end
mkdir(WORK); mkdir(fullfile(WORK, 'input')); mkdir(fullfile(WORK, 'output'));
copyfile(fullfile(QDF_CODE, '*.m'), WORK);
copyfile(fullfile(PROJ, 'toolbox_extras', '*.m'), WORK);   % GRM, RICH, linear-pool ensemble (not in upstream)
addpath(WORK, '-begin');

[~, sha] = system(sprintf('git -C "%s" rev-parse --short HEAD', QDF_ROOT));
sha = strtrim(sha);
fprintf('QuantDiffForecast commit: %s\n', sha);

% ── Manifest ──────────────────────────────────────────────────────────────
M = readtable(fullfile(PROJ, 'manifest.csv'), 'TextType', 'string', ...
              'Delimiter', ',');
keep = true(height(M), 1);
if ~isempty(ONLY_ANALYSES), keep = keep & ismember(M.analysis, ONLY_ANALYSES); end
if ~isempty(ONLY_STATES),   keep = keep & ismember(M.state,    ONLY_STATES);   end
M = M(keep, :);
fprintf('%d runs selected\n\n', height(M));

status = strings(height(M), 1);
tAll = tic;

for k = 1:height(M)
    r   = M(k, :);
    out = fullfile(PROJ, r.out_dir);

    if isfile(fullfile(out, 'DONE')) && ~FORCE
        fprintf('[%2d/%d] %-26s already done, skipped\n', k, height(M), r.run_id);
        status(k) = "skipped"; continue
    end
    fprintf('[%2d/%d] %-26s window=%d horizon=%d seed=%d ... ', ...
            k, height(M), r.run_id, r.windowsize, r.horizon, r.seed);

    if isfolder(out), rmdir(out, 's'); end
    mkdir(out);
    delete(fullfile(WORK, 'input',  '*'));
    delete(fullfile(WORK, 'output', '*'));
    copyfile(fullfile(PROJ, r.input_file), fullfile(WORK, 'input'));

    optDir = fullfile(PROJ, r.options_dir);
    addpath(optDir);
    rng(r.seed, 'twister');
    t0 = tic;
    try
        [AICcs, performanceC, performanceF, forecast_model1, forecast_model12, ...
            data1, datalatest] = Run_Forecasting_ODEModel( ...
            str2func(r.options_fn), 1, 1, r.windowsize, r.horizon);
        close all

        if ~isempty(dir(fullfile(WORK, 'output', '*')))
            movefile(fullfile(WORK, 'output', '*'), out);
        end
        writematrix(forecast_model1,  fullfile(out, 'curves_process.csv'));
        writematrix(forecast_model12, fullfile(out, 'curves_predictive.csv'));
        save(fullfile(out, 'run_result.mat'), 'AICcs', 'performanceC', ...
             'performanceF', 'data1', 'datalatest', '-v7');

        elapsed = toc(t0);
        fid = fopen(fullfile(out, 'run_info.txt'), 'w');
        fprintf(fid, 'run_id: %s\nanalysis: %s\nstate: %s\nmodel: %s\n', ...
                r.run_id, r.analysis, r.state, r.model);
        fprintf(fid, 'input: %s\noptions: %s\nwindowsize: %d\nhorizon: %d\n', ...
                r.input_file, r.options_fn, r.windowsize, r.horizon);
        fprintf(fid, 'seed: %d\nelapsed_min: %.1f\nfinished: %s\n', ...
                r.seed, elapsed/60, string(datetime('now')));
        fprintf(fid, 'toolbox_commit: %s\nmatlab: %s\n', sha, version);
        fclose(fid);
        fclose(fopen(fullfile(out, 'DONE'), 'w'));
        % keep a copy of the exact options file next to its outputs
        copyfile(fullfile(optDir, r.options_fn + ".m"), out);

        fprintf('done (%.1f min)\n', elapsed/60);
        status(k) = "ok";
    catch ME
        close all
        fid = fopen(fullfile(out, 'ERROR.txt'), 'w');
        fprintf(fid, '%s', getReport(ME, 'extended', 'hyperlinks', 'off'));
        fclose(fid);
        fprintf('FAILED: %s\n', ME.message);
        status(k) = "FAILED: " + ME.message;
    end
    rmpath(optDir);
end

% ── Equal-weight ensemble (linear pool) per analysis x state ──────────────
if BUILD_ENSEMBLE
    fprintf('\nBuilding ensembles...\n');
    groups = unique(M(:, {'analysis', 'state'}));
    models = ["GLM", "GRM", "RICH"];
    for g = 1:height(groups)
        an = groups.analysis(g); st = groups.state(g);
        base = fullfile(PROJ, 'runs', an, st);
        ok = all(arrayfun(@(m) isfile(fullfile(base, m, 'DONE')), models));
        if ~ok
            fprintf('  %s/%s: not all models done, ensemble skipped\n', an, st);
            continue
        end
        ensDir = fullfile(base, 'ENS');
        if ~isfolder(ensDir), mkdir(ensDir); end
        w = ones(1, numel(models)) / numel(models);
        for kind = ["predictive", "process"]
            A = arrayfun(@(m) readmatrix(fullfile(base, m, ...
                    "curves_" + kind + ".csv")), models, 'UniformOutput', false);
            [~, f1, ~, ~, curvesEns] = getensemble_linearpool_from_curves( ...
                w, A{:}, 'Seed', 1);
            writematrix(curvesEns, fullfile(ensDir, "curves_" + kind + ".csv"));
            writematrix(f1, fullfile(ensDir, "quantiles_" + kind + "_median_LB_UB.csv"));
        end
        fid = fopen(fullfile(ensDir, 'run_info.txt'), 'w');
        fprintf(fid, 'Equal-weight linear pool of %s\nfinished: %s\ntoolbox_commit: %s\n', ...
                strjoin(models, ', '), string(datetime('now')), sha);
        fclose(fid);
        fclose(fopen(fullfile(ensDir, 'DONE'), 'w'));
        fprintf('  %s/%s: ensemble done\n', an, st);
    end
end

fprintf('\n════════════════════════════════════════\n');
fprintf('Finished in %.1f min\n', toc(tAll)/60);
for k = 1:height(M), fprintf('  %-26s %s\n', M.run_id(k), status(k)); end
fprintf('════════════════════════════════════════\n');
