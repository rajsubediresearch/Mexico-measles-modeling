# Project-specific toolbox additions

QuantDiffForecast (upstream, `toolbox/` submodule) ships the generalised logistic
model (`GLM.m`) but no longer the Richards-type growth models or the curve-based
linear-pool ensemble. These files, from the authors' earlier fork, are copied into
the toolbox working copy by `run_forecasts.m` and `plot_toolbox_figures.m`:

- `GRM.m` - generalised Richards model, dC/dt = r C^p (1 - (C/K)^a)
- `RICH.m` - Richards model, dC/dt = r C (1 - (C/K)^a)
- `getensemble_linearpool_from_curves.m` - exact empirical linear pool of bootstrap curves

They do not modify any upstream file.
