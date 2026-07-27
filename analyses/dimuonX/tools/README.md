# tools/

One-off scripts used to build, validate, and plot the studies referenced in
the note. These are **not** part of the core fitting framework (that lives
two directories up, at the repo root:
`Fitter.py`, `DataCardMaker.py`, `Utils.py`, `doFit.py`, `fit_signalshapes.py`,
`interpolation.py`, `makePostfitPlot.py`, `run_fit_adaptive.py`).

**Run them from `analyses/dimuonX/`** (the parent of this directory), e.g.
`python3 tools/bias_test.py`, not from inside `tools/`. They read/write data
(`signal_fits/`, `dimuonX_config.json`, background MC, output dirs, ...)
relative to that directory. The framework scripts they shell out to
(`doFit.py`, `run_fit_adaptive.py`, `interpolation.py`, `fit_signalshapes.py`)
are located via a `FRAMEWORK_DIR` computed from `__file__`, so those calls
work regardless of cwd — only the *data* paths need the right cwd.

## Bias / spurious-signal validation (note §"Validation of the Fit Procedure")
- `bias_test.py` — single-fit signal-injection bias test.
- `bias_test_toys.py` — the current bias test: combine-generated smooth-truth
  toys (fit the full background at high stats, generate `GenerateOnly` toys
  from each functional-form family, fit back through the full discrete-profiling
  envelope). Free of the parent-MC-fluctuation artifact that `ensemble_sb.py`-style
  subsample tests have. Produces the note's `fig:bias_toys`.
- `plot_bias_toys.py` — plots `bias_test_toys.py`'s output (median pull per
  generating family / injected signal / mass, with the |median pull| < 0.5 band).
- `ensemble_test.py` — background-only GoF ensemble over K 3%-subsamples (env `EN_*`).
- `ensemble_sb.py` — signal+background injection ensemble (yield pull + S+B GoF).
- `scan_spurious.py` — driver scanning n-sigma / bin-fraction configs (calls `ensemble_test.py`).
- `spurious_test.py`, `run_bkg_test.py`, `run_pull_test.py`, `scan_pull_config.py` — pull / spurious-signal drivers.
- `smooth_toy_test.py` — smooth-template toy test (method bias free of parent-MC fluctuations).
- `inject_signal.py`, `make_3pct.py` — build injected / 3%-subsample inputs.
- `compare_binsize_unc.py` — compares sigma(N) across bin fractions.
- `plot_bias_summary.py`, `plot_signal_pulls.py` — summary plots.
- `replot_ensemble.py`, `replot_bkg_ensemble.py` — refresh postfit chi2 / re-aggregate an existing ensemble without re-fitting.
- `run_adaptive_batch.sh`, `run_mock_batch.sh`, `run_mock_fixed.sh` — batch wrappers (mock inputs were removed).

## Note figures
- `make_note_figures.py` — builds the Fitting-section figures that come
  directly from `signal_fits/2B_loosemass/` and the GoF ensemble output
  (`interp_params_vs_mass.png`, `gof_ensemble.png`) into `note/Figures/8Systematics_figures/Fitting/`.
- `test_interpolation_closure.py` — leave-one-mass-out interpolation closure
  test; produces `note/Figures/8Systematics_figures/Fitting/closure_test_M{mass}.png`.

## Old signal-shape pipeline (superseded by prep_signal_loosemass.py / run_signal_fits_loosemass.py)
- `extract_signal_masses.py` — extract dimuon mass from signal ROOT files to h5.
- `compare_signal_shapes.py` — overlay 2B/2G/Inv fitted DCB params vs mass.
- `draw_signals.py` — legacy Python-2 signal drawing script (does not run under python3; kept for reference).

## Misc
- `test.py` — ad-hoc ROOT/TGraph scratch script, not a real test.
