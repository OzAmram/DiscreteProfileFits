# tools/

One-off scripts used to build, validate, and plot the studies under `studies/`.
These are **not** part of the core fitting framework (that stays in the repo
root: `Fitter.py`, `DataCardMaker.py`, `Utils.py`, `doFit.py`,
`fit_signalshapes.py`, `interpolation.py`, `makePostfitPlot.py`,
`run_fit_adaptive.py`).

**Run them from the repo root**, e.g. `python3 tools/bias_test.py`, not from
inside `tools/`. They call `doFit.py`/`interpolation.py` by relative name and
write results relative to the current directory, both of which assume the cwd is
the repo root. The few scripts that `import` core modules add the repo root to
`sys.path` themselves so the imports resolve.

## Bias / spurious-signal validation (→ studies/bias_validation)
- `bias_test.py` — single-fit signal-injection bias test.
- `ensemble_test.py` — background-only GoF ensemble (env `EN_*`).
- `ensemble_sb.py` — signal+background injection ensemble.
- `scan_spurious.py` — driver scanning n-sigma / bin-fraction configs (calls `ensemble_test.py`).
- `spurious_test.py`, `run_bkg_test.py`, `run_pull_test.py`, `scan_pull_config.py` — pull / spurious-signal drivers.
- `smooth_toy_test.py` — smooth-template toy test (method bias free of parent-MC fluctuations).
- `inject_signal.py`, `make_3pct.py` — build injected / 3%-subsample inputs.
- `compare_binsize_unc.py` — compares sigma(N) across bin fractions.
- `plot_bias_summary.py`, `plot_signal_pulls.py` — summary plots.
- `replot_ensemble.py`, `replot_bkg_ensemble.py` — refresh postfit chi2 / re-aggregate an existing ensemble without re-fitting.
- `run_adaptive_batch.sh`, `run_mock_batch.sh`, `run_mock_fixed.sh` — batch wrappers (mock inputs were removed).

## Old signal-shape pipeline (superseded by prep_signal_uaf.py / run_signal_fits_uaf.py)
- `extract_signal_masses.py` — extract dimuon mass from signal ROOT files to h5.
- `compare_signal_shapes.py` — overlay 2B/2G/Inv fitted DCB params vs mass (→ signal_comparison/).
- `draw_signals.py` — legacy Python-2 signal drawing script (does not run under python3; kept for reference).

## Misc
- `test.py`, `test_interpolation_closure.py` — ad-hoc checks.
