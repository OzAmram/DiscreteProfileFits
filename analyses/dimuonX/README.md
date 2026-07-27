# DimuonX: opposite-sign di-muon resonance search

A model-independent bump hunt for a narrow resonance in the opposite-sign
di-muon mass spectrum, using the CATHODE anomaly-detection method (a
generative model trained on sidebands + a weakly-supervised classifier)
followed by the discrete-profiling background-plus-signal fit implemented by
the framework two directories up (`../../doFit.py` etc.).

This directory is a particular *application* of that framework: everything
here is specific to this search. See `../../README.md` for the framework
itself (fit mechanics, config format, background functional forms).

**Run scripts from here** (`analyses/dimuonX/`), not from `tools/` or the repo
root — data paths (`signal_fits/`, `dimuonX_config.json`, output dirs, ...)
are relative to this directory. Calls out to the framework scripts
(`doFit.py`, `run_fit_adaptive.py`, `interpolation.py`, `fit_signalshapes.py`)
are resolved relative to `__file__`, so those work regardless of cwd.

Scripts write their output locally (fit results, plots, intermediate data);
those outputs are not tracked in git (see `.gitignore`) since they're bulk and
fully regenerable. What *is* tracked — and all a fresh clone actually
contains — is the code below, the note, and the chosen benchmark's
interpolated signal shapes (`signal_fits/2B_loosemass/`).

## Signal shapes

The benchmark signal used for the model-independent scan is `TpTpTo2T2STo2Mu2B`
(the "2B" variant, T′→2T2S→2μ2b), fit and interpolated from the
`CATHODE_training_loosemass_v1` production (no hard 10/80 GeV mass cuts, so
the interpolation has real support at the edges of the search range).

- **`signal_fits/2B_loosemass/`** — the benchmark's shapes: `sig_fit_*.json`
  (DCB fit at each generated mass point) and `case_interpolation_M*.json`
  (interpolated template at every 0.1 GeV, 10–80 GeV). This is the only
  signal-shape directory tracked in git, and what `doFit.py`/
  `run_fit_adaptive.py` use by default throughout this directory.
- `signal_loosemass_manifest.json` — mass points / raw MC file paths per
  signal family in the `CATHODE_training_loosemass_v1` production, used by
  the scripts below to rebuild `signal_fits/2B_loosemass/` (or produce the
  same templates for a different family) from the raw MC on EOS.
- `prep_signal_loosemass.py` — extracts `dimu_mass` from the raw ROOT MC into
  per-mass, per-family HDF5.
- `run_signal_fits_loosemass.py` — fits a DCB per mass point per family
  (calls the framework's `fit_signalshapes.py`).
- `run_interpolation_loosemass.py` — interpolates each family to a 0.1 GeV
  grid (calls the framework's `interpolation.py`).
- `plot_signal_interp_loosemass.py`, `loo_interp_test_loosemass.py` — plot the
  interpolation and a leave-one-out closure test per family.
- `plot_signal_shape_overlay.py` — overlays fitted DCB parameters across
  signal families.
- `prep_signal_uaf.py`, `run_signal_fits_uaf.py` — the same prep/fit steps for
  an earlier, separate signal production (`UAF_training`), superseded by the
  loosemass production above for the current benchmark.

## Data scan automation

- `dimuonX_config.json` — copy of the framework's example config; this is the
  actual validated config used for every fit in this directory (func-forms,
  `sig_norm`, thresholds — see `../../README.md` for field docs).
- `samesign_dimuon_config.json` — a variant config for the same-sign scouting
  input (includes `bernPower`); not used by the current scan.
- `build_scan_plan.py` — builds the global mass-hypothesis grid from a
  collaborator's per-window CATHODE output directory (globs `<series>_<lo>_<hi>`
  window dirs on EOS, steps by local signal σ snapped to template nodes,
  assigns each hypothesis to its nearest window). Supports `--exclude-band` to
  drop hypotheses whose fit window overlaps a known-bad MC region.
- `run_data_scan.py` — drives `run_fit_adaptive.py` over every (window, mass)
  in a scan-plan JSON, in parallel.
- `summarize_data_scan.py` — rebuilds the summary JSON + 4 diagnostic plots
  (`scan_significance.png`, `scan_bkg_gof.png` — bkg-only *and* S+B goodness
  of fit overlaid, `scan_window_sizes.png`, `scan_func_forms.png`) from
  existing fit output, without re-fitting.

The scans referenced in the note (same-sign validation-region data, and the
CATHODE-selected signal-region background MC) are produced with these three
scripts, pointed at the relevant collaborator EOS directory.

## Signal-shape mismatch sensitivity study

Quantifies the loss in discovery sensitivity from using one signal template
(width class) to fit a different true signal shape — see
`note/appendix_shape_mismatch.tex`.

- `shape_mismatch_test.py` — injects real signal MC (narrow=VLL, medium=TP,
  wide=H2→H1H3) into the 3% anomaly-cut background at 4 masses, fits every
  injectant with every template (3×3 matrix), `--ntoys K` for a bootstrap
  ensemble with mean±SEM.
- `run_shape_mismatch.sh` — cmsenv wrapper.
- `plot_shape_mismatch.py` — aggregates the results into the heatmap and
  error-bar figures used in the appendix, plus a text summary.

## Bias / validation tooling

`tools/` holds the fit-procedure validation scripts referenced in the note's
"Validation of the Fit Procedure" section (combine-toy bias test, 3%-subsample
GoF/pull ensembles, spurious-signal tests) and the scripts that build the
note's figures. See `tools/README.md` for the full list and how to run them.

## The note

`note/` contains the analysis-note LaTeX sections and figures:
- `control_region.tex` — same-sign control-region validation (feature
  distributions, data/MC comparison, the CATHODE pipeline-validation scan).
- `fitting.tex` — the statistical procedure: signal shape modeling +
  interpolation, background discrete profiling, fit-procedure validation
  (3%-ensemble GoF, combine-toy bias test, and the post-selection
  background-MC scan).
- `appendix_shape_mismatch.tex` — the signal-shape mismatch study above.
- `Figures/` — one subdirectory per section (`6Control_Region_figures/`,
  `8Systematics_figures/`, `Appendix_ShapeMismatch_figures/`), populated by
  `tools/make_note_figures.py`, `tools/plot_bias_toys.py`,
  `tools/test_interpolation_closure.py`, and the scan/study scripts above.
