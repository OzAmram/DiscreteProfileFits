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

Most of what's listed below is **not tracked in git** (see `.gitignore`) — it's
bulk per-fit output or raw MC, fully regenerable from the tracked code. The
exceptions are the code itself, the note, and the chosen benchmark's
interpolated signal shapes (`signal_fits/2B_loosemass/*.json`), which are
tracked so the analysis is reproducible without re-deriving them from EOS.

## Signal shapes

The benchmark signal used for the model-independent scan is `TpTpTo2T2STo2Mu2B`
(the "2B" variant, T′→2T2S→2μ2b), fit and interpolated from the
`CATHODE_training_loosemass_v1` production (no hard 10/80 GeV mass cuts, so
the interpolation has real support at the edges of the search range).

- `signal_loosemass_manifest.json` — mass points / raw MC file paths per
  signal family in that production.
- `prep_signal_loosemass.py` — extracts `dimu_mass` from the raw ROOT MC into
  per-mass, per-family HDF5 (→ `signal_data_loosemass/<group>/mass_<M>.h5`,
  gitignored).
- `run_signal_fits_loosemass.py` — fits a DCB per mass point per family
  (→ `signal_fits_loosemass/<group>/sig_fit_<M>.json`, gitignored; calls the
  framework's `fit_signalshapes.py`).
- `run_interpolation_loosemass.py` — interpolates each family to a 0.1 GeV
  grid (→ `signal_fits_loosemass/<group>/case_interpolation_M<M>.json`,
  gitignored; calls the framework's `interpolation.py`).
- `plot_signal_interp_loosemass.py`, `loo_interp_test_loosemass.py` — plot the
  interpolation and a leave-one-out closure test per family
  (→ `signal_comparison_loosemass/`, `signal_loo_interp/`, both gitignored).
- **`signal_fits/2B_loosemass/`** — the chosen benchmark's shapes only
  (`sig_fit_*.json` + `case_interpolation_M*.json`, 10–80 GeV), copied out of
  `signal_fits_loosemass/TpTpTo2T2STo2Mu2B_MTp1000/`. **This is the only
  signal-shape directory tracked in git**; it's what `doFit.py`/
  `run_fit_adaptive.py` use by default throughout this directory.
- `signal_fits/2B/`, `signal_fits/2G/`, `signal_fits/Inv/`, `signal_data/` —
  an earlier signal production (hard 10/80 GeV mass cuts; 2G/Inv are the
  photon/invisible T′ variants, unused by the current default). Gitignored,
  superseded by the loosemass production above; kept locally for comparison.
- `signal_data_uaf/`, `signal_fits_uaf/`, `signal_uaf_manifest.json`,
  `prep_signal_uaf.py`, `run_signal_fits_uaf.py` — a separate, earlier signal
  production (`UAF_training`, 10 families × 19 masses). Gitignored, superseded
  by the loosemass production for the current default.
- `plot_signal_shape_overlay.py` — overlays fitted DCB parameters across
  signal families/productions.

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

### Scan outputs (all gitignored, regenerate via the scripts above)
- `sr_bkgmc_v3_fits/` — full-range scan (97 pts, 13.8–75 GeV) on the actual
  CATHODE-selected signal-region background MC (`Results_sr_dimuon_combined_v3`);
  validates the pipeline post-anomaly-selection. `[47.5, 52.5]` GeV excluded
  (known deficiency in that MC sample).
- `vr_v4_fits/` — full-range scan (119 pts) on the same-sign validation-region
  data (`Results_vr_dimuon_combined_v4`), the current `note/control_region.tex`
  results.
- `data_fits/`, `data_fits_separate_v2/` — **stale**, pre-grid-rebuild v1/v2
  same-sign scans (superseded by `vr_v4_fits/`; large, ~1.4–1.9 GB each).
- `scan_plan.json`, `data_scan_summary.json`, `data_scan_run.log` — outputs
  from the stale v1 scan above.

## Signal-shape mismatch sensitivity study

Quantifies the loss in discovery sensitivity from using one signal template
(width class) to fit a different true signal shape — see
`note/appendix_shape_mismatch.tex`.

- `shape_mismatch_test.py` — injects real signal MC (narrow=VLL, medium=TP,
  wide=H2→H1H3) into the 3% anomaly-cut background at 4 masses, fits every
  injectant with every template (3×3 matrix), `--ntoys K` for a bootstrap
  ensemble with mean±SEM.
- `run_shape_mismatch.sh` — cmsenv wrapper.
- `plot_shape_mismatch.py` — aggregates → `shape_mismatch_test/shape_mismatch_{Z,rbias}.png`
  (heatmaps) and `_errbars.png` variants, plus a text summary. Output dir and
  `*.log` files are gitignored.

## Bias / validation tooling

See `tools/README.md` for the full list. Key outputs referenced by the note
(`fitting.tex`, §"Validation of the Fit Procedure"): `bias_toys/`,
`bias_toys_ftestfix/` (combine-toy bias test), `ensemble_capped_adaptive/`,
`ensemble_ftestfix/`, `ensemble_sb_capped_adaptive/` (3%-subsample GoF/pull
ensembles). All gitignored, regenerable via `tools/bias_test_toys.py` /
`tools/ensemble_test.py` / `tools/ensemble_sb.py`.

## Raw data inputs (gitignored)

- `bkg_mc_masses.h5` — full background MC (pre-anomaly-selection), used as
  the "truth" input for the smooth-toy bias test.
- `bkg_mc_3pct.h5` — a fixed 3% subsample of the above, emulating the
  anomaly-cut statistics for validation.
- `scouting_2024_same_sign_dimuon_mass.h5` — raw 2024 scouting same-sign data.

## The note

`note/` contains the analysis-note LaTeX sections and figures:
- `control_region.tex` — same-sign control-region validation (feature
  distributions, data/MC comparison, the `vr_v4_fits/` pipeline-validation
  scan).
- `fitting.tex` — the statistical procedure: signal shape modeling +
  interpolation, background discrete profiling, fit-procedure validation
  (3%-ensemble GoF, combine-toy bias test, and the `sr_bkgmc_v3_fits/`
  post-selection background-MC scan).
- `appendix_shape_mismatch.tex` — the signal-shape mismatch study above.
- `Figures/` — one subdirectory per section (`6Control_Region_figures/`,
  `8Systematics_figures/`, `Appendix_ShapeMismatch_figures/`), populated by
  `tools/make_note_figures.py`, `tools/plot_bias_toys.py`,
  `tools/test_interpolation_closure.py`, and the scan/study scripts above.

`fitting.tex` also exists as a stray duplicate directly in this directory
(pre-dates `note/` being organized, superseded by `note/fitting.tex`, gitignored)
— safe to delete.

## Other

- `studies/` — archived one-off studies (gitignored).
- `datacardInputs_mass_test.root` — a stray scratch artifact from a manual
  `doFit.py` run (gitignored, harmless).
