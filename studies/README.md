# studies/

Archived results (fit outputs, aggregated plots, driver logs) of the validation
studies for the dimuon fitting framework. The scripts that produced them live in
`../tools/`. Active production outputs are **not** here — they stay in the repo
root (`signal_fits/`, `signal_fits_uaf/`, `signal_comparison/`, `data_fits/`).

## bias_validation/
The main validation that the fit method is unbiased and produces no spurious
signal. Result directories are named `<study>_N<nsigma>_b<binfrac>`:
- `ensemble_*` — background-only GoF ensembles at various N-sigma / bin-fraction.
- `ensemble_sb_*` — signal+background injection ensembles (incl. `_dcb_` DCB signal).
- `smooth_toy_*` — smooth-template toys (method bias with no parent-MC fluctuations).
- `spurious_test_*`, `pull_test*`, `bkg3pct_*`, `fullmc_spur/` — pull / spurious-signal scans.
- `logs/` — driver stdout for the runs above.
- `bias_summary.png`, `signal_pulls*.png` — summary plots.

**Conclusion (see the memory note):** the low-mass "bias" was a parent-MC
statistical dip inherited by all 3% subsamples, not a method bias; the smooth-toy
test shows the method is unbiased, so no spurious-signal systematic is needed.
Validated config: +/-7 sigma window, 0.5 sigma bins.

## config_scan/
Scan over fit configurations (functional forms / F-test thresholds) + its driver log.

## bkg_mc_validation/
Background-only fits on the full / new background MC (`bkg_mc_fits`,
`bkg_mc_fits_new`) used for the goodness-of-fit checks, + driver log.
