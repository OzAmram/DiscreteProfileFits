# Fitting framework for F-test + Discrete profile fits

The goal of this repo is to provide a universal, easy-to-use framework for performing parametric fits with multiple families of functional forms.
Signal events are fit to a parametric shape (double crystal-ball) and saved.
An input dataset along with a mass range to be fit is provided.
The data is then fit with the various choices of functional forms describing the background shape.
For each functional form an F-test is performed to select the optimal number of parameters.
Then, a signal plus background fit is performed using Combine. Discrete profiling is used to include all functional forms, each with their determined optimal number of parameters.

## Installation

This framework relies on Combine. Follow the latest recommendations at https://cms-analysis.github.io/HiggsAnalysis-CombinedLimit/latest/#within-cmssw-recommended-for-cms-users

The instructions for the current recommended version (v10.4.2) are:

```bash
cmsrel CMSSW_14_1_0_pre4
cd CMSSW_14_1_0_pre4/src
cmsenv
git -c advice.detachedHead=false clone --depth 1 --branch v10.4.2 \
    https://github.com/cms-analysis/HiggsAnalysis-CombinedLimit.git HiggsAnalysis/CombinedLimit
scramv1 b clean; scramv1 b -j$(nproc --ignore=2)
```

Clone this repo inside the `src/` directory:

```bash
git clone git@github.com:OzAmram/DiscreteProfileFits.git
```

## Inputs

Input files containing masses to be fit must be HDF5 files with a single dataset called `masses`.

Signal shapes are stored in `.json` files produced by `fit_signalshapes.py`.

## Fit configuration file

Both `doFit.py` and `fit_signalshapes.py` accept a JSON config file (`-c / --config`) that sets
all fit parameters. Command-line arguments override config values when both are given.

Example config (`dimuonX_config.json`):

```json
{
    "inputFile":    "bkg_mc_masses.h5",
    "m_min":        30.0,
    "m_max":        36.1,
    "bin_size":     0.05,
    "dcb_model":    true,
    "sqrts":        "13.6",
    "lumi":         "",
    "scale_j_unc":  0.01,
    "res_j_unc":    0.035,
    "sig_norm":     100.0,
    "ftest_thresh": 0.05,
    "err_thresh":   0.5,
    "func_forms": {
        "bern":    [2, 3, 4],
        "polyExp": [1, 2, 3],
        "exp":     [1, 2, 3],
        "expPoly": [1, 2, 3, 4]
    }
}
```

`func_forms` lists the polynomial **order** to try for each family, not the number
of free parameters — the two differ per family (`Utils.get_nPars`): `bern` → order,
`polyExp` → order+1, `expPoly` → order, `exp` → 2*order-1 (`order` slopes plus
order-1 coefficients, the last fraction being fixed by normalization). The orders
above cap every family at 4 free parameters, except `exp` at 5. Allowing too much
flexibility lets the background absorb genuine signal-like features, so prefer
capping these rather than adding orders.

Key fields:

| Field | Description |
|-------|-------------|
| `inputFile` | HDF5 file with a `masses` dataset |
| `m_min`, `m_max` | Fit window boundaries (GeV) |
| `bin_size` | Histogram bin width (GeV). A good rule of thumb is `0.25 × σ_signal` |
| `dcb_model` | Use double crystal-ball signal model (recommended) |
| `sig_norm` | Signal yield corresponding to signal strength `r = 1` |
| `ftest_thresh` | p-value threshold for the F-test to prefer a higher-order function |
| `err_thresh` | Fractional fit-error threshold; models with large errors are excluded from F-test |
| `func_forms` | Dict mapping functional-form name → list of orders to try. Available forms: `bern` (Bernstein polynomial), `exp` (sum of exponentials), `polyExp` (polynomial × exponential), `expPoly` (exponential of a polynomial) |

## Available background functional forms

| Name | Formula |
|------|---------|
| `bern` | Bernstein polynomial of degree `order−1` |
| `exp` | Sum of `order` exponentials |
| `polyExp` | `(1 + Σ pᵢ xⁱ) × exp(c x)` |
| `expPoly` | `exp(Σ pᵢ xⁱ)` |

## Signal shape fitting

Fit signal MC masses at each mass point:

```bash
python3 fit_signalshapes.py \
    -i signal_data/2G/signal_mass_40.h5 \
    -M 40 \
    --dcb-model \
    --no-error-band \
    --m-data-max 80.0 \
    -o signal_fits/2G/
```

This produces `signal_fits/2G/sig_fit_40.json` with the fitted DCB parameters.

### Signal shape interpolation

After fitting all mass points, interpolate to produce signal templates at arbitrary mass
spacing:

```bash
python3 interpolation.py \
    -s case \
    --json-dir signal_fits/2G/ \
    --masses $(python3 -c "import numpy as np; print(*np.round(np.arange(15,80.01,0.2),1))") \
    -m 15 -M 75 \
    -o signal_fits/2G/
```

This produces one JSON per mass point, e.g. `signal_fits/2G/case_interpolation_M40.0.json`.

## Running a single fit

```bash
python3 doFit.py \
    -c dimuonX_config.json \
    -M 40 \
    --m-min 34.7 --m-max 45.3 \
    --bin-size 0.166 \
    -s signal_fits/2G/case_interpolation_M40.0.json \
    -o fits/M40/
```

A good starting point for the window and bin size at mass `M` with signal width `σ`:
- Window: `[M − 8σ, M + 8σ]` (floor at the data lower edge)
- Bin size: `0.25 × σ`

Results and postfit plots are saved in the output directory. The JSON
`fit_results_{M}.json` contains fit quality metrics including `sbfit_prob` and
`bkgfit_prob` (chi²/ndof p-values for the s+b and background-only postfit plots).

## Adaptive fit wrapper

`run_fit_adaptive.py` automates window-size selection. It starts at ±8σ and shrinks
the window by 0.5σ per side if the s+b postfit chi²  p-value is below threshold,
stopping at ±4σ. The bin size is always set to 0.25σ.

```bash
python3 run_fit_adaptive.py \
    -M 72 \
    -c dimuonX_config.json \
    -s signal_fits/2G \
    -o fits/M72/
```

Options:

| Flag | Default | Description |
|------|---------|-------------|
| `-M` | required | Signal mass hypothesis (GeV) |
| `-c` | `dimuonX_config.json` | Config file |
| `-s` | `signal_fits/2G` | Directory containing `case_interpolation_M{mass}.json` |
| `-o` | `bkg_mc_fits/M{mass}` | Output directory |
| `--pval-thresh` | `0.05` | Minimum acceptable s+b chi² p-value |
| `--n-sigma-start` | `8.0` | Starting half-width in σ units |
| `--n-sigma-min` | `4.0` | Minimum half-width before declaring failure |
| `--n-sigma-step` | `0.5` | Step size for shrinking the window |

Each attempted window size is saved in a subdirectory `window_{N}sig/` for review.
The best passing attempt is copied to `best/`. The script exits with code 0 on success
and 1 on failure.
