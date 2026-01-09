# Fitting framework for F-test + Discrete profile fits

The goal of this repo is to provide a ~universal easy to use framework for performing parametric fits with multiple familiies of functional forms.
The different functional forms and their considered range of number of paramers are specified. 
Signal events are fit to a parametetric shape (double crystall-ball) and saved.

An input dataset along with a mass range to be fit is provided.
The data is then fit with the various choices of functional forms describing the background shape.
For each functional form an F-test is performed to select the optimal number of parameters. 
Then, a signal plus background fit is performed using combine. Discrete profiling is used to include all functional forms, each with their determined optimal number of parameters. 

## Installation:

This framework relies on Combine. You should follow the latest recommendations for Combine [link](https://cms-analysis.github.io/HiggsAnalysis-CombinedLimit/latest/#within-cmssw-recommended-for-cms-users)
The instructions for the current recommended version (v10.4.2) are:

```
cmsrel CMSSW_14_1_0_pre4
cd CMSSW_14_1_0_pre4/src
cmsenv
git -c advice.detachedHead=false clone --depth 1 --branch v10.4.2 https://github.com/cms-analysis/HiggsAnalysis-CombinedLimit.git HiggsAnalysis/CombinedLimit
cd HiggsAnalysis/CombinedLimit
scramv1 b clean; scramv1 b -j$(nproc --ignore=2) # always make a clean build, with n - 2 cores on the system
```

You can then clone this repo somewhere inside the `src/` directory of the CMSSW release. 
`git clone git@github.com:OzAmram/DiscreteProfileFits.git` 

## Inputs

The input files containing masses to be fit should be h5 files with a single field called `masses`. 

Signal shapes are stored in `.root` files (for now), obtained with the `fit_signalshapes.py` script. 

## Running the fit:

The signal shape parameters can be obtained by fitting the signal masses with a command like:

`python3 fit_signalshapes.py -i test_signal_masses.h5 --dcb-model -M 15 -o sig_test/`

The `-M` option gives the resonance mass of the signal being fit. Output is stored in the directory specified by `-o`.
This creates a file containing the signal shape parameters in `sig_test/sig_fit_15.root` (TODO: change this to a json for ease of use). 

(TODO: signal shape interpolation between mass points) 

This signal shape can then be used to perform the signal + background fit as:

`python3 doFit.py -i test_masses.h5 -M 15 -s sig_test/sig_fit_15.root --dcb-model --m-min 11 --m-max 19 --bin-size 0.4 -o fit_test/`

The fit results and plots get saved in `fit_test`. 

TODO:
- Signal + background fit plots
- Rescale x-axis on plots to reflect true mass values not 0 to 1
- More functional forms
