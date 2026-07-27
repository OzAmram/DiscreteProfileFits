#!/bin/bash
# Run the signal-shape mismatch sensitivity test under cmsenv (combine on PATH).
# Usage: ./run_shape_mismatch.sh [--smoke] [extra args to shape_mismatch_test.py]
set -u
source /cvmfs/cms.cern.ch/cmsset_default.sh
cd /uscms_data/d3/oamram/CMSSW_DiMuonX/src
eval `scramv1 runtime -sh`
cd /uscms_data/d3/oamram/CMSSW_DiMuonX/src/DimuonX/dimuonX_fits/analyses/dimuonX
python3 shape_mismatch_test.py "$@"
