#!/bin/bash
# Run adaptive fits on signal-injected mock datasets, throttled.
set -u
source /cvmfs/cms.cern.ch/cmsset_default.sh
cd /uscms_data/d3/oamram/CMSSW_DiMuonX/src
eval `scramv1 runtime -sh`
cd /uscms_data/d3/oamram/CMSSW_DiMuonX/src/DimuonX/dimuonX_fits
mkdir -p mock_fits

MAX_JOBS=3
MASSES=(15 35 55 75)

running=0
for M in "${MASSES[@]}"; do
  if [ "$running" -ge "$MAX_JOBS" ]; then
    wait -n
    running=$((running - 1))
  fi
  echo "Launching mock M$M"
  python3 run_fit_adaptive.py -M "$M" -c bkg_mc_config.json -s signal_fits/2B \
    -i "mock_data/mock_M$M.h5" \
    --m-data-min 10 --m-data-max 80 \
    -o "mock_fits/M$M" > "mock_fits/adaptive_M$M.log" 2>&1 &
  running=$((running + 1))
done

wait
echo "ALL DONE"
