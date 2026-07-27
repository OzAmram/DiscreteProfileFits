#!/bin/bash
# Fit signal-injected mocks at a FIXED +-8sigma window (no adaptive shrinking),
# reporting discovery significance. Throttled to MAX_JOBS.
set -u
source /cvmfs/cms.cern.ch/cmsset_default.sh
cd /uscms_data/d3/oamram/CMSSW_DiMuonX/src
eval `scramv1 runtime -sh`
cd /uscms_data/d3/oamram/CMSSW_DiMuonX/src/DimuonX/dimuonX_fits/analyses/dimuonX
mkdir -p mock_fits_fixed

MAX_JOBS=3
# mass m_min m_max bin_size  (+-8 sigma, capped at [10,80], bin=0.25 sigma)
ROWS=(
  "15 13.08 16.92 0.060"
  "35 30.39 39.61 0.144"
  "55 47.69 62.31 0.228"
  "75 65.23 80.00 0.305"
)

running=0
for row in "${ROWS[@]}"; do
  read M MMIN MMAX BS <<< "$row"
  if [ "$running" -ge "$MAX_JOBS" ]; then wait -n; running=$((running-1)); fi
  echo "Launching mock M$M  [$MMIN,$MMAX] bin=$BS"
  out="mock_fits_fixed/M$M"
  mkdir -p "$out"
  python3 ../../doFit.py -c dimuonX_config.json -M "$M" \
    --m-min "$MMIN" --m-max "$MMAX" --bin-size "$BS" \
    -i "mock_data/mock_M$M.h5" \
    -s "signal_fits/2B_loosemass/case_interpolation_M$M.0.json" \
    -o "$out" > "$out/fit_log.txt" 2>&1 &
  running=$((running+1))
done
wait
echo "ALL DONE"
