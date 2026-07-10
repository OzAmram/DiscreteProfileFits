#!/bin/bash
# Run adaptive fits across mass points with a concurrency cap.
set -u

# --- CMSSW environment (required for ROOT / combine) ---
source /cvmfs/cms.cern.ch/cmsset_default.sh
cd /uscms_data/d3/oamram/CMSSW_DiMuonX/src
eval `scramv1 runtime -sh`
cd /uscms_data/d3/oamram/CMSSW_DiMuonX/src/DimuonX/dimuonX_fits
mkdir -p bkg_mc_fits_new

MAX_JOBS=3

# mass:sig_dir pairs (M13 uses 2B which has coverage down to 10 GeV)
JOBS=(
  "13:signal_fits/2B"
  "17:signal_fits/2G"
  "22:signal_fits/2G"
  "27:signal_fits/2G"
  "32:signal_fits/2G"
  "40:signal_fits/2G"
  "50:signal_fits/2G"
  "57:signal_fits/2G"
  "62:signal_fits/2G"
  "72:signal_fits/2G"
  "77:signal_fits/2G"
)

running=0
for entry in "${JOBS[@]}"; do
  M="${entry%%:*}"
  SIG="${entry##*:}"

  if [ "$running" -ge "$MAX_JOBS" ]; then
    wait -n          # block until any one job finishes
    running=$((running - 1))
  fi

  echo "Launching M$M (sig=$SIG)"
  python3 run_fit_adaptive.py -M "$M" -c bkg_mc_config.json -s "$SIG" \
    --m-data-min 10 --m-data-max 80 \
    -o "bkg_mc_fits_new/M$M" > "bkg_mc_fits_new/adaptive_M$M.log" 2>&1 &
  running=$((running + 1))
done

wait
echo "ALL DONE"
