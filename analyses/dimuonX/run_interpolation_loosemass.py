"""
Interpolate DCB signal-shape parameters for every loosemass signal group.

For each group under signal_fits_loosemass/<group>/ (containing sig_fit_<mass>.json
from run_signal_fits_loosemass.py), call interpolation.py to build the DCB
parameter templates on a fine mass grid. Graphs are anchored on all available fit
points (10-100 GeV); the target grid spans the search range on a 0.1 GeV step, so
the low/high edges are interpolated (not extrapolated) thanks to the loose-mass
training. Outputs land next to the input JSONs:
  <group>/case_interpolation.json          (parameterization vs mass)
  <group>/case_interpolation_M<mass>.json  (per-target-mass DCB params)
  <group>/case_interpolation_<param>.png   (fit/spline diagnostic plots)

Runs groups sequentially (ROOT is not thread-safe for this). Set GRID_MIN/GRID_MAX/
GRID_STEP env vars to change the target grid; FIT_MIN/FIT_MAX for the pol-fit range.
"""
import glob
import json
import os
import subprocess
import sys

import numpy as np

# interpolation.py lives in the general framework, one level up from this
# analysis directory.
FRAMEWORK_DIR = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
OUTROOT = "signal_fits_loosemass"
GRID_MIN = float(os.environ.get("GRID_MIN", "10.0"))
GRID_MAX = float(os.environ.get("GRID_MAX", "80.0"))
GRID_STEP = float(os.environ.get("GRID_STEP", "0.1"))
FIT_MIN = os.environ.get("FIT_MIN", "10")   # pol-fit x-range (all fit points)
FIT_MAX = os.environ.get("FIT_MAX", "100")

masses = np.round(np.arange(GRID_MIN, GRID_MAX + 1e-6, GRID_STEP), 1)
mass_args = [("%g" % m) for m in masses]

groups = sorted(d for d in glob.glob(os.path.join(OUTROOT, "*"))
                if os.path.isdir(d) and glob.glob(os.path.join(d, "sig_fit_*.json")))
print("Interpolating %d groups over %d target masses (%.1f-%.1f step %.1f)"
      % (len(groups), len(masses), GRID_MIN, GRID_MAX, GRID_STEP), flush=True)

results = []
for gdir in groups:
    npts = len(glob.glob(os.path.join(gdir, "sig_fit_*.json")))
    cmd = [sys.executable, os.path.join(FRAMEWORK_DIR, "interpolation.py"), "-s", "case",
           "--json-dir", gdir, "-o", gdir,
           "-m", FIT_MIN, "-M", FIT_MAX, "--masses"] + mass_args
    log = os.path.join(gdir, "interp_log.txt")
    with open(log, "w") as f:
        rc = subprocess.run(cmd, stdout=f, stderr=f).returncode
    ok = os.path.exists(os.path.join(gdir, "case_interpolation.json"))
    print("  %-34s rc=%d  %d fit pts  interp=%s"
          % (os.path.basename(gdir), rc, npts, "ok" if ok else "MISSING"), flush=True)
    results.append((gdir, rc, ok))

nbad = sum(1 for _, rc, ok in results if rc != 0 or not ok)
print("\nDone: %d/%d groups interpolated%s"
      % (len(results) - nbad, len(results),
         "" if nbad == 0 else "  (%d FAILED)" % nbad))
