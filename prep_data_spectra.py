"""
Copy each window's masked (sideband-included) spectrum into a local h5 with the
dataset key 'masses' that the (validated) fit code expects. Source files use key
'mass'. Writes data_spectra/<window>.h5.
"""
import json
import os

import h5py
import numpy as np

BASE = "/eos/uscms/store/user/tvami/CATHODE/Flow_Learning/Results_vr"
SERIES = "nsplit_20260708_104130"
OUT = "data_spectra"

plan = json.load(open("scan_plan.json"))
os.makedirs(OUT, exist_ok=True)
for tag in plan:
    src = os.path.join(BASE, "%s_%s" % (SERIES, tag),
                       "stage3_classifier_masked", "global_anomalous_masses.h5")
    with h5py.File(src, "r") as f:
        key = "masses" if "masses" in f else "mass"
        arr = np.asarray(f[key][()]).ravel().astype(np.float32)
    dst = os.path.join(OUT, "%s.h5" % tag)
    with h5py.File(dst, "w") as f:
        f.create_dataset("masses", data=arr, compression="gzip", chunks=True)
    print("%s -> %s  (N=%d, %.2f-%.2f GeV)" % (tag, dst, arr.size, arr.min(), arr.max()))
print("done")
