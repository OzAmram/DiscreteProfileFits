"""
Preprocess the CATHODE_training_loosemass_v1 signal models for signal-shape fitting.

Same as prep_signal_uaf.py but for the loose-mass training (no hard dimuon-mass
cuts at 10/80 GeV), which should give better shape fits at the low- and high-mass
edges. Each postprocessed signal file stores the dimuon mass in
'Background_dimuon_mass'; fit_signalshapes.py wants a plain h5 with key 'masses'.
This extracts the mass column per file into
signal_data_loosemass/<group>/mass_<M>.h5, grouping files by process+heavy-parent
(the resonance mass is the last mass token in the filename). Writes
signal_loosemass_manifest.json: {group: [[mass, h5path], ...]}.
"""
import glob
import json
import os
import re

import h5py
import numpy as np

SRC = "/eos/uscms/store/user/tvami/CATHODE/CATHODE_training_loosemass_v1"
OUT = "signal_data_loosemass"
# resonance-mass token sits at the end: MH3-<n>, MH-<n>, MH2-<n>, or MS<n>
TOKEN = re.compile(r"[-_](MH3-|MH2-|MH-|MS)(\d+)$")

os.makedirs(OUT, exist_ok=True)
manifest = {}
files = sorted(glob.glob(os.path.join(SRC, "signal_*_postprocessed.h5")))
print("found %d signal files" % len(files))
for f in files:
    stem = os.path.basename(f)[len("signal_"):-len("_postprocessed.h5")]
    m = TOKEN.search(stem)
    if not m:
        print("SKIP (no mass token): %s" % stem)
        continue
    mass = int(m.group(2))
    group = stem[:m.start()]
    gdir = os.path.join(OUT, group)
    os.makedirs(gdir, exist_ok=True)
    dst = os.path.join(gdir, "mass_%d.h5" % mass)
    with h5py.File(f, "r") as d:
        arr = np.asarray(d["Background_dimuon_mass"][()]).ravel().astype(np.float32)
    arr = arr[np.isfinite(arr)]
    with h5py.File(dst, "w") as d:
        d.create_dataset("masses", data=arr, compression="gzip", chunks=True)
    manifest.setdefault(group, []).append([mass, dst])

for g in manifest:
    manifest[g].sort()
json.dump(manifest, open("signal_loosemass_manifest.json", "w"), indent=2)

print("\n%-34s %-6s %s" % ("group", "Nfit", "masses"))
print("-" * 70)
tot = 0
for g in sorted(manifest):
    ms = [x[0] for x in manifest[g]]
    tot += len(ms)
    print("%-34s %-6d %s" % (g, len(ms), " ".join(str(x) for x in ms)))
print("-" * 70)
print("total: %d files in %d groups" % (tot, len(manifest)))
print("wrote signal_loosemass_manifest.json")
