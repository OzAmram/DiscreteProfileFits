"""Make a random 3% subsample of the background MC, emulating the anomaly-score
cut (eff ~3%, ~uncorrelated with mass) expected in the real search."""
import h5py, numpy as np
rng = np.random.default_rng(42)
FRAC = 0.03
with h5py.File("bkg_mc_masses.h5") as f:
    m = f["masses"][:]
n = len(m)
keep = rng.random(n) < FRAC
sub = m[keep].astype(np.float32)
with h5py.File("bkg_mc_3pct.h5", "w") as f:
    f.create_dataset("masses", data=sub, compression="gzip", compression_opts=4, chunks=True)
print(f"full={n:,}  kept={len(sub):,}  ({100*len(sub)/n:.3f}%)")
print(f"range [{sub.min():.2f},{sub.max():.2f}]  mean={sub.mean():.2f}")
