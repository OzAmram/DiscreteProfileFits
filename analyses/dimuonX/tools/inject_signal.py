"""
Create mock datasets = background MC + injected 2B signal events, sized to give
approximately 5 sigma discovery significance.

Significance estimate (counting, peak on locally-flat background):
    Z ~ S_win / sqrt(B_win)   evaluated in a +-2 sigma window around the peak.
To reach a target Z, inject  N_inject = Z*sqrt(B_win) / f_win  total signal
events, where f_win is the fraction of signal MC events falling in +-2 sigma.
"""
import argparse
import json
import os

import h5py
import numpy as np


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("-M", "--mass", type=int, required=True)
    ap.add_argument("--bkg", default="bkg_mc_masses.h5")
    ap.add_argument("--sig-dir", default="signal_data/2B")
    ap.add_argument("--interp-dir", default="signal_fits/2B_loosemass")
    ap.add_argument("--target-z", type=float, default=5.0)
    ap.add_argument("--out-dir", default="mock_data")
    ap.add_argument("--seed", type=int, default=42)
    args = ap.parse_args()

    rng = np.random.default_rng(args.seed)
    M = args.mass

    # signal shape params (peak center + width)
    with open(os.path.join(args.interp_dir, f"case_interpolation_M{float(M)}.json")) as f:
        sp = json.load(f)
    mean, sigma = sp["mean"], sp["sigma"]

    # load background and signal masses
    with h5py.File(args.bkg) as f:
        bkg = f["masses"][:]
    with h5py.File(os.path.join(args.sig_dir, f"signal_mass_{M}.h5")) as f:
        sig = f["masses"][:]

    # +-2 sigma window centered on the fitted peak
    lo, hi = mean - 2 * sigma, mean + 2 * sigma
    B_win = np.count_nonzero((bkg >= lo) & (bkg <= hi))
    f_win = np.count_nonzero((sig >= lo) & (sig <= hi)) / len(sig)

    S_win   = args.target_z * np.sqrt(B_win)        # signal events needed in window
    N_inject = int(round(S_win / f_win))            # total signal events to inject

    replace = N_inject > len(sig)
    idx = rng.choice(len(sig), size=N_inject, replace=replace)
    inj = sig[idx]

    mock = np.concatenate([bkg, inj]).astype(np.float32)
    rng.shuffle(mock)

    os.makedirs(args.out_dir, exist_ok=True)
    out = os.path.join(args.out_dir, f"mock_M{M}.h5")
    with h5py.File(out, "w") as f:
        f.create_dataset("masses", data=mock, compression="gzip", compression_opts=4, chunks=True)

    # save the injected signal masses alone, so the in-window truth yield can be
    # counted later for any fit window (needed for sig_norm and the pull).
    inj_out = os.path.join(args.out_dir, f"inj_M{M}.h5")
    with h5py.File(inj_out, "w") as f:
        f.create_dataset("masses", data=inj.astype(np.float32),
                         compression="gzip", compression_opts=4, chunks=True)

    Z_pred = S_win / np.sqrt(B_win)
    print(f"M={M}: peak={mean:.2f} sigma={sigma:.3f}  window=[{lo:.2f},{hi:.2f}]")
    print(f"  B(+-2s)={B_win}   f_win={f_win:.3f}   sqrt(B)={np.sqrt(B_win):.1f}")
    print(f"  S_in_window={S_win:.0f}   N_inject={N_inject} (replace={replace})")
    print(f"  predicted Z ~ {Z_pred:.2f}")
    print(f"  wrote {out}  ({len(mock):,} events, {N_inject} signal)")


if __name__ == "__main__":
    main()
