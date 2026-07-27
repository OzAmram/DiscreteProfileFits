import argparse
import glob
import os
import re

import h5py
import numpy as np
import uproot

FILENAME_RE = re.compile(r"To2Mu(2G|Inv|2B).*-MS(\d+)_")


def parse_filename(fname):
    m = FILENAME_RE.search(os.path.basename(fname))
    if m:
        return m.group(1), int(m.group(2))
    return None, None


def extract_masses(root_file):
    with uproot.open(root_file) as f:
        arr = f["Events"]["dimu_mass"].array(library="np")
    return arr.astype(np.float32)


def main():
    parser = argparse.ArgumentParser(
        description="Extract dimu_mass from signal ROOT files and save to H5"
    )
    parser.add_argument("--input-dir", required=True,
                        help="Directory containing signal ROOT files")
    parser.add_argument("--output-dir", default="signal_data",
                        help="Root output directory (creates {output-dir}/{variant}/signal_mass_{mass}.h5)")
    parser.add_argument("--variants", nargs="+", default=["2G", "Inv", "2B"],
                        help="Which signal variants to process (default: 2G Inv 2B)")
    args = parser.parse_args()

    root_files = sorted(glob.glob(os.path.join(args.input_dir, "*.root")))
    print(f"Found {len(root_files)} ROOT files in {args.input_dir}")

    processed = 0
    for root_file in root_files:
        variant, mass = parse_filename(root_file)
        if variant is None:
            print(f"  Skipping (unrecognised filename): {os.path.basename(root_file)}")
            continue
        if variant not in args.variants:
            continue

        out_dir = os.path.join(args.output_dir, variant)
        os.makedirs(out_dir, exist_ok=True)
        out_h5 = os.path.join(out_dir, f"signal_mass_{mass}.h5")

        print(f"  {variant} MS{mass:3d}  {os.path.basename(root_file)}")
        masses = extract_masses(root_file)
        print(f"    {len(masses)} events  mass in [{masses.min():.2f}, {masses.max():.2f}] GeV")

        with h5py.File(out_h5, "w") as hf:
            hf.create_dataset("masses", data=masses)

        processed += 1

    print(f"\nDone — wrote {processed} H5 files to {args.output_dir}/")


if __name__ == "__main__":
    main()
