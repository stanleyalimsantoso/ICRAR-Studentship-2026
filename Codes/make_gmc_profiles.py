import argparse
import os
import glob
from gmc_profile import build_profile

def make_all_profiles(output_dir, profiles_root, pattern="*_gmc_*.dat"):
    os.makedirs(profiles_root, exist_ok=True)

    files = sorted(glob.glob(os.path.join(output_dir, pattern)))
    print(f"Found {len(files)} GMC .dat files in {output_dir}")

    highest_frac = -1
    highest_file = None

    for file in files:
        print(f"Profiling {os.path.basename(file)}")
        frac = build_profile(file, profiles_root)

        if frac > highest_frac:
            highest_frac = frac
            highest_file = file

    if highest_file is None:
        print("No GMC files processed.")
    else:
        print(f"Largest AGB mass fraction file: {os.path.basename(highest_file)}")
        print(f"AGB mass fraction: {highest_frac:.6e}")


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--output-dir", required=True)
    p.add_argument("--profiles-root", required=True)
    p.add_argument("--pattern", default="*_gmc_*.dat")
    args = p.parse_args()

    make_all_profiles(args.output_dir, args.profiles_root, args.pattern)

if __name__ == "__main__":
    main()
