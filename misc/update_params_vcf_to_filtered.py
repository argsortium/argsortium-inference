#!/usr/bin/env python3
"""Update vcf_file in *.params.csv files to point to preprocessed filtered VCFs.

Replaces the .vcf.gz suffix with .no_multiallelics.filtered.vcf.gz in the
vcf_file column of every matching params CSV. Modifies files in-place.

Usage:
    python update_params_vcf_to_filtered.py <glob_pattern> [<glob_pattern> ...]

Example:
    python update_params_vcf_to_filtered.py \
        "/grid/.../CEU_0_CHB_20_YRI_0/*.params.csv" \
        "/grid/.../CEU_0_CHB_125_YRI_125/*.params.csv"
"""

import csv
import glob
import io
import os
import re
import sys


RAW_SUFFIX = ".vcf.gz"
FILTERED_SUFFIX = ".no_multiallelics.filtered.vcf.gz"


def clean_vcf_path(raw: str) -> str:
    """Extract the VCF path up to .vcf.gz, stripping any trailing garbage."""
    m = re.match(r"(.+?\.vcf\.gz)", raw)
    if not m:
        raise ValueError(f"Cannot parse vcf path from: {raw!r}")
    return m.group(1)


def process_file(csv_path: str) -> None:
    with open(csv_path, newline="") as f:
        reader = csv.DictReader(f)
        fieldnames = reader.fieldnames
        rows = list(reader)

    if "vcf_file" not in fieldnames:
        print(f"  SKIP (no vcf_file column): {csv_path}")
        return

    changed = False
    for row in rows:
        vcf = clean_vcf_path(row["vcf_file"])
        filtered = vcf.replace(RAW_SUFFIX, FILTERED_SUFFIX)
        if not os.path.exists(filtered):
            print(f"  WARNING filtered VCF not found: {filtered}")
        if row["vcf_file"] != filtered:
            row["vcf_file"] = filtered
            changed = True

    if not changed:
        print(f"  already up to date: {os.path.basename(csv_path)}")
        return

    with open(csv_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    print(f"  updated: {os.path.basename(csv_path)}")


def main():
    if len(sys.argv) < 2:
        print(__doc__)
        sys.exit(1)

    paths = []
    for pattern in sys.argv[1:]:
        paths.extend(sorted(glob.glob(pattern)))

    if not paths:
        print("No files matched.")
        sys.exit(1)

    print(f"Processing {len(paths)} params files...")
    for path in paths:
        process_file(path)
    print("Done.")


if __name__ == "__main__":
    main()
