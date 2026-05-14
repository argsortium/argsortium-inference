#!/usr/bin/env python3
"""Convert a HapMap recombination map to the 4-column bedGraph format expected
by arg-sample's -R flag: chrom, start, end, rate (per bp per generation).

HapMap rate (cM/Mb) → per bp per generation: multiply by 1e-8.
"""

import argparse


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--hapmap", required=True, help="Input HapMap file")
    parser.add_argument("--out", required=True, help="Output bedGraph file")
    args = parser.parse_args()

    rows = []
    with open(args.hapmap) as f:
        next(f)  # skip header
        for line in f:
            fields = line.strip().split("\t")
            chrom = fields[0]
            pos = int(fields[1])
            rate_per_bp = float(fields[2]) * 1e-8  # cM/Mb → per bp per gen
            rows.append((chrom, pos, rate_per_bp))

    with open(args.out, "w") as f:
        for i, (chrom, start, rate) in enumerate(rows):
            end = rows[i + 1][1] if i + 1 < len(rows) else start + 1
            f.write(f"{chrom}\t{start}\t{end}\t{rate:.6e}\n")

    print(f"Wrote {len(rows)} intervals to {args.out}")


if __name__ == "__main__":
    main()
