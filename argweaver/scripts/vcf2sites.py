#!/usr/bin/env python3

import sys
import gzip
import argparse


def parse_args():
    parser = argparse.ArgumentParser(
        description="Convert a phased VCF to a .sites file format."
    )
    parser.add_argument("--vcf", required=True, help="Input phased VCF file.")
    parser.add_argument("--out", required=True, help="Output sites file name.")
    return parser.parse_args()


def open_vcf(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def main():
    args = parse_args()

    samples = []
    chrom = None
    positions = []
    data = []

    with open_vcf(args.vcf) as f:
        for line in f:
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                samples = line.strip().split("\t")[9:]
                continue
            fields = line.strip().split("\t")
            chrom = fields[0]
            pos = int(fields[1])
            ref = fields[3]
            alt = fields[4].split(",")[0]
            alleles = [ref, alt]

            fmt = fields[8].split(":")
            gt_idx = fmt.index("GT")

            bases = []
            for sample_field in fields[9:]:
                gt = sample_field.split(":")[gt_idx]
                a1, a2 = gt.split("|")
                if a1 not in ("0", "1") or a2 not in ("0", "1"):
                    raise ValueError(
                        f"Multi-allelic genotype at {chrom}:{pos} — "
                        f"GT={gt}. Run bcftools norm -m - first."
                    )
                bases.append(alleles[int(a1)])
                bases.append(alleles[int(a2)])

            positions.append(pos)
            data.append((pos, "".join(bases)))

    hap_names = [f"{s}_{i}" for s in samples for i in range(2)]
    region = [chrom, str(min(positions)), str(max(positions))]

    with open(args.out, "w") as f:
        f.write("NAMES\t{}\n".format("\t".join(hap_names)))
        f.write("REGION\t{}\n".format("\t".join(region)))
        for pos, bases in data:
            f.write(f"{pos}\t{bases}\n")

    print(f"Wrote {args.out}")


if __name__ == "__main__":
    main()
