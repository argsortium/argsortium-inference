#!/usr/bin/env python3

import sys
import argparse
from cyvcf2 import VCF


def parse_args():
    parser = argparse.ArgumentParser(
        description="Convert a phased VCF to a .sites file format."
    )
    parser.add_argument(
        "--vcf", required=True, help="Input phased VCF file."
    )
    parser.add_argument(
        "--out", required=True, help="Output sites file name."
    )
    return parser.parse_args()


def main():
    args = parse_args()
    vcf_file = args.vcf
    output_file = args.out

    # Read VCF
    vcf = VCF(vcf_file)

    # Get region info (first contig + min/max POS)
    chrom = None
    positions = []
    for variant in vcf:
        if chrom is None:
            chrom = variant.CHROM
        positions.append(variant.POS)
    region = [chrom, str(min(positions)), str(max(positions))]

    # Get samples
    samples = vcf.samples
    num_samples = len(samples)

    # Hap names: diploid individuals → two haplotypes each
    hap_names = [f"{name}_{i}" for name in samples for i in range(2)]

    # Extract haplotype sequences
    vcf = VCF(vcf_file)  # re-open to iterate again
    data = []
    for variant in vcf:
        pos = variant.POS
        phased = [gt[:2] for gt in variant.genotypes]
        hap1 = [variant.REF, variant.ALT[0]]
        bases = []
        for gt in phased:
            bases.append(hap1[gt[0]])
            bases.append(hap1[gt[1]])
        data.append([pos, "".join(bases)])

    # Write sites file
    with open(output_file, "w") as f:
        f.write("NAMES\t{}\n".format("\t".join(hap_names)))
        f.write("REGION\t{}\n".format("\t".join(region)))
        for row in data:
            f.write("{}\t{}\n".format(row[0], row[1]))

    print(f"Wrote {output_file}")


if __name__ == "__main__":
    main()

