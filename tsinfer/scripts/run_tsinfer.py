#!/usr/bin/env python3
"""Run tsinfer + tsdate on a VCF-Zarr.

Uses a CONSTANT recombination rate (no per-site recombination map). The ancestral
state is read from a per-simulation ancestral FASTA (the stdpopsim
`*.mutated.ancestral.fa.gz`) and attached to the zarr with
`tsinfer.add_ancestral_state_array`. Output is a dated .trees.
"""
import argparse

import pyfaidx
import zarr
import tsinfer
import tsdate


def main():
    parser = argparse.ArgumentParser(description="tsinfer + tsdate (constant recomb rate, ancestral FASTA)")
    parser.add_argument("--zarr", required=True, help="input VCF zarr (.vcz)")
    parser.add_argument("--fasta", required=True, help="ancestral reference FASTA (may be bgzipped)")
    parser.add_argument("--chrom", required=True, help="contig name in the FASTA (e.g. chr11)")
    parser.add_argument("--output", required=True, help="output dated .trees path")
    parser.add_argument("--recomb-rate", type=float, default=1.2e-8,
                        help="constant recombination rate (per bp per generation)")
    parser.add_argument("--mut-rate", type=float, required=True,
                        help="mutation rate (per bp per generation), used by tsdate")
    parser.add_argument("--threads", type=int, default=4,
                        help="number of threads for inference")
    args = parser.parse_args()

    # Ancestral state from the FASTA. add_ancestral_state_array() takes the full
    # ancestral sequence and aligns it to the zarr's variant_position internally.
    vcf_zarr = zarr.open(args.zarr, mode="a")
    reader = pyfaidx.Fasta(args.fasta)
    ancestral_str = str(reader[args.chrom])

    # Mask the final variant position with N (avoids an out-of-range edge case at
    # the very last site when building the ancestral array).
    last_pos = int(vcf_zarr["variant_position"][-1])
    ancestral_str = ancestral_str[:last_pos] + "N" + ancestral_str[(last_pos + 1):]

    if "ancestral_state" in vcf_zarr:
        del vcf_zarr["ancestral_state"]
    tsinfer.add_ancestral_state_array(vcf_zarr, ancestral_str)
    print(f"ancestral state from {args.fasta} [{args.chrom}] attached to {args.zarr}")

    vdata = tsinfer.VariantData(args.zarr, ancestral_state="ancestral_state")
    inferred_ts = tsinfer.infer(
        vdata,
        recombination_rate=args.recomb_rate,
        num_threads=args.threads,
    )
    print(f"Inferred a genealogy for {inferred_ts.num_samples} (haploid) genomes")

    # tsdate rejects unary nodes that tsinfer leaves in; simplify them away first.
    dated_ts = tsdate.date(inferred_ts.simplify(), mutation_rate=args.mut_rate)
    dated_ts.dump(args.output)
    print(f"Dated tree sequence written to {args.output}")


if __name__ == "__main__":
    main()
