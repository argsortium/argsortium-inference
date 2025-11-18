#!/usr/bin/env python3
import argparse
import tsinfer
import zarr
import pyfaidx
import tsdate
import os

def main():
    parser = argparse.ArgumentParser(
        description="Run tsinfer + tsdate with ancestral state assignment."
    )
    parser.add_argument(
        "--zarr", required=True,
        help="Path to input VCF zarr file (.vcz)"
    )
    parser.add_argument(
        "--fasta", required=True,
        help="Path to FASTA reference file"
    )
    parser.add_argument(
        "--outdir", default=None,
        help="Output directory (default: same folder as input .vcz)"
    )
    parser.add_argument(
        "--recomb-rate", type=float, default=1.36e-8,
        help="Per-site recombination rate (default: 1.36e-8)"
    )
    parser.add_argument(
        "--mut-rate", type=float, default=2.35e-8,
        help="Per-site mutation rate (default: 2.35e-8)"
    )
    parser.add_argument(
        "--threads", type=int, default=4,
        help="Number of threads to use for inference"
    )
    parser.add_argument(
        "--chrom", default="NC_000014.9",
        help="Chromosome ID in FASTA (default: NC_000014.9)"
    )

    args = parser.parse_args()

    # If outdir not given, default to input folder
    input_dir = os.path.dirname(os.path.abspath(args.zarr))
    outdir = args.outdir if args.outdir else input_dir
    os.makedirs(outdir, exist_ok=True)

    # derive output filename (replace .vcz with .trees)
    base = os.path.splitext(os.path.basename(args.zarr))[0]
    if base.endswith(".vcf"):  # safeguard if file is like sample.vcf.vcz
        base = os.path.splitext(base)[0]
    output_path = os.path.join(outdir, base + ".trees")

    # derive output name by stripping ".vcz" (if present)
    #base = os.path.splitext(args.zarr)[0]
    #output_path = base + ".tskit.trees"

    print(f"Input Zarr: {args.zarr}")
    print(f"Output Trees: {output_path}")

    # Open zarr file
    vcf_zarr = zarr.open(args.zarr)

    # Read reference FASTA for ancestral states
    reader = pyfaidx.Fasta(args.fasta)
    ancestral_str = str(reader[args.chrom])

    # Replace last known position with N
    last_pos = vcf_zarr['variant_position'][-1]
    ancestral_str = ancestral_str[:last_pos] + "N" + ancestral_str[(last_pos + 1):]

    # Remove old dataset if present
    if "ancestral_state" in vcf_zarr:
        del vcf_zarr["ancestral_state"]

    # Add ancestral states
    tsinfer.add_ancestral_state_array(vcf_zarr, ancestral_str)

    # Run tsinfer
    vdata = tsinfer.VariantData(args.zarr, ancestral_state="ancestral_state")
    inferred_ts = tsinfer.infer(
        vdata,
        recombination_rate=args.recomb_rate,
        num_threads=args.threads
    )
    print(f"Inferred a genetic genealogy for {inferred_ts.num_samples} (haploid) genomes")

    # Date with tsdate
    inferred_ts_w_dates = tsdate.date(
        inferred_ts.simplify(),
        mutation_rate=args.mut_rate,
        progress=True
    )
    inferred_ts_w_dates.dump(output_path)
    print(f"Dated tree sequence written to {output_path}")

if __name__ == "__main__":
    main()
