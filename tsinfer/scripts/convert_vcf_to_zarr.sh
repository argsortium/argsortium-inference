#!/bin/bash
# Usage: ./run_vcf2zarr.sh input.vcf.gz /path/to/output_dir


#zarr version should < 3
set -euo pipefail

vcf_file="$1"
out_path="$2"

# ensure output directory exists
mkdir -p "$out_path"

# remove directory and any .vcf* extension
filename=$(basename "$vcf_file")
base=${filename%%.vcf*}

echo "Processing $vcf_file"
echo "Base name: $base"

vcf2zarr explode "$vcf_file" "${out_path}/${base}.icf"
vcf2zarr encode "${out_path}/${base}.icf" "${out_path}/${base}.vcz"

echo "Finished: ${out_path}/${base}.vcz"
