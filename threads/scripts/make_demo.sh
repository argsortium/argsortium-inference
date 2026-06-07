#!/bin/bash
# Estimate an smc++ demography from a single VCF (on N random samples) and write
# a Threads .demo (time<TAB>2N). Meant to be run INSIDE the tskit-v1 container
# (needs smc++, bcftools, bgzip/tabix, python). Resumable: skips if OUT exists.
#
# Usage: make_demo.sh <vcf.gz> <mu> <out.demo> [n_samples] [n_distinguished] [missing_cutoff]
set -euo pipefail

VCF="$1"; MU="$2"; OUT="$3"
NS="${4:-10}"; ND="${5:-5}"; CUT="${6:-1000000}"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SEED_STR="$(basename "$VCF")"

if [ -s "$OUT" ]; then
    echo "skip (exists): $OUT"
    exit 0
fi

WD="$(mktemp -d /tmp/smcdemo.XXXXXX)"
trap 'rm -rf "$WD"' EXIT
cd "$WD"

# smc++ needs an indexed VCF (random access by contig)
ln -sf "$VCF" in.vcf.gz
tabix -f -p vcf in.vcf.gz

# N random samples, deterministic per VCF -> POP spec
SAMPLES=$(bcftools query -l in.vcf.gz | shuf --random-source=<(yes "$SEED_STR") -n "$NS" | paste -sd, -)
CONTIG=$(tabix -l in.vcf.gz | head -1)
echo "[$SEED_STR] contig=$CONTIG samples=$SAMPLES"

# one vcf2smc per distinguished individual (robust to per-individual degeneracy)
i=0
for d in $(echo "$SAMPLES" | tr ',' ' ' | head -"$ND"); do
    smc++ vcf2smc --missing-cutoff "$CUT" -d "$d" "$d" \
        in.vcf.gz data.$i.smc.gz "$CONTIG" "POP:$SAMPLES"
    i=$((i + 1))
done

smc++ estimate -o "$WD" "$MU" data.*.smc.gz
smc++ plot --csv plot.png model.final.json

mkdir -p "$(dirname "$OUT")"
python "$SCRIPT_DIR/smc_demo_to_threads.py" plot.csv "$OUT"
echo "wrote $OUT"
