#!/bin/bash
#SBATCH --job-name=argweaver_CEU0_CHB0_YRI20_chr14_100630758
#SBATCH --partition=cpuq
#SBATCH --time=48:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=20G
#SBATCH --output=/grid/siepel/home/khalid/argsortium-inference/argweaver/logs/argweaver_%j.out
#SBATCH --error=/grid/siepel/home/khalid/argsortium-inference/argweaver/logs/argweaver_%j.err

mkdir -p /grid/siepel/home/khalid/argsortium-inference/argweaver/logs

cd /grid/siepel/home/khalid/argsortium-inference/argweaver

uv run snakemake \
  --configfile myconfigs/2026_05_11_runs/config.CEU_0_CHB_0_YRI_20__chr14_100630758_5000000.yaml \
  --use-apptainer \
  --apptainer-args "--bind /grid" \
  --cores 10
