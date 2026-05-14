#!/bin/bash
#SBATCH --job-name=argweaver
#SBATCH --partition=cpuq
#SBATCH --time=48:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=20G
#SBATCH --output=/grid/siepel/home/khalid/argsortium-inference/argweaver/logs/argweaver_%j.out
#SBATCH --error=/grid/siepel/home/khalid/argsortium-inference/argweaver/logs/argweaver_%j.err

ARGWEAVER_DIR=/grid/siepel/home/khalid/argsortium-inference/argweaver

mkdir -p "${ARGWEAVER_DIR}/logs"

CONFIG_BASE=$(basename "${CONFIG_FILE}" .yaml | sed 's/^config\.//')
WORKDIR="${ARGWEAVER_DIR}/.workdirs/${CONFIG_BASE}"
mkdir -p "${WORKDIR}"

cd "${ARGWEAVER_DIR}"

uv run snakemake \
  --snakefile "${ARGWEAVER_DIR}/Snakefile" \
  --directory "${WORKDIR}" \
  --configfile "${CONFIG_FILE}" \
  --use-apptainer \
  --apptainer-args "--bind /grid" \
  --cores 10
