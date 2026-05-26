#!/bin/bash
# Submit a singer chunk to long-96core-shared (shared, 24h limit).
# Snakemake runs jobs in parallel internally up to --cores.
# Usage: bash submit_long96_shared_chunk.sh <chunk_csv> [cores]
# Example: bash submit_long96_shared_chunk.sh slurm/chunks/chunk9.csv 77

CHUNK_CSV="${1:?Usage: $0 <chunk_csv> [cores]}"
CHUNK_CSV=$(realpath "$CHUNK_CSV")
CHUNK_NAME=$(basename "$CHUNK_CSV" .csv)
# Number of cores = number of data rows in the chunk (header excluded)
CORES="${2:-$(( $(wc -l < "$CHUNK_CSV") - 1 ))}"

WORKFLOW_DIR="/gpfs/projects/SiepelVeeramahGroup/shkhalid/argsortium-inference/singer"
PROJECT_DIR="/gpfs/projects/SiepelVeeramahGroup/shkhalid/argsortium-inference"
CONFIG="${WORKFLOW_DIR}/config.yaml"
TASK_DIR="${WORKFLOW_DIR}/slurm/chunks/${CHUNK_NAME}"

# SLURM resource defaults
PARTITION="long-96core-shared"
MAX_TIME="24:00:00"
MEM=$(( CORES * 2 ))G
CPUS=${CORES}

mkdir -p "${WORKFLOW_DIR}/slurm/logs" "$TASK_DIR"

sbatch <<EOF
#!/bin/bash
#SBATCH --job-name=singer_${CHUNK_NAME}
#SBATCH --partition=${PARTITION}
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=${CPUS}
#SBATCH --mem=${MEM}
#SBATCH --time=${MAX_TIME}
#SBATCH --output=${WORKFLOW_DIR}/slurm/logs/singer_${CHUNK_NAME}_%j.out
#SBATCH --error=${WORKFLOW_DIR}/slurm/logs/singer_${CHUNK_NAME}_%j.err

set -euo pipefail

echo "Running chunk: ${CHUNK_CSV}"
echo "Working dir: ${TASK_DIR}"
echo "Cores: ${CPUS}"

uv run --project "${PROJECT_DIR}" snakemake \\
    --snakefile "${WORKFLOW_DIR}/Snakefile" \\
    --configfile "${CONFIG}" \\
    --config params_csv="${CHUNK_CSV}" \\
    --directory "${TASK_DIR}" \\
    --cores ${CPUS} \\
    --use-singularity \\
    --singularity-args "--bind /gpfs:/gpfs"
EOF
