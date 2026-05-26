#!/bin/bash
# Submit a singer chunk to extended-96core (dedicated, non-shared, 7d limit).
# Snakemake runs up to 96 jobs in parallel internally.
# Usage: bash submit_extended96_dedicated.sh <chunk_csv>
# Example: bash submit_extended96_dedicated.sh slurm/chunks/chunk1.csv

CHUNK_CSV="${1:?Usage: $0 <chunk_csv>}"
CHUNK_CSV=$(realpath "$CHUNK_CSV")
CHUNK_NAME=$(basename "$CHUNK_CSV" .csv)

WORKFLOW_DIR="/gpfs/projects/SiepelVeeramahGroup/shkhalid/argsortium-inference/singer"
PROJECT_DIR="/gpfs/projects/SiepelVeeramahGroup/shkhalid/argsortium-inference"
CONFIG="${WORKFLOW_DIR}/config.yaml"
TASK_DIR="${WORKFLOW_DIR}/slurm/chunks/${CHUNK_NAME}"

# SLURM resource defaults
PARTITION="extended-96core"
MAX_TIME="168:00:00"  # 7 days
MEM="192G"            # 96 tasks x 2GB
CPUS=96

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

uv run --project "${PROJECT_DIR}" snakemake \\
    --snakefile "${WORKFLOW_DIR}/Snakefile" \\
    --configfile "${CONFIG}" \\
    --config params_csv="${CHUNK_CSV}" \\
    --directory "${TASK_DIR}" \\
    --cores ${CPUS} \\
    --use-singularity \\
    --singularity-args "--bind /gpfs:/gpfs"
EOF
