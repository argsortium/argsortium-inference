#!/bin/bash
# Submit a singer chunk to extended-40core dedicated (192GB, 7d limit).
# Limits concurrency to 19 jobs (19 x 10GB = 190GB).
# Usage: bash submit_extended40_dedicated_10g.sh <chunk_csv>

CHUNK_CSV="${1:?Usage: $0 <chunk_csv>}"
CHUNK_CSV=$(realpath "$CHUNK_CSV")
CHUNK_NAME=$(basename "$CHUNK_CSV" .csv)

WORKFLOW_DIR="/gpfs/projects/SiepelVeeramahGroup/shkhalid/argsortium-inference/singer"
PROJECT_DIR="/gpfs/projects/SiepelVeeramahGroup/shkhalid/argsortium-inference"
CONFIG="${WORKFLOW_DIR}/config.yaml"
TASK_DIR="${WORKFLOW_DIR}/slurm/chunks/${CHUNK_NAME}"

PARTITION="extended-40core"
MAX_TIME="168:00:00"
MEM="185G"
CPUS=40
CORES=18  # 185GB / 10GB per job

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
    --cores ${CORES} \\
    --use-singularity \\
    --singularity-args "--bind /gpfs:/gpfs"
EOF
