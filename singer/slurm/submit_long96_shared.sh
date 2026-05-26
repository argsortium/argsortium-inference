#!/bin/bash
# Submit singer runs as a SLURM array on long-96core-shared (Milan, 24h limit).
# Each array task runs snakemake for one row of the params CSV.
# Usage: bash submit_long96.sh [start] [end]
#   defaults: 501 1000

WORKFLOW_DIR="/gpfs/projects/SiepelVeeramahGroup/shkhalid/argsortium-inference/singer"
PROJECT_DIR="/gpfs/projects/SiepelVeeramahGroup/shkhalid/argsortium-inference"
CONFIG="${WORKFLOW_DIR}/config.yaml"

# SLURM resource defaults
PARTITION="long-96core-shared"
MAX_TIME="24:00:00"
MEM="2G"
CPUS=1

START="${1:-501}"
END="${2:-1000}"

# Pull CSV path and MCMC params from config.yaml
CSV=$(awk -F': *' '/^params_csv:/ {gsub(/"/, "", $2); print $2}' "$CONFIG")
N=$(awk -F': *' '/^mcmc_samples:/ {print $2}' "$CONFIG")
STEP=$(awk -F': *' '/^step_size:/ {print $2}' "$CONFIG")
FINAL_TSKIT=$(python3 -c "print(list(range(0, $N, $STEP))[-1])")

mkdir -p "${WORKFLOW_DIR}/slurm/logs"

sbatch <<EOF
#!/bin/bash
#SBATCH --job-name=singer_long96
#SBATCH --partition=${PARTITION}
#SBATCH --array=${START}-${END}
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=${CPUS}
#SBATCH --mem=${MEM}
#SBATCH --time=${MAX_TIME}
#SBATCH --output=${WORKFLOW_DIR}/slurm/logs/singer_%A_%a.out
#SBATCH --error=${WORKFLOW_DIR}/slurm/logs/singer_%A_%a.err

set -euo pipefail

ROW=\$(sed -n "\$((SLURM_ARRAY_TASK_ID + 1))p" "${CSV}")
UID_VAL=\$(echo "\$ROW" | cut -d',' -f1)
OUTDIR=\$(echo "\$ROW" | cut -d',' -f5)
TARGET="\${OUTDIR}/\${UID_VAL}/\${UID_VAL}.singer.tskit_${FINAL_TSKIT}.trees"

echo "Task \$SLURM_ARRAY_TASK_ID: \$UID_VAL"
echo "Target: \$TARGET"

TASK_DIR="\${OUTDIR}/\${UID_VAL}"
mkdir -p "\$TASK_DIR"

uv run --project "${PROJECT_DIR}" snakemake \\
    --snakefile "${WORKFLOW_DIR}/Snakefile" \\
    --configfile "${CONFIG}" \\
    --directory "\$TASK_DIR" \\
    --cores ${CPUS} \\
    --use-singularity \\
    --singularity-args "--bind /gpfs:/gpfs" \\
    "\$TARGET"
EOF
