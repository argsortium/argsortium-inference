# Singer Pipeline — Context for Claude

## Architecture: SLURM array + per-task Snakemake

Each singer run = one SLURM array task that invokes its own Snakemake instance targeting one output file. Snakemake is NOT used as a cluster orchestrator — SLURM handles all scheduling.

**Why:** Avoids a long-running Snakemake controller on the login node, survives disconnects, simpler.

## Key files

- **Snakefile:** `Snakefile`
- **config.yaml:** `config.yaml`
- **params CSV:** `myconfigs/2026_05_11_runs/master_params.csv` (1000 rows, header on row 1)
- **Submit (Milan):** `slurm/submit_milan.sh`
- **Submit (Skylake):** `slurm/submit_skylake.sh`
- **pyproject.toml:** `/gpfs/projects/SiepelVeeramahGroup/shkhalid/argsortium-inference/pyproject.toml` (snakemake 9.14.1, run via `uv run`)

## Params CSV columns

`uid, vcf_file, mu, recomb_map, output_dir, inference_start, inference_end`

## SLURM queues

- **Milan:** `extended-96core-shared` (3.5d limit, AMD EPYC, 96 cores/256GB per node, shared)
- **Skylake:** `extended-40core-shared` (3.5d limit, Intel Skylake, 40 cores/192GB per node, shared)
- Each singer run: 1 core, 2GB, up to 36h — use shared queues to avoid wasting full nodes
- DO NOT use `long-96core-shared` — its max time is 24h, not enough for 36h runs

## Submit script usage

```bash
module load slurm/seawulf3/21.08.8
bash slurm/submit_milan.sh [start] [end]    # default: 1 700
bash slurm/submit_skylake.sh [start] [end]  # default: 1 700
```

CSV path, N, step_size, and FINAL_TSKIT are all read from `config.yaml` at submission time.

## Critical pipeline fixes

### 1. Per-task working directory (fixes .snakemake/ race condition)
Each Snakemake invocation uses `--directory ${OUTDIR}/${UID_VAL}` (the UID-specific output dir). This keeps `.snakemake/` metadata isolated so concurrent tasks don't collide on shared files like `iocache/latest.pkl`.

### 2. VCF written to UID-specific directory (fixes concurrent write conflict)
`unzip_vcf` rule writes the temp VCF to `{outdir}/{uid}/{uid}.vcf` — inside the UID output dir, which is unique per task. The old path `{outdir}/{uid}.vcf` was shared across UIDs in the same region/chr.

### 3. VCF subsetting by inference region
Subsets to `inference_start:inference_end` using bgzip + awk (no tabix index required — VCFs are not indexed):
```bash
bgzip -d -c {input.vcf} | awk -v s={params.start} -v e={params.end} '/^#/ || ($2 >= s && $2 <= e)' > {output.vcf}
```

### 4. Singularity bind mount
All singularity calls require `--bind /gpfs:/gpfs` to access the GPFS filesystem inside the container.


## Pending TODOs

- **Older pipeline version uses per-run config files** instead of the CSV — will be updated later to match this CSV-based approach
- **CSV splitting:** Both submit scripts currently use the full `master_params.csv`. Plan is to split into milan/skylake halves when submitting all 1000 jobs concurrently
