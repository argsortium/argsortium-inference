# argsortium-inference — Context for Claude

Repo of Snakemake workflows wrapping ARG-inference tools (singer, argweaver, relate, threads, tsinfer, …). Each workflow is its own subdirectory with `Snakefile`, `config.yaml`, optional `scripts/`, and `slurm/` (or `slurm_example/`).

## Always use these

- **Python / Snakemake:** invoke via `uv run --project <repo-root> snakemake …`. `pyproject.toml` pins `snakemake==9.14.1`; the venv is `.venv/` at the repo root. **Do not** `pip install` or run a system snakemake.
- **Containers:** all tool calls run inside `shared/container/arg_inference_tools.sif`. Pass `--use-singularity --singularity-args "--bind <fs>"` to snakemake (the bind path depends on the cluster — see below). Snakefiles already declare `container: config.get("singularity")` on the relevant rules.
- **Params CSV** (one row per inference task): canonical example `*/myconfigs/2026_05_11_runs/master_params.csv`. Columns: `uid, vcf_file, mu, recomb_map, output_dir, inference_start, inference_end`. Header on row 1. Path is set per workflow in `config.yaml` under `params_csv:`. All Snakefiles read this CSV at parse time into a `params[uid]` dict.

## Architecture: SLURM array + per-task Snakemake

Each inference run = one SLURM array task that invokes its own Snakemake instance targeting one output file. Snakemake is **not** used as a cluster orchestrator — SLURM handles all scheduling.

**Why:** avoids a long-running Snakemake controller on the login node, survives disconnects, simpler to reason about and resume.

Per-task contract:

- Each task uses `--directory ${OUTDIR}/${UID_VAL}` so its `.snakemake/` metadata is isolated. **Critical** — without this, concurrent tasks collide on shared files like `.snakemake/iocache/latest.pkl`.
- Intermediate per-task files (temp VCFs, etc.) must live under the UID-specific output dir. E.g. singer's `unzip_vcf` writes `{outdir}/{uid}/{uid}.vcf`, not the shared `{outdir}/{uid}.vcf` (which was reused across UIDs in the same region and caused concurrent-write corruption).
- Per-workflow submit scripts live in `<workflow>/slurm/submit_*.sh` (or `slurm_example/` in `singer/`). They read CSV path + MCMC params from `config.yaml`, build the TARGET from one CSV row, then call `uv run snakemake … <TARGET>`.

## Snakefile conventions (important — wildcard ambiguity bites)

Follow the **singer** pattern (`singer/Snakefile`) — proven not to trip snakemake 9.14.1:

- **Only `{outdir}` and `{uid}` as path wildcards.** Constrain `uid = r"[^/]+"` and `outdir = r".+"`.
- **Bake iteration / step numbers as f-string literals at parse time**, not as `{wildcard}`s:
  ```python
  output: f"{{outdir}}/{{uid}}/{{uid}}.{PENULTIMATE}.smc.gz"
  ```
- **Never put `start`/`end` (or other CSV-derived numeric fields) into output filenames.** Pass them through `params:` via lambdas:
  ```python
  start = lambda wc: params[wc.uid]["inference_start"],
  end   = lambda wc: params[wc.uid]["inference_end"],
  ```
  and reference them in the shell.

**Why:** `uid` values in `master_params.csv` often already contain `__<start>_<end>` (e.g. `…mutated__55039445_60039445`). If `{start}` / `{end}` also appear in an output pattern with `\d+` constraints, snakemake's regex has multiple valid splits → it tries one, fails the `params[wc.uid]` lookup, raises `MissingRuleException` internally to backtrack, and in 9.14.1 a separate formatter bug (`fmt_iofile` calling `.is_storage` on a `str`) crashes the whole process. The argweaver Snakefile previously had this bug; it was rewritten to match singer's convention.

**How to apply:** any new workflow should mirror `singer/Snakefile`'s structure. Don't add wildcards for things you can read from `params[wc.uid]`.

## Other recurring patterns

- **VCF subsetting by inference region** (singer): `bgzip + awk`, no tabix index required (VCFs are not indexed):
  ```bash
  bgzip -d -c {input.vcf} \
    | awk -v s={params.start} -v e={params.end} '/^#/ || ($2 >= s && $2 <= e)' \
    > {output.vcf}
  ```

## Clusters

Two clusters are in active use, with **different filesystems and partitions**. Match `--singularity-args "--bind <fs>"`, the `pyproject.toml` path, and the submit script to the cluster you're on.

### CSHL HPC (current default for argweaver runs)

- Repo path: `/grid/siepel/home/khalid/argsortium-inference/`
- Filesystem: `/grid` → use `--singularity-args "--bind /grid"`
- Partition: `cpuq` (default CPU)
- QoS: `cpu_snice` — 7-day max walltime, 100 concurrent jobs, 1000 submissions per user (also: `cpuq_base` 48h/20 concurrent; `cpu_fill` 1d/400 concurrent for fillers; `slow_nice` 30d/20 concurrent low-priority)
- Per task: 1 CPU, 2GB, walltime in `D-HH:MM:SS` form (e.g. `3-12:00:00`)
- Throttle arrays with `%100`. Specify `--qos=cpu_snice` in the sbatch directives.
- Submit: `bash argweaver/slurm/submit_array.sh [start] [end]` (defaults 1 500).
- Doc: `~/slurm/Requesting Resources … _ CSHL HPC.html`.

### Seawulf3 (singer's submit scripts target this)

- Repo path: `/gpfs/projects/SiepelVeeramahGroup/shkhalid/argsortium-inference/`
- Filesystem: `/gpfs` → use `--singularity-args "--bind /gpfs:/gpfs"`
- Partitions:
  - **Milan:** `extended-96core-shared` (3.5d limit, AMD EPYC, 96c/256GB, shared)
  - **Skylake:** `extended-40core-shared` (3.5d limit, Intel Skylake, 40c/192GB, shared)
- Per task: 1 core, 2GB, up to 36h — use shared queues to avoid wasting full nodes
- **Do not** use `long-96core-shared` (24h max, too short for 36h runs)
- Submit:
  ```bash
  module load slurm/seawulf3/21.08.8
  bash singer/slurm/submit_milan.sh   [start] [end]
  bash singer/slurm/submit_skylake.sh [start] [end]
  ```

## Per-workflow notes

- `singer/` — canonical reference for the CSV + array pattern. `singer/CLAUDE.md` has additional pipeline-fix notes (per-task workdir, VCF region subsetting, bind mount).
- `argweaver/` — mirrors singer's convention. Submit via `argweaver/slurm/submit_array.sh`. Was previously broken by the wildcard ambiguity above; now fixed.
- `relate/`, `threads/`, `tsinfer/` — older code; to be migrated to the CSV + array pattern.
- `preprocessing/` — runs before inference; output paths feed `vcf_file` / `recomb_map` columns in the params CSV.

## Quick reference

```bash
# Dry-run one task locally:
uv run --project <repo-root> snakemake \
    --snakefile <workflow>/Snakefile \
    --configfile <workflow>/config.yaml \
    --directory "$TASK_DIR" --cores 1 --dry-run \
    --use-singularity --singularity-args "--bind <fs>" \
    "$TARGET"

# Submit array:
bash <workflow>/slurm/submit_array.sh [start] [end]

# Monitor:
squeue -u $USER
tail -f <workflow>/slurm/logs/<jobname>_<JOBID>_<TASK>.{out,err}
```

## Pending TODOs

- Migrate `relate/`, `threads/`, `tsinfer/` to the CSV + per-task array pattern (currently use per-run config files).
- CSV splitting on Seawulf: singer submit scripts currently use the full `master_params.csv`; plan is to split into Milan/Skylake halves when running all 1000 jobs concurrently.
