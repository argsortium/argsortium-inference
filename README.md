# ARGsortium Inference Workflows

Snakemake workflows for ARG inference and preprocessing. Each inference pipeline consumes a params CSV (one row per contig/region) so runs can be parallelized cleanly. Current Snakemake implementations cover preprocessing, Singer, and tsinfer+tsdate; helper scripts for ArgWeaver, Relate, and ArgNeedle are included for future wiring.

## Repository layout
- `preprocessing/` – VCF cleanup (normalize, frequency calc, variant picking) before inference.
- `tsinfer/` – Snakemake pipeline for tsinfer + tsdate (`scripts/run_tsinfer.py`).
- `singer/` – Snakemake pipeline wrapping Singer and conversion to tskit.
- `argweaver/`, `relate/`, `argneedle/` – helper scripts (conversion, format prep); no Snakefiles yet.
- `shared/container/arg_inference_tools.def|.sif` – Singularity recipe/image with toolchain (htslib/bgzip, bcftools, plink2, argweaver, Singer, Relate, Python deps).
- `pyproject.toml` – pins Snakemake (`>=3.12` Python).

## Prerequisites
- Singularity/Apptainer is required; build the bundled image:
  ```bash
  singularity build shared/container/arg_inference_tools.sif shared/container/arg_inference_tools.def
  ```
- Python 3.12+; Snakemake version is pinned via `uv` in this repo (`pyproject.toml`).
- Input data: gzipped, phased VCFs for inference workflows; reference FASTA, recombination maps, and params CSV as described below.

## Params CSV schema (used by Singer + tsinfer workflows)
Provide one row per contig/region you want to infer. Columns used today:

| column          | description                                   | used by      |
|-----------------|-----------------------------------------------|--------------|
| `uid`           | unique ID used in output filenames            | all          |
| `vcf_file`      | path to the (bgzipped) VCF for this region    | all          |
| `mu`            | per-site mutation rate                        | all          |
| `recomb_map`    | path to recombination map (HapMap-style)      | all          |
| `contig`        | contig/chr name (for FASTA lookup)            | tsinfer      |
| `ancestral_fasta` | FASTA with ancestral states                  | tsinfer      |
| `seed`          | random seed (not yet consumed in code)        | tsinfer      |
| `start`, `end`  | region bounds (passed to Singer; not used elsewhere) | Singer |

Add any extra columns you need; Snakemake will ignore unused fields.

## Preprocessing workflow
Normalizes VCFs, computes frequencies, selects one biallelic variant per position (with optional exclusions), and outputs filtered VCFs ready for inference.

Key config (`preprocessing/config.yaml`):
- `vcf_pattern`: glob for input VCFs.
- `output_dir`: where cleaned files are written.
- `add_chr`: whether to prefix `chr` in variant IDs when writing inclusion lists.
- `exclusion_list`: optional file with variants (`chr:pos:ref:alt`, one per line) to filter out for QC/MAF concerns.
- `singularity`: path to the container image.

Run (use `uv` to respect the pinned Snakemake version):
```bash
uv run snakemake -s preprocessing/Snakefile --configfile preprocessing/config.yaml --use-singularity --singularity-args '--bind /path:/path' --cores 4
```
Outputs: `{output_dir}/{sample}.no_multiallelics.filtered.vcf.gz` plus frequency and inclusion list files.

## tsinfer + tsdate workflow
Two rules: convert each VCF to Zarr with `vcf2zarr`, then run `scripts/run_tsinfer.py` to infer and date the tree sequence.

Key config (`tsinfer/config.yaml`):
- `params_file`: path to the params CSV (see schema).
- `output_dir`: destination for Zarrs and `.tsinfer.trees`.
- `threads`: threads for tsinfer (default 4).
- `singularity`: container path.

Run (via `uv`):
```bash
uv run snakemake -s tsinfer/Snakefile --configfile tsinfer/config.yaml --use-singularity --singularity-args '--bind /path:/path' --cores 8
```
Outputs: `{output_dir}/{uid}.vcz` and `{output_dir}/{uid}.tsinfer.trees`.

## Singer workflow
Unzips VCFs, runs `singer_master` for MCMC samples, then converts to tskit trees with `convert_to_tskit`.

Key config (`singer/config.yaml`):
- `params_file`, `output_dir`, `singularity` as above.
- `mcmc_samples`: number of Singer iterations.
- `step_size`: stride used when choosing which iteration to convert to tskit.
- `resume`: pass `-resume` to Singer if continuing a run.
- `Ne`, `polar`: Singer population size and polarization parameters.

Run (via `uv`):
```bash
uv run snakemake -s singer/Snakefile --configfile singer/config.yaml --use-singularity --singularity-args '--bind /path:/path' --cores 8
```
Outputs: Singer text outputs and `{output_dir}/{uid}.singer.tskit_<iter>.trees`.

## Other tools
- ArgWeaver helpers: `argweaver/scripts/vcf2sites.py` (VCF→.sites) and `argweaver/scripts/argweaver_to_tskit.py` (conversion; marked experimental).
- Relate helpers: `relate/scripts/runRelateVcf2tskit.sh` and `runRelateWBranchResampling.sh` for end-to-end runs and branch resampling.
- ArgNeedle helpers: `argneedle/scripts/create_map_file.sh`, `haps2tskit.sh`, `argn_to_tskit.py`.

These scripts are provided as references and are not yet wired into Snakemake.

## Tips
- All workflows expect bgzipped VCFs; preprocessing will normalize IDs and ensure biallelic variants.
- Set `output_dir` in each config to keep outputs separate per tool/run.
- Use `--cores` to parallelize across `uid` values from the params CSV.
