#!/bin/bash
# Re-run the fixed smc_to_tskit converter on 6 argweaver runs
# (CEU_0_CHB_20_YRI_0 and ADMIX_20_AFR_0_EUR_0, seeds 1-3; chr11 region; iter 3980).
# Backs up the existing .trees to .trees.bak (once) before overwriting.
set -uo pipefail

SCRIPT=/grid/siepel/home/khalid/argsortium-inference/argweaver/scripts/smc_to_tskit.py
AWBASE=/grid/siepel/home/khalid/argsortium-outputs/2026_05_11_simulation_sets/inferred_args/argweaver
contig=chr11_116827019_120727019
ITER=3980

run_one() {
    local pop=$1 seed=$2
    local cdir="$AWBASE/$pop/$contig"
    local smc sites out
    smc=$(ls "$cdir"/*seed${seed}.mutated/*.${ITER}.smc.gz 2>/dev/null | head -1)
    sites=$(ls "$cdir"/*seed${seed}.mutated.sites 2>/dev/null | head -1)
    if [[ -z "$smc" || -z "$sites" ]]; then
        echo "!! SKIP $pop seed$seed  (smc='$smc' sites='$sites')"; return 1
    fi
    out="${smc%.smc.gz}.trees"
    [[ -f "$out" && ! -f "$out.bak" ]] && cp -p "$out" "$out.bak"
    echo "================================================================"
    echo ">>> $pop seed$seed"
    echo "    smc:   $smc"
    echo "    sites: $sites"
    echo "    out:   $out"
    python3 "$SCRIPT" "$smc" "$out" "$sites"
}

rc=0
for seed in 1 2 3; do run_one CEU_0_CHB_20_YRI_0    "$seed" || rc=1; done
for seed in 1 2 3; do run_one ADMIX_20_AFR_0_EUR_0  "$seed" || rc=1; done
echo "================================================================"
echo "done (rc=$rc)"
exit $rc
