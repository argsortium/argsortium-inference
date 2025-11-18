#!/bin/bash
prefix=$1; iterstart=$2; iterend=$3; mu=2.35e-8; Ne=40000;

dir=$(dirname "$prefix"); base=$(basename "$prefix" .vcf)

~/relate/relate_v1.2.2_x86_64_static/bin/RelateFileFormats --mode ConvertFromVcf --haps "$dir/$base.haps" --sample "$dir/$base.sample" -i "$dir/$base"

~/relate/relate_v1.2.2_x86_64_static/bin/Relate --mode All -m $mu -N $Ne --haps "$dir/$base.haps" --sample "$dir/$base.sample" \
	--map ~/arg_utils/genetic_map_chr14.27572179.28322179 -o "$base.inferred_args.sample"$iterstart

# loop from iterstart to iterend
for i in $(seq $iterstart $iterend); do
	echo "iteration number $i"
	~/relate/relate_lib/bin/Convert --mode ConvertToTreeSequence --anc ${base}.inferred_args.sample${i}.anc -o ${base}.inferred_args.sample${i}.tskit \
		--mut ${base}.inferred_args.sample${i}.mut

	~/relate/relate_v1.2.2_x86_64_static/bin/RelateCoalescentRate --mode EstimatePopulationSize -i ${base}.inferred_args.sample${i} -o ${base}.inferred_args.sample${i}

	~/relate/relate_v1.2.2_x86_64_static/bin/RelateMutationRate --mode Avg -i ${base}.inferred_args.sample${i} -o ${base}.inferred_args.sample${i}

	~/relate/relate_v1.2.2_x86_64_static/bin/RelateCoalescentRate --mode ReEstimateBranchLengths \
		-i ${base}.inferred_args.sample${i} --mrate ${base}.inferred_args.sample${i}_avg.rate \
		-m $mu --coal ${base}.inferred_args.sample${i}.coal -o ${base}.inferred_args.sample$((i+1))

	#Thin files if iteration is NOT 0 or divisible by 20
	if [[ $i -ne 0 && $((i % 20)) -ne 0 ]]; then
		rm -f ${base}.inferred_args.sample${i}.tskit.trees \
			${base}.inferred_args.sample${i}.anc
			${base}.inferred_args.sample${i}.mut
	fi
	#remove other intermediate files otherwise too many files are created
	rm -f ${base}.inferred_args.sample${i}_avg.rate \
		${base}.inferred_args.sample${i}.coal \
		${base}.inferred_args.sample${i}.popsize \
		${base}.inferred_args.sample${i}.mutrate \
		${base}.inferred_args.sample${i}.bin 

done
