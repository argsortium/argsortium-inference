#!/bin/bash
prefix=$1; mu=2.35e-8; Ne=40000

dir=$(dirname "$prefix"); base=$(basename "$prefix" .vcf)

~/relate/relate_v1.2.2_x86_64_static/bin/RelateFileFormats --mode ConvertFromVcf --haps "$dir/$base.haps" --sample "$dir/$base.sample" -i "$dir/$base"
~/relate/relate_v1.2.2_x86_64_static/bin/Relate --mode All -m $mu -N $Ne --haps "$dir/$base.haps" --sample "$dir/$base.sample" --map genetic_map_chr14.27572179.28322179 -o "$base.inferred_arg"
~/relate/relate_lib/bin/Convert --mode ConvertToTreeSequence --anc "$base.inferred_arg.anc" --mut "$base.inferred_arg.mut" -o "$base.inferred_arg.tskit"

echo "Done: ${base}.inferred_arg.tskit"
