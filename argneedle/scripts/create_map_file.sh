#!/bin/bash

input=$1  # e.g. chr14.27572179.28322179.CHB_301.YRI_0.seed_32.no_multi_allelics.haps
filename=${input%.haps}.map

awk 'BEGIN{OFS="\t"} {print $1, $2, 0, $3}' "$input" > "$filename"

