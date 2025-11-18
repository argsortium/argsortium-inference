input_hap=$1
demo_file=$2
decoding_quant=$3
samples=$4

chrom=1
mapfile=${input_hap%.haps}.map
outfile=$(basename "${input_hap%.haps}")
echo $outfile
echo $mapfile

PATH_2_CONVERTER="/grid/siepel/home/khalid/arg_metrics/argneedle/"

arg_needle \
  --mode sequence \
  --hap_gz $input_hap \
  --normalize 1 \
  --normalize_demography $demo_file \
  --asmc_decoding_file $decoding_quant \
  --chromosome $chrom \
  --map  $mapfile \
  --out  $outfile \
  --num_sequence_samples $samples \
  --verbose 1

${PATH_2_CONVERTER}argn_to_tskit.py ${outfile}.argn
