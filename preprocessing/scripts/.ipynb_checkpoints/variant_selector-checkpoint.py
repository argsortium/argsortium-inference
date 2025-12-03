import pandas as pd
import argparse
import os

def main():
    """
    Parses command-line arguments and selects variants based on allele frequency
    and an optional exclusion list.
    """
    # 1. Set up argument parsing
    parser = argparse.ArgumentParser(
        description="Selects the variant with the highest allele frequency per unique chromosome:position and excludes specified variants."
    )
    
    parser.add_argument("allele_freq_file", help="Path to the input allele frequency file (PLINK2 .afreq format).")
    parser.add_argument("outfile", help="Path to the output file (list of variant IDs to include).")
    parser.add_argument("--add_chr", default = True, help = "add 'chr' before the variant ID")
    parser.add_argument("--exclusion_list", default="_NA",
                        help="Path to a file containing variant IDs to exclude (one ID per line). Default: _NA (no exclusion)."
    )
    args = parser.parse_args()

    print(f"Reading allele frequency file: {args.allele_freq_file}")
    print(f"Writing variant inclusion list to: {args.outfile}")
    print(f"Exclusion list path: {args.exclusion_list}")

    # 2. Read in list of variants to exclude
    EXCLUSION_IDS = []
    if args.exclusion_list != "_NA":
        # Check if the exclusion file exists before attempting to open
        if not os.path.exists(args.exclusion_list):
             raise FileNotFoundError(
                 f"Exclusion list file not found at: {args.exclusion_list}"
             )
        
        with open(args.exclusion_list) as f:
            EXCLUSION_IDS = [line.strip() for line in f if line.strip()]
    
    # 3. Process the allele frequency data
    try:
        # Read the .afreq file
        afreq_df = pd.read_csv(args.allele_freq_file, sep="\t")
    except FileNotFoundError:
        # Re-raise the error with a more specific message if pandas fails
        raise FileNotFoundError(
            f"Input allele frequency file not found or inaccessible: {args.allele_freq_file}"
        )
    
    # Filter out excluded variants
    afreq_filtered = afreq_df.query("ID not in @EXCLUSION_IDS")
    # Create the 'chr_pos' column for grouping
    # ID format is expected to be CHR:POS:REF:ALT or similar, splitting on ':'
    afreq_filtered["chr_pos"] = afreq_filtered["ID"].apply(lambda x: ":".join(str(x).split(":")[0:2]))
    
    # Select the top variant (highest ALT_FREQS) for each chr:pos
    afreq_top_freq = (afreq_filtered.sort_values("ALT_FREQS", ascending=False).groupby("chr_pos", as_index=False).first())
    
    # 4. Write the output file
    if args.add_chr:
        afreq_top_freq["ID"] = afreq_top_freq["ID"].apply(lambda x : f"chr{x}") #add chr before variant id
    afreq_top_freq["ID"].to_csv(args.outfile, sep=" ", index=False, header=False)
    print(f"Created variant inclusion list with {len(afreq_top_freq)} variants.")

if __name__ == "__main__":
    main()