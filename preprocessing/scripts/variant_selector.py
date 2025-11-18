import pandas as pd
import io
from typing import List, Optional

"""
This is a placeholder script to create list of variants to keep
1. Remove variants in the exclusion list given in the config
2. Group by chromosome and base pair position and only keep most common variant (i.e. closest to 0.5)
"""

def read_exclusion_list(exclusion_filepath: Optional[str]) -> List[str]:
    """
    Reads a list of variant IDs from an optional exclusion file.
    Assumes the file has one VAR_ID per line.
    Returns an empty list if the path is None or the file is not found/cannot be read.
    
    In a real-world scenario, you would replace the simulation below 
    with standard file reading using 'with open(...)'.
    """
    if not exclusion_filepath:
        print("No exclusion file path provided. Skipping exclusion.")
        return []

    try:
        # --- SIMULATED FILE READING START ---
        # Replace this simulation block with actual file I/O in production.
        if exclusion_filepath == 'exclusion_list.txt':
             # Mock content: 'rsA' and 'rsZ' are excluded
             exclusion_data = "rsA\nrsZ\n"
        else:
             print(f"Warning: Simulated exclusion file '{exclusion_filepath}' not found. Using empty list.")
             return []
        
        # Reading the simulated content
        exclusion_list = [line.strip() for line in exclusion_data.split('\n') if line.strip()]
        # --- SIMULATED FILE READING END ---

        return exclusion_list

    except Exception as e:
        print(f"Error reading exclusion file {exclusion_filepath}: {e}")
        return []

def select_highest_frequency_variants(afreq_filepath: str, exclusion_filepath: Optional[str] = None) -> list:
    """
    Reads an AFREQ file and an optional exclusion list file.
    Removes excluded variants, groups by CHROM and POS, and returns VAR_IDs
    for the variant with the highest FREQ at each remaining position.

    Args:
        afreq_filepath: Path to the main allele frequency file (TSV format).
        exclusion_filepath: Optional path to the file containing VAR_IDs to exclude (one per line).

    Returns:
        A list of VAR_IDs for the selected, non-excluded variants.
    """
    # 1. Read the main AFREQ file (Data Loading)
    try:
        # --- SIMULATED FILE READING START ---
        # In a real environment, replace this with: df = pd.read_csv(afreq_filepath, sep='\t')
        data = """
CHROM\tPOS\tVAR_ID\tA1\tFREQ
1\t1000\trsA\tA\t0.45
1\t1000\trsB\tC\t0.15
1\t1500\trsC\tG\t0.80
2\t2000\trsX\tT\t0.30
2\t2000\trsY\tA\t0.60
2\t2500\trsZ\tG\t0.50
3\t3000\trsD\tC\t0.90
"""
        df = pd.read_csv(io.StringIO(data), sep='\t')
        print(f"Successfully loaded AFREQ data from simulated file: '{afreq_filepath}'")
        # --- SIMULATED FILE READING END ---
    except Exception as e:
        print(f"Error loading AFREQ file {afreq_filepath}: {e}")
        return []

    # 2. Read the optional exclusion list
    exclusion_list = read_exclusion_list(exclusion_filepath)

    print("\n--- Processing Summary ---")
    print(f"Initial number of variants: {len(df)}")
    print(f"Variants to exclude: {exclusion_list}")

    # 3. Exclusion Step: Remove variants present in the exclusion list
    initial_len = len(df)
    df_filtered = df[~df['VAR_ID'].isin(exclusion_list)].copy()
    excluded_count = initial_len - len(df_filtered)
    print(f"Number of variants removed due to exclusion list: {excluded_count}")

    if df_filtered.empty:
        print("After exclusion, the DataFrame is empty. Returning empty list.")
        return []

    print("\n--- Filtered Data (After Exclusion, first 5 rows) ---")
    print(df_filtered.head().to_string())

    # 4. Group by Chromosome (CHROM) and Position (POS) and find the index 
    #    of the row with the maximum Allele Frequency (FREQ).
    idx_max_freq = df_filtered.groupby(['CHROM', 'POS'])['FREQ'].idxmax()

    # 5. Select the corresponding rows.
    selected_variants_df = df_filtered.loc[idx_max_freq]

    print("\n--- Selected Variants (Highest FREQ per CHROM:POS) ---")
    print(selected_variants_df.to_string())

    # 6. Extract and return the VAR_ID list.
    return selected_variants_df['VAR_ID'].tolist()

# Define simulated file paths for demonstration
AFREQ_FILE = 'data.afreq'
EXCLUSION_FILE = 'exclusion_list.txt'

# =====================================================================
print("============== TEST 1: With Optional Exclusion File ==============")
# We simulate reading the exclusion file containing ['rsA', 'rsZ']
final_var_ids_test1 = select_highest_frequency_variants(AFREQ_FILE, exclusion_filepath=EXCLUSION_FILE)

print("\n--- Final List of VAR_IDs (Test 1 Result) ---")
print(final_var_ids_test1)

# =====================================================================
print("\n\n============== TEST 2: Without Exclusion File (Optional feature check) ==============")
# We pass None, so no exclusions are performed.
final_var_ids_test2 = select_highest_frequency_variants(AFREQ_FILE, exclusion_filepath=None)

print("\n--- Final List of VAR_IDs (Test 2 Result) ---")
print(final_var_ids_test2)
