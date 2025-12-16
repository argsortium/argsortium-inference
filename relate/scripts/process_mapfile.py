#!/usr/bin/env python3
"""
Convert HapMap genetic map format to Relate genetic map format.

HapMap format example:
Chromosome  Position(bp)  Rate(cM/Mb)  Map(cM)
chr17       550394445     0.004866000  0.000000000

Relate format example:
pos COMBINED_rate Genetic_Map
27572179 1.36 0.0
28322179 1.36 1.02
"""

import argparse
import sys


def convert_hapmap_to_relate(input_file, pos_col, rate_col, map_col, output_file=None):
    """
    Convert HapMap format to Relate format.
    
    Args:
        input_file: Path to input HapMap file
        pos_col: Column number for position (0-indexed)
        rate_col: Column number for rate (0-indexed)
        map_col: Column number for genetic map (0-indexed)
        output_file: Path to output file (optional, defaults to input + .processed.relate)
    """
    
    # If no output file specified, add .processed.relate to input filename
    if output_file is None:
        output_file = input_file + ".processed.relate"
    elif not output_file.endswith('.processed.relate'):
        # Ensure output file has correct suffix
        if output_file.endswith('.relate'):
            output_file = output_file.replace('.relate', '.processed.relate')
        else:
            output_file = output_file + '.processed.relate'
    
    with open(input_file, 'r') as infile:
        # Skip header
        header = infile.readline().strip()
        
        # Write output file
        with open(output_file, 'w') as outfile:
            # Write Relate format header
            outfile.write("pos COMBINED_rate Genetic_Map\n")
            
            # Process each line
            for line in infile:
                line = line.strip()
                if not line:
                    continue
                    
                cols = line.split()
                
                # Check if all required columns exist
                max_col = max(pos_col, rate_col, map_col)
                if len(cols) <= max_col:
                    print(f"Warning: Skipping line with insufficient columns: {line}", file=sys.stderr)
                    continue
                
                pos = cols[pos_col]
                rate = cols[rate_col]
                genetic_map = cols[map_col]
                
                outfile.write(f"{pos} {rate} {genetic_map}\n")
    
    print(f"Conversion complete. Output written to: {output_file}")


def main():
    parser = argparse.ArgumentParser(
        description='Convert HapMap genetic map format to Relate genetic map format',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Example usage:
  # Using default column numbers (pos=1, rate=2, map=3) and auto-generated output name
  python process_mapfile.py input.hapmap
  
  # Specifying custom output file
  python process_mapfile.py input.hapmap -o output.txt
  
  # Specifying custom column numbers (0-indexed)
  python process_mapfile.py input.hapmap --pos 1 --rate 2 --map 3
        """
    )
    
    parser.add_argument(
        'input',
        help='Input HapMap genetic map file'
    )
    
    parser.add_argument(
        '-o', '--output',
        default=None,
        help='Output file path (optional, defaults to input + .processed.relate)'
    )
    
    parser.add_argument(
        '--pos',
        type=int,
        default=1,
        help='Column number for position, 0-indexed (default: 1)'
    )
    
    parser.add_argument(
        '--rate',
        type=int,
        default=2,
        help='Column number for recombination rate, 0-indexed (default: 2)'
    )
    
    parser.add_argument(
        '--map',
        type=int,
        default=3,
        help='Column number for genetic map distance, 0-indexed (default: 3)'
    )
    
    args = parser.parse_args()
    
    try:
        convert_hapmap_to_relate(
            args.input,
            args.pos,
            args.rate,
            args.map,
            args.output
        )
    except FileNotFoundError:
        print(f"Error: Input file '{args.input}' not found", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == '__main__':
    main()
