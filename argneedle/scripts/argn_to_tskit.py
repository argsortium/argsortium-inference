#!/usr/bin/env python3

import argparse
import os
import arg_needle_lib

def main():
    parser = argparse.ArgumentParser(description="Convert an ArgNeedle .argn file to a tskit .trees file")
    parser.add_argument("argn_file", help="Path to the .argn file to convert")
    args = parser.parse_args()

    input_path = args.argn_file
    if not input_path.endswith(".argn"):
        raise ValueError("Input file must end with .argn")

    # Create output path by replacing .argn with .trees
    output_path = os.path.splitext(input_path)[0] + ".trees"

    # Deserialize and convert
    myarg = arg_needle_lib.deserialize_arg(input_path)
    myarg_tskit = arg_needle_lib.arg_to_tskit(myarg)

    # Dump to .trees
    myarg_tskit.dump(output_path)
    print(f"Converted {input_path} → {output_path}")

if __name__ == "__main__":
    main()

