#!/usr/bin/env python3

import argparse
import collections
import pandas as pd
import networkx as nx
import numpy as np
import tskit
import os


'''
Currently this function doesn't work
alternative approach take bed file and create a list of local trees using tskit function .from_newick
'''

def convert_argweaver(infile, filename=None):
    """
    Convert an ARGweaver .arg file to a tskit tree sequence.
    Function from Yan Wong: https://github.com/tskit-dev/what-is-an-arg-paper/blob/3d3b4309182a6f85ec98c74f23aa4a72bbb0dda9/argutils/convert.py#L12
    Parameters
    ----------
    infile : file handle
        Opened .arg file (text mode).
    """
    start, end = next(infile).strip().split()
    assert start.startswith("start=")
    start = int(start[len("start="):])
    assert end.startswith("end=")
    end = int(end[len("end="):])

    df = pd.read_csv(infile, header=0, sep="\t", dtype={"name": str, "parents": str})
    for col in ("name", "parents", "age"):
        if col not in df.columns:
            raise ValueError(f"Column {col} not found in ARGweaver file")

    name_to_record = {str(row["name"]): dict(row) for _, row in df.iterrows()}

    parent_map = collections.defaultdict(list)
    G = nx.DiGraph()
    time_map = {}

    for row in name_to_record.values():
        child = row["name"]
        parents = row["parents"]
        time_map[row["age"]] = row["age"]
        G.add_node(child)
        if isinstance(parents, str):
            for parent in parents.split(","):
                G.add_edge(child, parent)
                parent_map[child].append(parent)

    tables = tskit.TableCollection(sequence_length=end)
    tables.nodes.metadata_schema = tskit.MetadataSchema.permissive_json()
    breakpoints = np.full(len(G), tables.sequence_length)
    aw_to_tsk_id = {}

    min_time_diff = min(np.diff(sorted(time_map.keys())))
    epsilon = min_time_diff / 1e6

    try:
        for node in nx.lexicographical_topological_sort(G):
            record = name_to_record[node]
            flags = 0
            if record["event"] == "gene":
                flags = tskit.NODE_IS_SAMPLE
                assert record["age"] == 0
                time = record["age"]
            else:
                if record["age"] == 0:
                    time_map[record["age"]] += epsilon
                time = time_map[record["age"]]
                time_map[record["age"]] += epsilon
            tsk_id = tables.nodes.add_row(flags=flags, time=time, metadata=record)
            aw_to_tsk_id[node] = tsk_id
            if record["event"] == "recomb":
                breakpoints[tsk_id] = record["pos"]
    except nx.exception.NetworkXUnfeasible:
        bad_edges = nx.find_cycle(G, orientation="original")
        raise nx.exception.NetworkXUnfeasible(
            f"Cycle found in ARGweaver graph ({filename}): {bad_edges}"
        )

    L = tables.sequence_length
    for aw_node in G:
        child = aw_to_tsk_id[aw_node]
        parents = [aw_to_tsk_id[p] for p in parent_map[aw_node]]
        if len(parents) == 1:
            tables.edges.add_row(0, L, parents[0], child)
        elif len(parents) == 2:
            x = breakpoints[child]
            tables.edges.add_row(0, x, parents[0], child)
            tables.edges.add_row(x, L, parents[1], child)
        else:
            assert len(parents) == 0

    tables.sort()
    ts = tables.tree_sequence()
    return ts.simplify(keep_unary=True)


def main():
    parser = argparse.ArgumentParser(
        description="Convert ARGweaver .arg file to .tskit.trees format"
    )
    parser.add_argument(
        "arg_file",
        help="Input ARGweaver .arg file"
    )
    parser.add_argument(
        "-o", "--output",
        help="Output filename (defaults to replacing .arg with .tskit.trees)"
    )

    args = parser.parse_args()

    if args.output:
        output_file = args.output
    else:
        if not args.arg_file.endswith(".arg"):
            raise ValueError("Input file must end with .arg or specify --output")
        output_file = args.arg_file[:-4] + ".tskit.trees"

    with open(args.arg_file, "r") as infile:
        ts = convert_argweaver(infile, filename=os.path.basename(args.arg_file))

    ts.dump(output_file)
    print(f"Tree sequence written to {output_file}")


if __name__ == "__main__":
    main()
