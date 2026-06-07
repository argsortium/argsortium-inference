#!/usr/bin/env python3
"""Convert an smc++ plot CSV (time / Ne) into a Threads `.demo` file.

Threads `.demo` format: whitespace-separated rows of
    <time_in_generations>\t<2N>
where 2N is the haploid effective size (= 2 x diploid Ne), piecewise-constant,
youngest-first, starting at time 0.0 (cf. Ne10000.demo == "0.0  20000").

smc++ writes the plotted points with `smc++ plot --csv`; the CSV has at least
`x` (generations) and `y` (diploid Ne) columns. We sort by time, force the
youngest epoch to time 0.0, drop duplicate times, and emit `time  2*Ne`.

Usage: smc_demo_to_threads.py <smcpp_plot.csv> <out.demo>
"""

import csv
import sys


def load_points(csv_path):
    """Return list of (time_generations, Ne) from an smc++ plot CSV."""
    points = []
    with open(csv_path) as handle:
        reader = csv.DictReader(handle)
        columns = {name.lower(): name for name in (reader.fieldnames or [])}
        x_col = columns.get("x")
        y_col = columns.get("y")
        if x_col is None or y_col is None:
            raise SystemExit(
                f"expected 'x' and 'y' columns in {csv_path}, got {reader.fieldnames}")
        for row in reader:
            points.append((float(row[x_col]), float(row[y_col])))
    if not points:
        raise SystemExit(f"no points found in {csv_path}")
    return points


def to_demo(points):
    """Return sorted, de-duplicated (time, 2N) rows starting at time 0.0."""
    points.sort(key=lambda p: p[0])
    points[0] = (0.0, points[0][1])      # youngest epoch anchored at present
    seen_times = set()
    rows = []
    for time, ne in points:
        if time in seen_times:
            continue
        seen_times.add(time)
        rows.append((time, 2.0 * ne))    # diploid Ne -> haploid 2N
    return rows


def main():
    if len(sys.argv) != 3:
        print(__doc__)
        sys.exit(1)
    rows = to_demo(load_points(sys.argv[1]))
    with open(sys.argv[2], "w") as out:
        for time, two_n in rows:
            out.write(f"{time:.6g}\t{two_n:.6g}\n")


if __name__ == "__main__":
    main()
