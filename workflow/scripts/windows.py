#!/usr/bin/env python
import argparse
import os
import sys
import pysam
from numba import njit
from functools import partial

# import multiprocessing


@njit
def make_dynamic_windows(name, seq, window):
    window_count = 0
    window_start = 0
    window_end = 0
    out = []
    for char in seq:
        if char.upper() != "N":
            window_count += 1
        window_end += 1

        if window_count == window:
            out.append(
                name
                + "\t"
                + str(window_start)
                + "\t"
                + str(window_end)
                + "\t"
                + str(window_count)
                + "\n"
            )
            window_start = window_end
            window_count = 0
            # print("here")
    out.append(
        name
        + "\t"
        + str(window_start)
        + "\t"
        + str(window_end)
        + "\t"
        + str(window_count)
        + "\n"
    )
    return out


# global var for inputs
args = None

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="", formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(
        "--infile", help="positional input", default=snakemake.input.fasta
    )
    parser.add_argument(
        "--window",
        help="length of window in non masked bases",
        type=int,
        default=snakemake.params.window,
    )
    parser.add_argument(
        "--outfile", help="positional output bed", default=snakemake.output.bed
    )
    args = parser.parse_args()
    out = open(args.outfile, "a")
    for rec in pysam.FastxFile(args.infile):
        rtn = make_dynamic_windows(rec.name, rec.sequence, args.window)
        for line in rtn:
            out.write(line)
    out.close()
