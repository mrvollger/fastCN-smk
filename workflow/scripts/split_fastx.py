#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Author: Mitchell R. Vollger
import argparse
import pysam
import gzip

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="", formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--infile", help="input fastq file", default="/dev/stdin")
    parser.add_argument(
        "--outputs", nargs="+", help="list of output files will be compressed output"
    )
    args = parser.parse_args()
    N_IDS = len(args.outputs)

    outs = [gzip.open(f, "wb") for f in args.outputs]
    # outs = [open(f, "w+") for f in args.outputs]
    out_idx = 0
    for rec in pysam.FastxFile(args.infile, persist=False):
        outs[out_idx].write((str(rec) + "\n").encode())
        # outs[out_idx].write((str(rec) + "\n"))

        out_idx += 1
        if out_idx == N_IDS:
            out_idx = 0

    for out in outs:
        out.close()
