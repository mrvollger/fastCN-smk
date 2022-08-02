#!/bin/env python

import numpy as np
import pandas as pd
import argparse

parser = argparse.ArgumentParser()

parser.add_argument(
    "--bed",
    "-b",
    type=str,
    required=True,
    help="Bed file with overlap of satellite regions with columns: (['chr', 'start','end', 'sample', 'cn', 'sat'])",
)
parser.add_argument("--output", "-o", type=str, required=True, help="Output file")

args = parser.parse_args()

wssd_df = pd.read_csv(
    args.bed,
    sep="\t",
    header=None,
    names=["chr", "start", "end", "sample", "cn", "sat"],
)


wssd_df["len"] = wssd_df["end"] - wssd_df["start"]
# min_len = min(wssd_df['len'])
min_len = 500


chrom_dict = {}
for chrom in wssd_df["chr"].unique():
    chrom_dict[chrom] = []
    chrom_df = wssd_df[wssd_df["chr"] == chrom].reset_index(drop=True)
    dup = [
        -7
        if (chrom_df.at[index, "len"] > 1000000)
        else 1
        if ((chrom_df.at[index, "sat"] <= 0.7) and (chrom_df.at[index, "cn"] >= 3.5))
        else 0
        for index in chrom_df.index
    ]
    start = False
    for i in range(len(dup) - 7):
        if (np.sum(dup[i : i + 7]) >= 6) and (start == True):
            window_end = chrom_df.at[i + 6, "end"]
        elif (
            (np.sum(dup[i : i + 7]) >= 6)
            and (start == False)
            and (min_len / chrom_df.at[i, "len"] > 0.2)
        ):
            window_start = chrom_df.at[i, "start"]
            window_end = chrom_df.at[i + 6, "end"]
            start = True
        elif (np.sum(dup[i : i + 7]) < 6) and (start == True):
            chrom_dict[chrom].append(
                "\t".join([str(chrom), str(window_start), str(window_end)])
            )
            start = False
        elif (np.sum(dup[i : i + 7]) >= 6) and (start == True):
            chrom_dict[chrom].append(
                "\t".join([str(chrom), str(window_start), str(window_end)])
            )
            start = False
        else:
            continue


with open(args.output, "w") as outFile:
    for contig in chrom_dict:
        if len(chrom_dict[contig]) != 0:
            outFile.write("%s\n" % ("\n".join(chrom_dict[contig])))
