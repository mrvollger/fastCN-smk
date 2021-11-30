#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Author: William T. Harvey

import gzip

file = gzip.open(snakemake.input.sam, "rt")

count = 0

line_cols = []


with gzip.open(snakemake.output.comp, "wt") as outfile:
    outfile.write("qname,flag,rname,pos,mapq\n")
    while True:
        line = file.readline().rstrip()
        if not line:
            break
        if line[0] == "@":
            continue
        if line.split("\t")[1:] == line_cols:
            count += 1
        elif line_cols == []:
            line_cols = line.split("\t")[1:]
            count = 1
        else:
            out_string = ",".join([str(count)] + line_cols)
            outfile.write(f"{out_string}\n")
            line_cols = line.split("\t")[1:]
            count = 1

file.close()
