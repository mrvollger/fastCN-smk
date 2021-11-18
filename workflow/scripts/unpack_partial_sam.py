#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Author: William T. Harvey


import gzip

file = gzip.open(snakemake.input.comp, 'rt')

# Consume header
line = file.readline()

with open(snakemake.output.exp, 'w') as out_file:
	out_file.write('@\n')
	while True:
		line = file.readline().rstrip()
		if not line:
			break
		else:
			line = line.split(',')
			out_file.write('\n'.join(['\t'.join(line)] * int(line[0])))
			out_file.write('\n')

file.close()


