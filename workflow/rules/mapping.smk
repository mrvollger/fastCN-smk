import pandas as pd
import numpy as np
import math
import os
import sys


scattergather:
    split=8,


rule split_reads:
    input:
        reads=lambda wc: config["reads"][wc.sm],
    output:
        reads=temp(scatter.split("temp/reads/{{sm}}/{scatteritem}.fq.gz")),
    resources:
        mem=4,
        hrs=8,
    params:
        sdir=SDIR,
    threads: 1
    shell:
        """
        cat {input.reads} \
            | seqtk seq -F '#' \
            | {params.sdir}/scripts/split_fastx.py --outputs {output.reads}
        """


rule mrsfast_index:
    input:
        ref=config.get("masked_ref", rules.masked_reference.output.fasta),
    output:
        index=config.get("masked_ref", rules.masked_reference.output.fasta) + ".index",
    resources:
        mem=8,
        hrs=24,
    threads: 1
    shell:
        """
        mrsfast --index {input.ref}
        """


rule mrsfast_alignment:
    input:
        reads="temp/reads/{sm}/{scatteritem}.fq.gz",
        ref=rules.mrsfast_index.output.index,
    output:
        sam=pipe("temp/mrsfast/{sample}/{sm}/{scatteritem}.sam.gz"),
    resources:
        mem=4,
        hrs=24,
    threads: 4
    shell:
        """
        extract-from-fastq36.py --in {input.reads} \
            | mrsfast --search {input.ref} --seq /dev/fd/0 \
                --disable-nohits --mem {resources.mem} --threads {threads} \
                -e 2 --outcomp \
                -o $(dirname {output.sam})/{wildcards.scatteritem}
        """


rule mrsfast_sort:
    input:
        sam=rules.mrsfast_alignment.output.sam,
    output:
        bam=temp("temp/mrsfast/{sample}/{sm}/{scatteritem}.bam"),
    resources:
        mem=4,
        hrs=24,
    threads: 4
    shell:
        """
        samtools view -b {input.sam} \
            | samtools sort -@ {threads} -o {output.bam} -
        """


rule merge_bam:
    input:
        bams=gather.split("temp/mrsfast/{{sample}}/{{sm}}/{scatteritem}.bam"),
    output:
        merged="results/{sample}/mapping/{sm}_merged.out.gz",
    resources:
        mem=8,
        hrs=24,
    threads: 4
    priority: 50
    shell:
        """
        samtools merge -@ {threads} - {input.bams} -u \
            | samtools view -h - \
            | awk 'BEGIN {{OFS="\\t"}} {{print $1,$2,$3,$4,$5}}' \
            | pigz -p {threads} \
        > {output.merged}
        """
