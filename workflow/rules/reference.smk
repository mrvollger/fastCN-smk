mask_config = {"samples": {config["reference"]: config["fasta"]}}


wildcard_constraints:
    sample=config["reference"],


module Rhodonite:
    snakefile:
        "https://github.com/mrvollger/Rhodonite/raw/master/workflow/Snakefile"
    config:
        mask_config


# import the rules from Rhodonite
use rule * from Rhodonite as Rhodonite_*


shell.prefix(f"export PATH=$PWD/fastCN:$PATH; set -eo pipefail;")


rule mask_file:
    input:
        rm=rules.Rhodonite_RepeatMasker.output.bed,
        trf=rules.Rhodonite_trf.output.bed,
        gaps=rules.Rhodonite_gaps.output.bed,
        fai=f'{config["fasta"]}.fai',
    output:
        bed="results/{sample}/{sample}.mask.bed",
    conda:
        "../envs/env.yml"
    shell:
        """
        zcat -f -- {input.rm} {input.trf} {input.gaps} \
            | cut -f 1-3 \
            | bedtools sort -i - -g {input.fai} \
            | bedtools slop -i - -g {input.fai} -b 36 \
            | bedtools merge -i - \
        > {output.bed}
        """


rule exclude_file:
    input:
        sd=config["sd"],
        wm=rules.Rhodonite_windowmasker.output.bed,
        fai=f'{config["fasta"]}.fai',
    output:
        bed="results/{sample}/{sample}.exclude.bed",
    conda:
        "../envs/env.yml"
    params:
        window=400,
    shell:
        """
        cat \
            <(bedtools slop -b 36 -g {input.fai} -i {input.wm}) \
            <(bedtools slop -b {params.window} -g {input.fai} -i {input.sd}) \
            | cut -f 1-3 \
            | bedtools sort -i - -g {input.fai} \
            | bedtools merge -i - \
        > {output.bed}
        """


rule fastcn_GC_bin:
    input:
        mask=rules.mask_file.output.bed,
        exclude=rules.exclude_file.output.bed,
        fasta=config["fasta"],
        fai=f'{config["fasta"]}.fai',
    output:
        bin="results/{sample}/{sample}.bin",
    log:
        "logs/{sample}.GC_mask.log",
    conda:
        "../envs/env.yml"
    params:
        window=400,
    shell:
        """
        GC_control_gen \
            {input.fasta} \
            {input.exclude} \
            {input.mask} \
            {params.window} \
            {output.bin} \
            > {log} 2>&1
        """


rule masked_reference:
    input:
        mask=rules.mask_file.output.bed,
        fasta=config["fasta"],
        fai=f'{config["fasta"]}.fai',
    output:
        fasta="results/{sample}/{sample}.masked.fasta",
        fai="results/{sample}/{sample}.masked.fasta.fai",
    conda:
        "../envs/env.yml"
    log:
        "logs/{sample}/mask_reference.log",
    shell:
        """
        seqtk seq -M {input.mask} -n N \
            {input.fasta} -l 60 \
            > {output.fasta}
        samtools faidx {output.fasta}
        """


rule make_windows:
    input:
        fasta=rules.masked_reference.output.fasta,
    output:
        bed="results/{sample}/{sample}.windows.bed",
    threads: 1
    conda:
        "../envs/env.yml"
    log:
        "logs/{sample}/make_windows.log",
    params:
        window=config.get("window", 1000),
    script:
        "../scripts/windows.py"


rule autosome_control_windows:
    input:
        mask=rules.mask_file.output.bed,
        exclude=rules.exclude_file.output.bed,
        windows=rules.make_windows.output.bed,
    output:
        bed="results/{sample}/{sample}_auto_control.bed",
    conda:
        "../envs/env.yml"
    resources:
        mem=8,
        hrs=24,
    threads: 1
    shell:
        """
        less {input.mask} {input.exclude} \
            | bedtools sort -i - \
            | bedtools merge -i - \
            | bedtools subtract -A -a {input.windows} -b - \
            | grep -v chrX \
            | grep -v chrY \
            > {output.bed}
        """


rule chrX_control_windows:
    input:
        mask=rules.mask_file.output.bed,
        exclude=rules.exclude_file.output.bed,
        windows=rules.make_windows.output.bed,
    output:
        bed="results/{sample}/{sample}_chrX_control.bed",
    resources:
        mem=8,
        hrs=24,
    threads: 1
    shell:
        """
        (less {input.mask} {input.exclude} \
            | bedtools sort -i - \
            | bedtools merge -i - \
            | bedtools subtract -A -a {input.windows} -b - \
            | grep -w chrX  || true ) \
            > {output.bed}
        """
