
module Rhodonite:
    snakefile:
        "https://github.com/mrvollger/Rhodonite/raw/master/workflow/Snakefile"
    config:
        mask_config


# import the rules from Rhodonite
use rule * from Rhodonite as Rhodonite_*


rule mask_file:
    input:
        rm=rules.Rhodonite_RepeatMasker.output.bed,
        trf=rules.Rhodonite_trf.output.bed,
        gaps=rules.Rhodonite_gaps.output.bed,
        fai=f'{config["fasta"]}.fai',
    output:
        bed="results/{sample}.mask.bed",
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
        bed="results/{sample}.exclude.bed",
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
    params:
        window=400,
    shell:
        """
        export PATH=$PWD/fastCN:$PATH
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
        fa="results/{sample}/{sample}.masked.fasta",
        fai="results/{sample}/{sample}.masked.fasta.fai",
    conda:
        "envs/env.yml"
    log:
        "logs/{sample}/mask_reference.log",
    shell:
        """
        seqtk seq -M {input.mask} -n N \
            {input.fasta} -l 60 \
            > {output.fa}
        samtools faidx {output.fa}
        """


rule make_windows:
    input:
        fa=rules.masked_reference.output.fa,
    output:
        bed="results/{sample}/{sample}.windows.bed",
    threads: 1
    conda:
        "envs/env.yml"
    log:
        "logs/{sample}/make_windows.log",
    params:
        window=config.get("window", 1000),
    script:
        "scripts/make_windows.py"


rule autosome_control_windows:
    input:
        mask=rules.mask_file.output.bed,
        exclude=rules.exclude_file.output.bed,
        windows=rules.make_windows.output.bed,
    output:
        bed="results/{sample}/{sample}_auto_control.bed",
    resources:
        mem=8,
        hrs=24,
    threads: 1
    shell:
        """
        cat {input.mask} {input.exclude} \
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
        cat {input.mask} {input.exclude} \
            | bedtools sort -i - \
            | bedtools merge -i - \
            | bedtools subtract -A -a {input.windows} -b - \
            | grep -w chrX \
            > {output.bed}
        """
