
rule GC_correct:
    input:
        merged=rules.merged_mrsfast_bam.output.merged,
        fai=config.get("masked_ref", rules.masked_reference.output.fasta) + ".fai",
        bin=config.get("gc_control", rules.fastcn_GC_bin.output.bin),
    output:
        binary=pipe("results/{sample}/binary/{sm}.bin"),
    conda:
        "../envs/env.yml"
    log:
        "logs/{sample}/binary/{sm}.log",
    resources:
        mem=8,
        hrs=24,
    threads: 1
    shell:
        """
        zcat {input.merged} \
            | SAM_GC_correction \
                {input.fai} {input.bin} /dev/stdin \
                $(dirname {output.binary})/{wildcards.sm}
        """


rule gzip_bin:
    input:
        binary=rules.GC_correct.output.binary,
    output:
        zipped="results/{sample}/binary/{sm}.bin.gz",
    conda:
        "../envs/env.yml"
    resources:
        mem=8,
        hrs=24,
    threads: 4
    log:
        "logs/{sample}/binary/{sm}.gz.log",
    shell:
        """
        pigz -p {threads} -c {input.binary} > {output.zipped}
        """


rule convert_windows:
    input:
        fai=config.get("masked_ref", rules.masked_reference.output.fasta) + ".fai",
        binary=rules.gzip_bin.output.zipped,
        ref_windows=config.get("reference_windows", rules.make_windows.output.bed),
    output:
        windows=temp("temp/{sample}/windows/{sm}.depth.bed"),
    conda:
        "../envs/env.yml"
    log:
        "logs/{sample}/windows/{sm}.log",
    resources:
        mem=16,
        hrs=24,
    threads: 1
    shell:
        """
        perbp-to-windows.py \
            --depth {input.binary} \
            --out {output.windows} \
            --chromlen {input.fai} \
            --windows {input.ref_windows}
        """


rule copy_number_call:
    input:
        windows=rules.convert_windows.output.windows,
        control_bed=config.get(
            "autosome_control", rules.autosome_control_windows.output.bed
        ),
        chrX_control_bed=config.get(
            "chrX_control", rules.chrX_control_windows.output.bed
        ),
    output:
        cn_bed=temp("temp/{sample}/windows/{sm}.depth.bed.CN.bed"),
    log:
        "logs/{sample}/windows/{sm}.cn.log",
    conda:
        "../envs/env.yml"
    resources:
        mem=16,
        hrs=24,
    threads: 1
    shell:
        """
        depth-to-cn.py --in {input.windows} --autocontrol {input.control_bed} --chrX {input.chrX_control_bed}
        """


rule bed_to_bed9:
    input:
        cn_bed=rules.copy_number_call.output.cn_bed,
    output:
        bed9="results/{sample}/tracks/bed9/{sm}.bed.gz",
    log:
        "logs/{sample}/windows/{sm}.bed9.log",
    conda:
        "../envs/env.yml"
    resources:
        mem=4,
        hrs=24,
    threads: 1
    script:
        "../scripts/make_bed9.py"


'''
rule wssd_binary:
    input:
        bed=rules.copy_number_call.output.cn_bed,
        sat_bed=SAT_BED,
        gap_bed=GAP_BED,
        cen_bed=CEN_BED,
    output:
        wssd_bin="bed/{sample}_wssd_binary.bed",
        sat_bed=temp("bed/{sample}_wssd_sat.bed"),
        temp_sat=temp("bed/{sample}_wssd_sat.bed.tmp"),
    resources:
        mem=2,
        hrs=24,
    threads: 1
    shell:
        """
        bedtools coverage -a {input.bed} -b {input.sat_bed} | cut -f 1-5,8 > {output.temp_sat}
        {SNAKEMAKE_DIR}/utils/wssd_binary.py -b {output.temp_sat} -o {output.sat_bed}
        bedtools subtract -a {output.sat_bed} -b {input.gap_bed} | bedtools subtract -a - -b {input.cen_bed} > {output.wssd_bin}
        """
'''
