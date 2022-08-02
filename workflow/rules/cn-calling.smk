
rule unpack_mrsfast:
    input:
        comp=rules.compress_mrsfast_further.output.comp,
    output:
        exp=pipe("results/{sample}/mapping/{sm}_merged_exp.out.gz"),
    conda:
        "../envs/env.yml"
    resources:
        mem=2,
        hrs=24,
        load=1,
    threads: 1
    benchmark:
        "benchmarks/{sample}/exp_mrsfast/{sm}.tbl"
    log:
        "logs/mrsfast/{sample}/{sm}.merged_exp.log",
    script:
        "../scripts/unpack_partial_sam.py"


rule GC_correct:
    input:
        exp=rules.unpack_mrsfast.output.exp,
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
        SAM_GC_correction \
                {input.fai} {input.bin} {input.exp} \
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
        windows=temp("temp/{sample}/windows/wssd/{sm}.depth.bed"),
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
        cn_bed=temp("temp/{sample}/windows/wssd/{sm}.depth.bed.CN.bed"),
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




