rule bed_to_bed9:
    input:
        cn_bed="temp/{sample}/windows/{type}/{sm}.depth.bed.CN.bed",
    output:
        bed9="results/{sample}/tracks/bed9/{type}/{sm}.bed.gz",
    log:
        "logs/{sample}/windows/{sm}.{type}.bed9.log",
    conda:
        "../envs/env.yml"
    resources:
        mem_mb=1024 * 8,
        runtime=60 * 24,
    threads: 1
    script:
        "../scripts/make_bed9.py"


rule make_bb:
    input:
        bed=rules.bed_to_bed9.output.bed9,
        fai=config.get("masked_ref", rules.masked_reference.output.fasta) + ".fai",
    output:
        bed=temp("temp/{sample}/tracks/{type}/{sm}.bed"),
        bigbed="results/{sample}/tracks/{type}/bigbed/{sm}_{type}.bb",
    conda:
        "../envs/env.yml"
    resources:
        mem_mb=1024 * 2,
        runtime=60 * 24,
    threads: 1
    log:
        "logs/{sample}/tracks/{type}/{sm}.bigbed.log",
    params:
        as_file=workflow.source_path("../utils/track.as"),
    shell:
        """
        zcat {input.bed} > {output.bed}
        bedToBigBed -tab -type=bed9+1 \
            -as={params.as_file} \
            {output.bed} {input.fai} {output.bigbed}
        """


rule make_trackdb:
    input:
        bigwig=expand(
            rules.make_bb.output.bigbed, sm=config["reads"].keys(), allow_missing=True
        ),
    output:
        track="results/{sample}/tracks/{type}/trackDb.{sample}.txt",
        hub="results/{sample}/tracks/{type}/hub.txt",
        genomes="results/{sample}/tracks/{type}/genomes.txt",
        html="results/{sample}/tracks/{type}/bigbed/description.html",
    conda:
        "../envs/env.yml"
    threads: 1
    resources:
        mem_mb=1024 * 2,
        runtime=60 * 24,
    log:
        "logs/{sample}/tracks/{type}/trackHub.log",
    params:
        samples=list(config["reads"].keys()),
        reads=list(config["reads"].values()),
    script:
        "../scripts/make_trackdb.py"


rule wssd_binary:
    input:
        bed="results/{sample}/tracks/bed9/wssd/{sm}.bed.gz",
        sat_bed=SAT_BED,
        gap_bed=GAP_BED,
        cen_bed=CEN_BED,
    output:
        sat_bed=temp("results/{sample}/wssd/{sm}_wssd_sat.bed"),
        temp_sat=temp("results/{sample}/wssd/{sm}_wssd_sat.bed.tmp"),
        wssd_bin="results/{sample}/wssd/{sm}_wssd_binary.bed",
    params:
        script=workflow.source_path("../scripts/wssd_binary.py"),
    conda:
        "../envs/env.yml"
    log:
        "logs/{sample}/wssd/{sm}_binary.log",
    resources:
        mem_mb=1024 * 2,
        runtime=60 * 24,
    threads: 1
    shell:
        """
        bedtools coverage -a {input.bed} -b {input.sat_bed} | cut -f 1-4,10,14 > {output.temp_sat}
        python {params.script} -b {output.temp_sat} -o {output.sat_bed}
        bedtools subtract -a {output.sat_bed} -b {input.gap_bed} | bedtools subtract -a - -b {input.cen_bed} > {output.wssd_bin}
        """
