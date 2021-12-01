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
        mem=4,
        hrs=24,
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
        mem=2,
        hrs=24,
    threads: 1
    log:
        "logs/{sample}/tracks/{type}/{sm}.bigbed.log",
    params:
        as_file=f"{SDIR}/utils/track.as",
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
            rules.make_bb.output, sm=config["reads"].keys(), allow_missing=True
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
        mem=2,
        hrs=24,
    log:
        "logs/{sample}/tracks/{type}/trackHub.log",
    params:
        samples=list(config["reads"].keys()),
        reads=list(config["reads"].values())
    script:
        "../scripts/make_trackdb.py"
