rule make_bb:
    input:
        bed=rules.bed_to_bed9.output.bed9,
        fai=config.get("masked_ref", rules.masked_reference.output.fasta) + ".fai",
    output:
        bigbed="results/{sample}/tracks/wssd/{sm}_wssd.bb",
    conda:
        "../envs/env.yml"
    resources:
        mem=2,
        hrs=24,
    threads: 1
    log:
        "logs/{sample}/tracks/{sm}.bigbed.log",
    params:
        as_file=f"{SDIR}/utils/track.as",
    shell:
        """
        bedToBigBed -tab -type=bed9+1 \
            -as={params.as_file} \
            {input.bed} {input.fai} {output.bigbed}
        """


rule make_trackdb:
    input:
        bigwig=expand(
            rules.make_bb.output, sm=config["reads"].keys(), allow_missing=True
        ),
    output:
        track="results/{sample}/tracks/trackDb.{sample}.txt",
        hub="results/{sample}/tracks/hub.txt",
        genomes="results/{sample}/tracks/genomes.txt",
    conda:
        "../envs/env.yml"
    threads: 1
    log:
        "logs/{sample}/tracks/trackHub.log",
    params:
        samples=list(config["reads"].keys()),
    script:
        "../scripts/make_trackdb.py"
