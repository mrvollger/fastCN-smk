track_db_header = """
track wssd_cn
compositeTrack off
shortLabel WSSD CN
longLabel  WSSD copy number estimates
visibility hide
priority 30
type bigBed 9 +
itemRgb on
maxItems 100000
"""

hub = """
hub WSSD_CN
shortLabel WSSD CN
longLabel WSSD copy number estimates
genomesFile genomes.txt
email mvollger.edu
"""
genomes = """
genome {sample}
trackDb trackDb.{sample}.txt
"""

track = """
    track wssd_{sm}
    parent wssd_cn
    bigDataUrl wssd/{sm}_wssd.bb
    shortLabel {sm} wssd CN
    longLabel {sm} Copy Number
    type bigBed 9 +
    itemRgb on
    visibility dense
"""


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
    threads: 1
    log:
        "logs/{sample}/tracks/trackHub.log",
    run:
        out = open(output.track, "w")
        out.write(track_db_header)
        for sm in config["reads"].keys():
            out.write(track.format(sm=sm))
        out.close()
        open(output.hub, "w").write(hub)
        open(output.genomes, "w").write(genomes.format(sample=wildcards.sample))
