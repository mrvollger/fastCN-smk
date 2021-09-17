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

with open(snakemake.output.track, "w") as out:
    out.write(track_db_header)
    [out.write(track.format(sm=sm)) for sm in snakemake.params.samples]

open(snakemake.output.hub, "w").write(hub)
open(snakemake.output.genomes, "w").write(
    genomes.format(sample=snakemake.wildcards.sample)
)
