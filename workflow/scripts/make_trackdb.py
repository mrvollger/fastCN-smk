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
html wssd/description.html
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

html = """
<html>
<h3>Description</h3>
This track represents copy number estimates. Copy number is estimated over 500 bp windows of uniquely mappable sequence. Sequences are colored from cold to hot (0 - 120+) and exact copy can be found by clicking on the region of interest. \n\n 

<h3>Code Availability</h3>
<a href="https://github.com/mrvollger/fastCN-smk">GitHub</a>

<h3>Copy Number Key</h3>
{cn_key}

<h3>Credits</h3>
Please feel free to contact <a href=mailto:wharvey@uw.edu >William Harvey</a> or <a href=mailto:mvollger@uw.edu>Mitchell Vollger</a> with any questions and/or concerns regarding this track.

<h3>References</h3>
Bailey JA, Gu Z, Clark RA, Reinert K, Samonte RV, Schwartz S, Adams MD, Myers EW, Li PW, Eichler EE. Recent segmental duplications in the human genome. Science 2002
<br><br>

Pendleton AL, Shen F, Taravella AM, Emery S, Veeramah KR, Boyko AR, Kidd JM. Comparison of village dog and wolf genomes highlights the role of the neural crest in dog domestication. BMC Biol. 2018 
<br><br>

Sudmant PH, Mallick S, Nelson BJ, Hormozdiari F, Krumm N, Huddleston J, et al. Global diversity, population stratification, and selection of human copy-number variation. Science. 2015
<br><br>

Sudmant PH, Kitzman JO, Antonacci F, Alkan C, Malig M, Tsalenko A, et al. Diversity of human copy number. Science. 2010

<h3>Sample table</h3>
{sample_table}

<br><br>
</html>
"""

color_hash = {
    0: "229,229,229",
    1: "196,196,196",
    2: "0,0,0",
    3: "0,0,205",
    4: "65,105,225",
    5: "100,149,237",
    6: "180,238,180",
    7: "255,255,0",
    8: "255,165,0",
    9: "139,26,26",
    10: "255,0,0",
    20: "0,255,255",
    30: "32,178,170",
    40: "0,255,0",
    50: "34,139,34",
    60: "0,100,0",
    70: "75,0,130",
    80: "139,0,139",
    90: "148,0,211",
    100: "199,21,133",
    110: "255,105,180",
    120: "255,192,203",
}


def html_table(color_dict, h1, h2, rgb=True):
    rtn = '<table border="1">'
    rtn += f"<tr><th>{h1}</th><th>{h2}</th></tr>"
    for key, value in color_dict.items():
        if rgb:
            second = f'<div style="font-size:15px; color:rgb({value})">&#9632;</div>'
            # bgcolor= td
        else:
            second = value
        rtn += "<tr>"
        rtn += f'<td style="text-align: center; vertical-align: middle;">{key}</td>'
        rtn += f'<td style="background: rgb({value}); text-align: center; vertical-align: middle;">{second}</td>'
        rtn += "</tr>"
    rtn += "</table>"
    return rtn


with open(snakemake.output.track, "w") as out:
    out.write(track_db_header)
    [out.write(track.format(sm=sm)) for sm in snakemake.params.samples]

open(snakemake.output.hub, "w").write(hub)
open(snakemake.output.genomes, "w").write(
    genomes.format(sample=snakemake.wildcards.sample)
)

sample_table = dict(zip(snakemake.params.samples, snakemake.params.reads))

open(snakemake.output.html, "w").write(
    html.format(
        cn_key=html_table(color_hash, h1="Copy number", h2="Color"),
        sample_table=html_table(sample_table, h1="Sample", h2="Read file", rgb=False),
    )
)
