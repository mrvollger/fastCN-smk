import pandas as pd
import math

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

df = pd.read_csv(
    snakemake.input.cn_bed, sep="\t", header=None, names=["chr", "pos", "end", "cn"]
)
df["cn_round"] = df["cn"].round(0).astype("int32")
# df.apply(lambda row: int(row["cn"]), axis=1)
df["mapq"] = "0"
df["strand"] = "+"
df["big_start"] = "0"
df["big_end"] = "0"
df["color"] = df.apply(
    lambda row: color_hash[max(0, int(round(row["cn"])))]
    if row["cn"] < 10
    else color_hash[min(120, int(math.floor(round(row["cn"]) / 10.0) * 10))],
    axis=1,
)
df.sort_values(["chr", "pos"], inplace=True)
df[
    [
        "chr",
        "pos",
        "end",
        "cn_round",
        "mapq",
        "strand",
        "big_start",
        "big_end",
        "color",
        "cn",
    ]
].to_csv(snakemake.output.bed9, sep="\t", index=False, header=False, compression="gzip")
