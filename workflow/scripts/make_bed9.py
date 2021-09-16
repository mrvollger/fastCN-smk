import pandas as pd

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
    120: "255,192,203",
    110: "255,105,180",
    100: "199,21,133",
    90: "148,0,211",
    80: "139,0,139",
    70: "75,0,130",
    60: "0,100,0",
    50: "34,139,34",
    40: "0,255,0",
    30: "32,178,170",
    20: "0,255,255",
}

df = pd.read_csv(
    snakemake.input.cn_bed, sep="\t", header=None, names=["chr", "pos", "end", "cn"]
)
df["cn_round"] = df.apply(lambda row: int(row["cn"]), axis=1)
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
].to_csv(snakeamke.output.bed9, sep="\t", index=False, header=False)