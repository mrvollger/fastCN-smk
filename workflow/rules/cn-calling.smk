
rule GC_correct:
    input:
        merged="mappings/{sm}_bed.gz",
        fai=config.get("masked_ref", rules.masked_reference.output.fai),
        gc_control=config.get("gc_control", rules.fastcn_GC_bin.output.bin),
    output:
        binary=temp("results/binary/{sample}/{sm}.bin"),
    resources:
        mem=8,
        hrs=24,
    threads: 1
    shell:
        """
        zcat {input.merged} \
            | SAM_GC_correction \
                {input.fai} {input.bin} /dev/stdin \
                results/binary/{wildcards.sample}/{wildcards.sm}
        """


rule gzip_bin:
    input:
        binary=rules.GC_correct.output.binary,
    output:
        zipped="results/binary/{sample}/{sm}.bin.gz",
    resources:
        mem=8,
        hrs=24,
    threads: 4
    shell:
        """
        pigz -p {threads} -c {input.binary} > {output.zipped}
        """


rule convert_windows:
    input:
        fai=config.get("masked_ref", rules.masked_reference.output.fai),
        binary=rules.gzip_bin.output.zipped,
        ref_windows=config.get("masked_ref", rules.make_windows.output.bed),
    output:
        windows="results/windows/{sample}/{sm}.depth.1kb.bed",
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
        cn_bed="windows/{sample}/{sm}.depth.1kb.bed.CN.bed",
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
        bed9=temp("results/tracks/{sample}/{sm}.bed9"),
    resources:
        mem=4,
        hrs=24,
    threads: 1
    run:
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
            input.cn_bed, sep="\t", header=None, names=["chr", "pos", "end", "cn"]
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
        ].to_csv(output.bed9, sep="\t", index=False, header=False)


rule make_bb:
    input:
        bed=rules.bed_to_bed9.output.bed9,
        fai=config.get("masked_ref", rules.masked_reference.output.fasta) + ".fai",
    output:
        track="results/tracks/{sample}/{sm}_wssd.bb",
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
            {input.bed} {input.fai} {output.track}
        """


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
