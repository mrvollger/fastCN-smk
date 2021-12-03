
rule copy_ref:
    input:
        ref=config["fasta"]
    output:
        ref_keep="results/{sample}/{sample}.fasta"
    log:
        "logs/{sample}/quickmer/ref.log"
    resources:
        hrs=1,
        mem=8
    threads: 1
    run:
        if input.ref.endswith(".gz"):
            shell('gunzip -c {input.ref} > {output.ref_keep}')
        else:
            shell('ln -s $(readlink -f {input.ref}) {output.ref_keep}')


rule generate_ref_file:
    input:
        ref=rules.copy_ref.output.ref_keep,
        include_bed=config.get("include_bed", rules.include_file.output.include)
    output:
        bed="results/{sample}/{sample}.fasta.bed",
        qgc="results/{sample}/{sample}.fasta.qgc",
        qm="results/{sample}/{sample}.fasta.qm"
    conda:
        "../envs/env.yml"
    log:
        "logs/{sample}/quickmer/search.log"
    params:
        kmer = config.get('quickmer_kmer_size', "30"),
        window_size = config.get('quickmer_window_size',"1000")
    threads: 16
    resources:
        hrs=12,
        mem=5,
    shell:
        """
        quicKmer2 search -k {params.kmer} -t 20 -s 3G -e 2 -d 100 -w {params.window_size} -c {input.include_bed} {input.ref}
        """


rule index_ref:
    input:
        ref="results/{sample}/{sample}.fasta"
    output:
        index="results/{sample}/{sample}.fasta.fai"
    log:
        "logs/{sample}/quickmer/fai.log"
    conda:
        "../envs/env.yml"
    resources:
        mem=12,
        hrs=2,
    threads: 1
    shell:
        """
        samtools faidx {input.ref}
        """


rule quickmer_count:
    input:
        fastq=gather.split("temp/reads/{{sm}}/{scatteritem}.fq.gz"),
        ref_qm=config.get("fasta_qm", rules.generate_ref_file.output.qm),
        ref_qgc=config.get("fasta_qgc", rules.generate_ref_file.output.qgc),
        ref_bed=config.get("fasta_bed", rules.generate_ref_file.output.bed),
        ref=config.get("quickmer_ref", rules.copy_ref.output.ref_keep)
    output:
        qm2= "temp/{sample}/sunk/{sm}/{sm}.qm2.txt",
        qm2_bin= "temp/{sample}/sunk/{sm}/{sm}.qm2.bin"
    conda:
        "../envs/env.yml"
    log:
        "logs/{sample}/quickmer/{sm}/count.log"
    resources:
        mem=12,
        hrs=12,
    threads: 8
    shell:
        """
        mkdir -p $(dirname {output.qm2})
        zcat {input.fastq} | seqtk seq -A -U -l 60 - | quicKmer2 count -t {threads} {input.ref} /dev/stdin $( echo {output.qm2_bin} | sed 's/.bin//' )
        """


rule quickmer_est:
    input:
        qm2_bin=rules.quickmer_count.output.qm2_bin,
        ref=config.get("quickmer_ref", rules.copy_ref.output.ref_keep)
    output:
        bed=temp("temp/{sample}/windows/sunk/{sm}.depth.bed.CN.bed"),
        png="temp/{sample}/sunk/{sm}/{sm}.qm2.png"
    conda:
        "../envs/env.yml"
    log:
        "logs/{sample}/quickmer/{sm}/est.log"
    resources:
        mem=25,
        hrs=4
    threads: 1
    shell:
        """
        quicKmer2 est {input.ref} $( echo {input.qm2_bin} | sed 's/.bin//' ) {output.bed}
        """

'''
rule quickmer_browser:
    input:
        qm2=rules.quickmer_est.output.bed
    output:
        browser_file = temp("temp/{sample}/windows/quickmer2/{sm}.depth.bed.CN.bed")
    conda:
        "../envs/env.yml"
    log:
        "logs/{sample}/quickmer/{sm}/bed.log"
    resources:
        mem=12,
        hrs=3
    threads: 1
    shell:
        """
        grep -v decoy {input.qm2} | grep -v chrEBV > {output.browser_file}
        """
'''
