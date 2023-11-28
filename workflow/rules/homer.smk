rule homer_annotate_peaks:
    input:
        genome="reference/{species}.{build}.{release}.dna.fasta",
        genome_index="reference/{species}.{build}.{release}.dna.fasta.fai",
        gtf="reference/{species}.{build}.{release}.gtf",
        wig="results/Coverage/{species}.{build}.{release}.{datatype}/{sample}.bw",
        peak="tmp/xsv_select/{species}.{build}.{release}.dna/{macs2_preak_type}/{sample}_peaks.{macs2_preak_type}.bed",
    output:
        annotations=protected(
            "results/{species}.{build}.{release}.dna/PeakCalling/{sample}.{macs2_preak_type}.bed"
        ),
    log:
        "logs/homer/annotatepeaks/{species}.{build}.{release}.dna/{sample}.{macs2_preak_type}.log",
    benchmark:
        "benchmark/homer/annotatepeaks/{species}.{build}.{release}.dna/{sample}.{macs2_preak_type}.tsv",
    params:
        mode="tss",
        extra=config.get("params", {}).get("homer", {}).get("annotatepeaks", "")
    wrapper:
        f"{snakemake_wrappers_version}/bio/homer/annotatePeaks"
    