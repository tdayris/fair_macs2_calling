rule homer_annotate_peaks:
    input:
        genome="reference/{species}.{build}.{release}.{datatype}.fasta",
        genome_index="reference/{species}.{build}.{release}.{datatype}.fasta.fai",
        gtf="reference/{species}.{build}.{release}.gtf",
        wig="results/{species}.{build}.{release}.{datatype}/Coverage/{sample}.bw",
        peak="tmp/xsv_select/{species}.{build}.{release}.{datatype}/{macs2_preak_type}/{sample}_peaks.{macs2_preak_type}.bed",
    output:
        annotations=protected(
            "results/{species}.{build}.{release}.{datatype}/PeakCalling/{sample}.{macs2_preak_type}.bed"
        ),
    log:
        "logs/homer/annotatepeaks/{species}.{build}.{release}.{datatype}/{sample}.{macs2_preak_type}.log",
    benchmark:
        "benchmark/homer/annotatepeaks/{species}.{build}.{release}.{datatype}/{sample}.{macs2_preak_type}.tsv"
    params:
        mode="tss",
        extra=config.get("params", {}).get("homer", {}).get("annotatepeaks", ""),
    wrapper:
        f"{snakemake_wrappers_version}/bio/homer/annotatePeaks"
