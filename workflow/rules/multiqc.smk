rule fair_macs2_calling_multiqc_config:
    input:
        "tmp/fair_fastqc_multiqc_bigr_logo.png",
    output:
        temp("tmp/fair_macs2_calling_multiqc_config.yaml"),
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 512,
        runtime=lambda wildcards, attempt: attempt * 5,
        tmpdir=tmp,
    localrule: True
    log:
        "logs/fair_macs2_calling_multiqc_config.log",
    benchmark:
        "benchmark/fair_macs2_calling_multiqc_config.tsv"
    params:
        extra=lookup_config(
            dpath="params/fair_macs2_calling_multiqc_config", default=None
        ),
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/fair_macs2_calling_multiqc_config.py"


rule fair_macs2_calling_multiqc_report:
    input:
        config="tmp/fair_macs2_calling_multiqc_config.yaml",
        picard_qc=collect(
            "tmp/fair_bowtie2_mapping_picard_create_multiple_metrics/{sample.species}.{sample.build}.{sample.release}.{datatype}/stats/{sample.sample_id}{ext}",
            sample=lookup(
                query="species == '{species}' & release == '{release}' & build == '{build}'",
                within=samples,
            ),
            ext=[
                ".alignment_summary_metrics",
                ".insert_size_metrics",
                ".insert_size_histogram.pdf",
                ".base_distribution_by_cycle_metrics",
                ".base_distribution_by_cycle.pdf",
                ".gc_bias.detail_metrics",
                ".gc_bias.summary_metrics",
                ".gc_bias.pdf",
            ],
            datatype="{datatype}",
        ),
        fastp_pair_ended=collect(
            "tmp/fair_bowtie2_mapping_fastp_trimming_pair_ended/{sample.sample_id}.fastp.json",
            sample=lookup(
                query="downstream_file == downstream_file & species == '{species}' & release == '{release}' & build == '{build}'",
                within=samples,
            ),
        ),
        fastp_single_ended=collect(
            "tmp/fair_bowtie2_mapping_fastp_trimming_single_ended/{sample.sample_id}.fastp.json",
            sample=lookup(
                query="downstream_file != downstream_file & species == '{species}' & release == '{release}' & build == '{build}'",
                within=samples,
            ),
        ),
        fastqc_pair_ended=collect(
            "results/QC/report_pe/{sample.sample_id}.{stream}_fastqc.zip",
            sample=lookup(
                query="downstream_file == downstream_file & species == '{species}' & release == '{release}' & build == '{build}'",
                within=samples,
            ),
            stream=stream_tuple,
        ),
        fastqc_single_ended=collect(
            "results/QC/report_pe/{sample.sample_id}_fastqc.zip",
            sample=lookup(
                query="downstream_file != downstream_file & species == '{species}' & release == '{release}' & build == '{build}'",
                within=samples,
            ),
        ),
        bowtie2=collect(
            "logs/fair_bowtie2_mapping_bowtie2_alignment/{sample.species}.{sample.build}.{sample.release}.{datatype}/{sample.sample_id}.log",
            sample=lookup(
                query="species == '{species}' & release == '{release}' & build == '{build}'",
                within=samples,
            ),
            datatype="{datatype}",
        ),
        samtools=collect(
            "tmp/fair_bowtie2_mapping_samtools_stats/{sample.species}.{sample.build}.{sample.release}.{datatype}/{sample.sample_id}.txt",
            sample=lookup(
                query="species == '{species}' & release == '{release}' & build == '{build}'",
                within=samples,
            ),
            datatype="{datatype}",
        ),
        rseqc_infer_experiment=collect(
            "tmp/fair_bowtie2_mapping_rseqc_infer_experiment/{sample.species}.{sample.build}.{sample.release}.{datatype}/{sample.sample_id}.infer_experiment.txt",
            sample=lookup(
                query="species == '{species}' & release == '{release}' & build == '{build}'",
                within=samples,
            ),
            datatype="{datatype}",
        ),
        rseqc_bamstat=collect(
            "tmp/fair_bowtie2_mapping_rseqc_bamstat/{sample.species}.{sample.build}.{sample.release}.{datatype}/{sample.sample_id}.bamstat.txt",
            sample=lookup(
                query="species == '{species}' & release == '{release}' & build == '{build}'",
                within=samples,
            ),
            datatype="{datatype}",
        ),
        rseqc_read_gc=collect(
            "tmp/fair_bowtie2_mapping_rseqc_read_gc/{sample.species}.{sample.build}.{sample.release}.{datatype}/{sample.sample_id}.GC.xls",
            sample=lookup(
                query="species == '{species}' & release == '{release}' & build == '{build}'",
                within=samples,
            ),
            datatype="{datatype}",
        ),
        rseqc_read_distribution=collect(
            "tmp/fair_bowtie2_mapping_rseqc_read_distribution/{sample.species}.{sample.build}.{sample.release}.{datatype}/{sample.sample_id}.txt",
            sample=lookup(
                query="species == '{species}' & release == '{release}' & build == '{build}'",
                within=samples,
            ),
            datatype="{datatype}",
        ),
        rseqc_inner_distance=collect(
            "tmp/fair_bowtie2_mapping_rseqc_inner_distance/{sample.species}.{sample.build}.{sample.release}.{datatype}/{sample.sample_id}.inner_distance_freq.txt",
            sample=lookup(
                query="species == '{species}' & release == '{release}' & build == '{build}'",
                within=samples,
            ),
            datatype="{datatype}",
        ),
        goleft_indexcov_ped=collect(
            "tmp/fair_bowtie2_mapping_goleft_indexcov/{sample.species}.{sample.release}.{sample.build}.{datatype}/{sample.sample_id}-indexcov.ped",
            sample=lookup(
                query="species == '{species}' & release == '{release}' & build == '{build}'",
                within=samples,
            ),
            datatype="{datatype}",
        ),
        goleft_indexcov_roc=collect(
            "tmp/fair_bowtie2_mapping_goleft_indexcov/{sample.species}.{sample.release}.{sample.build}.{datatype}/{sample.sample_id}-indexcov.roc",
            sample=lookup(
                query="species == '{species}' & release == '{release}' & build == '{build}'",
                within=samples,
            ),
            datatype="{datatype}",
        ),
        ngsderive_readlen=collect(
            "tmp/fair_bowtie2_mapping_ngsderive_readlen/{sample.species}.{sample.build}/{sample.release}.{datatype}/{sample.sample_id}.readlen.tsv",
            sample=lookup(
                query="species == '{species}' & release == '{release}' & build == '{build}'",
                within=samples,
            ),
            datatype="{datatype}",
        ),
        ngsderive_instrument=collect(
            "tmp/fair_bowtie2_mapping_ngsderive_instrument/{sample.species}.{sample.build}/{sample.release}.{datatype}/{sample.sample_id}.instrument.tsv",
            sample=lookup(
                query="species == '{species}' & release == '{release}' & build == '{build}'",
                within=samples,
            ),
            datatype="{datatype}",
        ),
        ngsderive_encoding=collect(
            "tmp/fair_bowtie2_mapping_ngsderive_encoding/{sample.species}.{sample.build}/{sample.release}.{datatype}/{sample.sample_id}.encoding.tsv",
            sample=lookup(
                query="species == '{species}' & release == '{release}' & build == '{build}'",
                within=samples,
            ),
            datatype="{datatype}",
        ),
        fastq_screen_single_ended=branch(
            lookup_config(
                dpath="params/fair_fastqc_multiqc_fastq_screen/fastq_screen_config"
            ),
            then=collect(
                "tmp/fair_fastqc_multiqc/fastq_screen_single_ended/{single_ended_data.sample_id}.fastq_screen.txt",
                single_ended_data=get_single_ended_samples(),
            ),
            otherwise=[],
        ),
        fastq_screen_pair_ended=branch(
            lookup_config(
                dpath="params/fair_fastqc_multiqc_fastq_screen/fastq_screen_config"
            ),
            then=collect(
                "tmp/fair_fastqc_multiqc_fastq_screen_pair_ended/{pair_ended_data.sample_id}.{stream}.fastq_screen.txt",
                pair_ended_data=get_pair_ended_samples(),
                stream=stream_tuple,
            ),
            otherwise=[],
        ),
        macs2=collect(
            "tmp/fair_macs2_calling_macs2_callpeak_{macs2_peak_type}/{sample.species}.{sample.build}.{sample.release}.{datatype}/{macs2_peak_type}/{sample.sample_id}_peaks.xls",
            sample=lookup(
                query="species == '{species}' & release == '{release}' & build == '{build}'",
                within=samples,
            ),
            datatype="{datatype}",
            macs2_peak_type=macs2_peak_types,
        ),
        deeptools_coverage_raw="tmp/fair_macs2_calling_deeptools_plotcoverage/{species}.{build}.{release}.{datatype}/Coverage.raw",
        deeptools_coverage_metrics="tmp/fair_macs2_calling_deeptools_plotcoverage/{species}.{build}.{release}.{datatype}/Coverage.metrics",
        deeptools_fingerprint="tmp/fair_macs2_calling_deeptools_fingerprint/{species}.{build}.{release}.{datatype}/raw_counts.tab",
    output:
        report(
            "results/{species}.{build}.{release}.{datatype}/QC/MultiQC_PeakCalling.html",
            caption="../report/multiqc.rst",
            category="Quality Controls",
            subcategory="General",
            labels={
                "report": "html",
                "step": "Mapping",
                "organism": "{species}.{build}.{release}",
            },
        ),
        "results/{species}.{build}.{release}.{datatype}/QC/MultiQC_Mapping_data.zip",
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: (2 * 1024) * attempt,
        runtime=lambda wildcards, attempt: int(60 * 0.5) * attempt,
        tmpdir=tmp,
    params:
        extra=lookup_config(
            dpath="params/fair_macs2_calling/multiqc",
            default="--verbose --no-megaqc-upload --no-ansi --force",
        ),
        use_input_files_only=True,
    log:
        "logs/fair_macs2_calling/multiqc_report/{species}.{build}.{release}.{datatype}.log",
    benchmark:
        "benchmark/fair_macs2_calling/multiqc_report/{species}.{build}.{release}.{datatype}.tsv"
    wrapper:
        f"{snakemake_wrappers_prefix}/bio/multiqc"
