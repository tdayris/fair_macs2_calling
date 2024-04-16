rule fair_macs2_calling_multiqc_config:
    input:
        "tmp/fair_fastqc_multiqc/bigr_logo.png",
    output:
        temp("tmp/fair_macs2_calling/multiqc_config.yaml"),
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 512,
        runtime=lambda wildcards, attempt: attempt * 5,
        tmpdir=tmp,
    localrule: True
    log:
        "logs/fair_macs2_calling/multiqc_config.log",
    benchmark:
        "benchmark/fair_macs2_calling/multiqc_config.tsv"
    params:
        extra=lambda wildcards, input: {
            "title": "Macs2 Peak-Calling quality control report",
            "subtitle": "Produced on raw fastq recieved from sequencer",
            "intro_text": (
                "This pipeline building this report expects "
                "user-provided informations about sequencing "
                "protocol, samples organism, and wet-lab "
                "experimental design."
            ),
            "report_comment": (
                "This report was generated using: "
                "https://github.com/tdayris/fair_macs2_calling"
            ),
            "show_analysis_paths": False,
            "show_analysis_time": False,
            "custom_logo": input[0],
            "custom_logo_url": "https://bioinfo_gustaveroussy.gitlab.io/bigr/webpage/",
            "custom_logo_title": "Bioinformatics Platform @ Gustave Roussy",
            "report_header_info": [
                {"Contact E-mail": "bigr@gustaveroussy.fr"},
                {"Application type": "Short-gapped reads"},
                {"Project Type": "Peak-Calling"},
            ],
            "software_versions": {
                "Quality controls": {
                    "fastqc": "1.12.1",
                    "fastq_screen": "0.15.3",
                    "bowtie2": "1.3.1",
                    "multiqc": "1.20.0",
                },
                "Mapping": {
                    "bowtie2": "2.5.3",
                    "sambamba": "1.0",
                    "samtools": "1.19.2",
                    "picard": "3.1.1",
                    "rseqc": "5.0.3",
                    "fastp": "0.23.4",
                    "ngsderive": "3.3.2",
                    "goleft": "0.2.4",
                },
                "PeakCalling": {
                    "macs2": "2.2.9.1",
                    "homer": "4.11",
                },
                "Pipeline": {
                    "snakemake": "8.5.3",
                    "fair_macs2_calling": "2.1.0",
                    "fair_bowtie2_mapping": "3.3.2",
                    "fair_fastqc_multiqc": "2.2.6",
                    "fair_genome_indexer": "3.4.4",
                },
            },
            "disable_version_detection": True,
            "run_modules": [
                "fastqc",
                "fastq_screen",
                "fastp",
                "bowtie2",
                "samtools",
                "picard",
                "rseqc",
                "ngsderive",
                "goleft_indexcov",
                "macs2",
                "homer",
            ],
            "report_section_order": {
                "fastq_screen": {"order": 1000},
                "ngsderive": {"order": 950},
                "fastqc": {"order": 900},
                "fastp": {"order": 890},
                "bowtie2": {"order": 880},
                "picard": {"order": 870},
                "samtools": {"order": 860},
                "rseqc": {"order": 850},
                "goleft_indexcov": {"order": 840},
                "macs2": {"order": 830},
                "homer": {"order": 820},
                "software_versions": {"order": -1000},
            },
        },
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/fair_macs2_calling_multiqc_config.py"


rule fair_macs2_calling_multiqc_report:
    input:
        config="tmp/fair_macs2_calling/multiqc_config.yaml",
        picard_qc=collect(
            "tmp/fair_bowtie2_mapping/picard_create_multiple_metrics/{sample.species}.{sample.build}.{sample.release}.dna/stats/{sample.sample_id}{ext}",
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
        ),
        fastp_pair_ended=collect(
            "tmp/fair_bowtie2_mapping/fastp_trimming_pair_ended/{sample.sample_id}.fastp.json",
            sample=lookup(
                query="downstream_file == downstream_file & species == '{species}' & release == '{release}' & build == '{build}'",
                within=samples,
            ),
        ),
        fastp_single_ended=collect(
            "tmp/fair_bowtie2_mapping/fastp_trimming_single_ended/{sample.sample_id}.fastp.json",
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
            stream=stream_list,
        ),
        fastqc_single_ended=collect(
            "results/QC/report_pe/{sample.sample_id}_fastqc.zip",
            sample=lookup(
                query="downstream_file != downstream_file & species == '{species}' & release == '{release}' & build == '{build}'",
                within=samples,
            ),
        ),
        bowtie2=collect(
            "logs/fair_bowtie2_mapping/bowtie2_alignment/{sample.species}.{sample.build}.{sample.release}.dna/{sample.sample_id}.log",
            sample=lookup(
                query="species == '{species}' & release == '{release}' & build == '{build}'",
                within=samples,
            ),
        ),
        samtools=collect(
            "tmp/fair_bowtie2_mapping/samtools_stats/{sample.species}.{sample.build}.{sample.release}.dna/{sample.sample_id}.txt",
            sample=lookup(
                query="species == '{species}' & release == '{release}' & build == '{build}'",
                within=samples,
            ),
        ),
        rseqc_infer_experiment=collect(
            "tmp/fair_bowtie2_mapping/rseqc_infer_experiment/{sample.species}.{sample.build}.{sample.release}.dna/{sample.sample_id}.infer_experiment.txt",
            sample=lookup(
                query="species == '{species}' & release == '{release}' & build == '{build}'",
                within=samples,
            ),
        ),
        rseqc_bamstat=collect(
            "tmp/fair_bowtie2_mapping/rseqc_bamstat/{sample.species}.{sample.build}.{sample.release}.dna/{sample.sample_id}.bamstat.txt",
            sample=lookup(
                query="species == '{species}' & release == '{release}' & build == '{build}'",
                within=samples,
            ),
        ),
        rseqc_read_gc=collect(
            "tmp/fair_bowtie2_mapping/rseqc_read_gc/{sample.species}.{sample.build}.{sample.release}.dna/{sample.sample_id}.GC.xls",
            sample=lookup(
                query="species == '{species}' & release == '{release}' & build == '{build}'",
                within=samples,
            ),
        ),
        rseqc_read_distribution=collect(
            "tmp/fair_bowtie2_mapping/rseqc_read_distribution/{sample.species}.{sample.build}.{sample.release}.dna/{sample.sample_id}.txt",
            sample=lookup(
                query="species == '{species}' & release == '{release}' & build == '{build}'",
                within=samples,
            ),
        ),
        rseqc_inner_distance=collect(
            "tmp/fair_bowtie2_mapping/rseqc_inner_distance/{sample.species}.{sample.build}.{sample.release}.dna/{sample.sample_id}.inner_distance_freq.txt",
            sample=lookup(
                query="species == '{species}' & release == '{release}' & build == '{build}'",
                within=samples,
            ),
        ),
        goleft_indexcov_ped=collect(
            "tmp/fair_bowtie2_mapping/goleft/indexcov/{sample.species}.{sample.release}.{sample.build}.dna/{sample.sample_id}-indexcov.ped",
            sample=lookup(
                query="species == '{species}' & release == '{release}' & build == '{build}'",
                within=samples,
            ),
        ),
        goleft_indexcov_roc=collect(
            "tmp/fair_bowtie2_mapping/goleft/indexcov/{sample.species}.{sample.release}.{sample.build}.dna/{sample.sample_id}-indexcov.roc",
            sample=lookup(
                query="species == '{species}' & release == '{release}' & build == '{build}'",
                within=samples,
            ),
        ),
        ngsderive_readlen=collect(
            "tmp/fair_bowtie2_mapping/ngsderive/readlen/{sample.species}.{sample.build}/{sample.release}.dna/{sample.sample_id}.readlen.tsv",
            sample=lookup(
                query="species == '{species}' & release == '{release}' & build == '{build}'",
                within=samples,
            ),
        ),
        ngsderive_instrument=collect(
            "tmp/fair_bowtie2_mapping/ngsderive/instrument/{sample.species}.{sample.build}/{sample.release}.dna/{sample.sample_id}.instrument.tsv",
            sample=lookup(
                query="species == '{species}' & release == '{release}' & build == '{build}'",
                within=samples,
            ),
        ),
        ngsderive_encoding=collect(
            "tmp/fair_bowtie2_mapping/ngsderive/encoding/{sample.species}.{sample.build}/{sample.release}.dna/{sample.sample_id}.encoding.tsv",
            sample=lookup(
                query="species == '{species}' & release == '{release}' & build == '{build}'",
                within=samples,
            ),
        ),
        macs2=collect(
            "tmp/fair_macs2_calling/macs2/{sample.species}.{sample.build}.{sample.release}.dna/broadPeak/{sample.sample_id}_peaks.xls",
            sample=lookup(
                query="species == '{species}' & release == '{release}' & build == '{build}'",
                within=samples,
            ),
        ),
        deeptools_coverage_raw="tmp/fair_macs2_calling/deeptools/plot_coverage/{species}.{build}.{release}.dna/Coverage.raw",
        deeptools_coverage_metrics="tmp/fair_macs2_calling/deeptools/plot_coverage/{species}.{build}.{release}.dna/Coverage.metrics",
        deeptools_fingerprint="tmp/fair_macs2_calling/deeptools/plot_fingerprint/{species}.{build}.{release}.dna/raw_counts.tab",
    output:
        report(
            "results/{species}.{build}.{release}.dna/QC/MultiQC_PeakCalling.html",
            caption="../report/multiqc.rst",
            category="Quality Controls",
            subcategory="General",
            labels={
                "report": "html",
                "step": "Mapping",
                "organism": "{species}.{build}.{release}",
            },
        ),
        "results/{species}.{build}.{release}.dna/QC/MultiQC_Mapping_data.zip",
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
        "logs/fair_macs2_calling/multiqc_report/{species}.{build}.{release}.dna.log",
    benchmark:
        "benchmark/fair_macs2_calling/multiqc_report/{species}.{build}.{release}.dna.tsv"
    wrapper:
        f"{snakemake_wrappers_prefix}/bio/multiqc"
