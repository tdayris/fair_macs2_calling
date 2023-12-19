import csv
import pandas
import snakemake
import snakemake.utils

from collections import defaultdict
from typing import Any, Dict, List, Optional, Union

snakemake.utils.min_version("7.29.0")

# containerized: "docker://snakemake/snakemake:v7.32.4"
# containerized: "docker://mambaorg/micromamba:git-8440cec-jammy-cuda-12.2.0"
# containerized: "docker://condaforge/mambaforge:23.3.1-1"


# Load and check configuration file
configfile: "config/config.yaml"


snakemake.utils.validate(config, "../schemas/config.schema.yaml")

# Load and check samples properties table
sample_table_path: str = config.get("samples", "config/samples.csv")
with open(sample_table_path, "r") as sample_table_stream:
    dialect: csv.Dialect = csv.Sniffer().sniff(sample_table_stream.read(1024))
    sample_table_stream.seek(0)

samples: pandas.DataFrame = pandas.read_csv(
    filepath_or_buffer=sample_table_path,
    sep=dialect.delimiter,
    header=0,
    index_col=None,
    comment="#",
    dtype=str,
)
samples = samples.where(samples.notnull(), None)
snakemake.utils.validate(samples, "../schemas/samples.schema.yaml")

# This is here for compatibility with
genome_table_path: str = config.get("genomes")
if genome_table_path:
    with open(genome_table_path, "r") as genome_table_stream:
        dialect: csv.Dialect = csv.Sniffer().sniff(genome_table_stream.read(1024))
        genome_table_stream.seek(0)

    genomes: pandas.DataFrame = pandas.read_csv(
        filepath_or_buffer=genome_table_path,
        sep=dialect.delimiter,
        header=0,
        index_col=None,
        comment="#",
        dtype=str,
    )
    genomes = genomes.where(genomes.notnull(), None)
else:
    genomes: pandas.DataFrame = samples[
        ["species", "build", "release"]
    ].drop_duplicates(keep="first", ignore_index=True)
    genomes.to_csv("genomes.csv", sep=",", index=False, header=True)
    config["genomes"] = "genomes.csv"

snakemake.utils.validate(genomes, "../schemas/genomes.schema.yaml")

snakemake_wrappers_version: str = "v3.0.0"


report: "../report/workflow.rst"


release_list: list[str] = list(set(genomes.release.tolist()))
build_list: list[str] = list(set(genomes.build.tolist()))
species_list: list[str] = list(set(genomes.species.tolist()))
datatypes: list[str] = ["dna", "cdna"]
macs2_peak_types: list[str] = ["broadPeak", "narrowPeak"]


wildcard_constraints:
    sample=r"|".join(samples.sample_id),
    release=r"|".join(release_list),
    build=r"|".join(build_list),
    species=r"|".join(species_list),
    datatype=r"|".join(datatypes),
    macs2_peak_type=r"|".join(macs2_peak_types),


def get_sample_information(
    wildcards: snakemake.io.Wildcards, samples: pandas.DataFrame
) -> dict[str, str | None]:
    """
    Return sample information for a given {sample} wildcards

    Parameters:
    wildcards (snakemake.io.Wildcards): Required for snakemake unpacking function
    samples   (pandas.DataFrame)      : Describe samples and their input files

    Return (dict[str, str | None]):
    Sample information
    """
    result: str | None = samples.loc[(samples["sample_id"] == str(wildcards.sample))]
    if len(result) > 0:
        return next(iter(result.to_dict(orient="index").values()))
    return defaultdict(lambda: None)


def get_reference_genome_data(
    wildcards: snakemake.io.Wildcards, genomes: pandas.DataFrame
) -> dict[str, str | None]:
    """
    Return genome information for a given set of {species, build, release} wildcards

    Parameters:
    wildcards (snakemake.io.Wildcards): Required for snakemake unpacking function
    genomes   (pandas.DataFrame)      : Describe genomes and reference file(s)

    Return (dict[str, str | None]):
    Genome information
    """
    result: str | None = genomes.loc[
        (genomes["species"] == str(wildcards.species))
        & (genomes["build"] == str(wildcards.build))
        & (genomes["release"] == str(wildcards.release))
    ]
    if len(result) > 0:
        return next(iter(result.to_dict(orient="index").values()))
    return defaultdict(lambda: None)


def get_effective_genome_size(
    wildcards: snakemake.io.Wildcards, genomes: pandas.DataFrame = genomes
) -> str | None:
    """
    Return effective genome size if available in genomes table

    Parameters:
    wildcards (snakemake.io.Wildcards): Required for snakemake unpacking function
    genomes   (pandas.DataFrame)      : Genome description and reference file paths

    Return    (str | None)
    effective genome size if available else None
    """
    genome_data: dict[str, str | None] = genomes.loc[
        (genomes.species == str(wildcards.species))
        & (genomes.build == str(wildcards.build))
        & (genomes.release == str(wildcards.release))
    ]
    if len(genome_data) > 0:
        return next(iter(genome_data.to_dict(orient="index").values())).get(
            "effective_genome_size"
        )
    return None


def get_blacklist(wildcards: snakemake.io.Wildcards, genomes: pandas.DataFrame = genomes) -> str | None:
    """
    Return blacklist file, if any.

    Parameters:
    wildcards (snakemake.io.Wildcards): Required for snakemake unpacking function
    genomes   (pandas.DataFrame)      : Genome description and reference file paths

    Return: str | None
    If blacklist exists, then return it. Else, return None
    """
    blacklist: str | None = None
    species: str = str(wildcards.species)
    build: str = str(wildcards.build)
    release: str = str(wildcards.release)
    genome_data: dict[str, str | None] = genomes.loc[
        (genomes["species"] == species)
        & (genomes["build"] == build)
        & (genomes["release"] == release)
    ]
    if len(genome_data) > 0:
        return next(iter(genome_data.to_dict(orient="index").values())).get(
            "blacklist", f"reference/blacklist/{species}.{build}.{release}.merged.bed"
        )


def get_deeptools_bamcoverage_input(
    wildcards: snakemake.io.Wildcards, genomes: pandas.DataFrame = genomes
) -> dict[str, str]:
    """
    Return expected input files for DeepTools bamCoverage, according to user-input,
    and snakemake-wrapper requirements

    Parameters:
    wildcards (snakemake.io.Wildcards): Required for snakemake unpacking function
    genomes   (pandas.DataFrame)      : Genome description and reference file paths

    Return (dict[str, str]):
    Dictionnary of all input files as required by DeepTools bamCoverage's snakemake-wrapper
    """
    species: str = str(wildcards.species)
    build: str = str(wildcards.build)
    release: str = str(wildcards.release)
    sample_id: str = str(wildcards.sample)
    datatype: str = "dna"

    results: dict[str, str] = {
        "bam": f"results/{species}.{build}.{release}.{datatype}/Mapping/{sample_id}.bam",
        "bai": f"results/{species}.{build}.{release}.{datatype}/Mapping/{sample_id}.bam.bai",
    }

    blacklist: str | None = get_blacklist(wildcards, genomes=genomes)
    if blacklist:
        results["blacklist"] = blacklist

    return results


def get_macs2_callpeak_input(
    wildcards: snakemake.io.Wildcards, samples: pandas.DataFrame = samples
) -> dict[str, str]:
    """
    Return expected input files for Macs2 peak calling, according to user-input,
    and snakemake-wrapper requirements

    Parameters:
    wildcards (snakemake.io.Wildcards): Required for snakemake unpacking function
    samples   (pandas.DataFrame)      : Describe sample names and related paths/genome

    Return (dict[str, str]):
    Dictionnary of all input files as required by Macs2's snakemake-wrapper
    """
    sample_data = get_sample_information(wildcards, samples=samples)
    species: str = str(wildcards.species)
    build: str = str(wildcards.build)
    release: str = str(wildcards.release)
    sample_id: str = str(wildcards.sample)
    datatype: str = "dna"

    input_sample: str | None = sample_data.get("input")
    treatment: str = f"results/{species}.{build}.{release}.{datatype}/Mapping/{sample_id}.bam"
    treatment_index: str = f"{treatment}.bai"
    if input_sample:
        input_data = get_sample_information(wildcards, samples=samples)
        input_id: str = str(input_data["sample_id"])
        control: str = f"results/{species}.{build}.{release}.{datatype}/Mapping/{input_id}.bam"
        control_index: str = f"{control}.bai"

        return {
            "treatment": treatment,
            "treatment_index": treatment_index,
            "control": control,
            "control_index": control_index,
        }

    return {
        "treatment": treatment,
        "treatment_index": treatment_index,
    }


def get_macs2_callpeak_params(
    wildcards: snakemake.io.Wildcards,
    samples: pandas.DataFrame = samples,
    config: dict[str, Any] = config,
    genomes: pandas.DataFrame = genomes,
) -> str:
    """
    Return expected parameters for Macs2 peak calling, according to user-input,
    and snakemake-wrapper requirements

    Parameters:
    wildcards (snakemake.io.Wildcards): Required for snakemake unpacking function
    samples   (pandas.DataFrame)      : Describe sample names and related paths/genome
    genomes   (pandas.DataFrame)      : Describe genomes

    Return (str):
    Parameters as string, as required by Macs2's snakemake-wrapper
    """
    results: str = config.get("params", {}).get("macs2", {}).get("callpeak", "")

    sample_data = get_sample_information(wildcards, samples=samples)
    if sample_data.get("downstream_file", False):
        if " BAMPE " not in results:
            results += " --format BAMPE "

        if "--nomodel" not in results:
            results += " --nomodel "
    elif "--nolambda" not in results:
        results += " --nolambda "

    species: str = str(sample_data["species"])
    if species.lower() == "homo_sapiens" and " hs " not in results:
        results += " --gsize hs "
    elif species.lower() == "mus_musculus" and " mm " not in results:
        results += " --gsize mm "
    elif species.lower() == "drosophila_melanogaster" and " dm " not in results:
        results += " --gsize dm "
    elif species.lower() == "caenorhabditis_elegans" and " ce " not in results:
        results += " --gsize ce "

    return results


def get_homer_annotate_peaks_input(
    wildcards: snakemake.io.Wildcards,
    genomes: pandas.DataFrame = genomes,
    samples: pandas.DataFrame = samples,
) -> dict[str, str]:
    """
    Return expected input files for Homer annotate peaks,
    according to user-input, and snakemake-wrapper requirements

    Parameters:
    wildcards (snakemake.io.Wildcards): Required for snakemake unpacking function
    samples   (pandas.DataFrame)      : Describe sample names and related paths/genome
    genomes   (pandas.DataFrame)      : Describe sample names and related paths/genome

    Return (dict[str, str]):
    Input files dict, as required by homer annotatepeaks's snakemake-wrapper
    """
    datatype: str = "dna"
    species: str = str(wildcards.species)
    build: str = str(wildcards.build)
    release: str = str(wildcards.release)
    sample: str = str(wildcards.sample)
    macs2_peak_type: str = str(wildcards.macs2_peak_type)

    reference: dict[str, str] = get_reference_genome_data(wildcards, genomes)

    wig: str =  f"results/{species}.{build}.{release}.{datatype}/Coverage/{sample}.bw"
    peaks: str = f"results/{species}.{build}.{release}.{datatype}/PeakCalling/{macs2_peak_type}/{sample}.{macs2_peak_type}.bed"
    fasta: str = reference.get("fasta", f"reference/{species}.{build}.{release}.{datatype}.fasta")
    fai: str = reference.get("fasta_index", f"reference/{species}.{build}.{release}.{datatype}.fasta.fai")
    gtf: str = reference.get("gtf", f"reference/{species}.{build}.{release}.gtf")

    return {"wig": wig, "peaks": peaks, "fasta": fasta, "fai": fai, "gtf": gtf}


def get_homer_annotate_peaks_params(
    wildcards: snakemake.io.Wildcards,
    samples: pandas.DataFrame = samples,
    genomes: pandas.DataFrame = genomes,
    config: dict[str, Any] = config,
) -> dict[str, str | None]:
    """
    Return expected parameters for Homer annotate peaks,
    according to user-input, and snakemake-wrapper requirements

    Parameters:
    wildcards (snakemake.io.Wildcards): Required for snakemake unpacking function
    samples   (pandas.DataFrame)      : Describe sample names and related paths/genome
    genomes   (pandas.DataFrame)      : Describe sample names and related paths/genome
    config    (dict[str, Any])        : User defined configuration

    Return (dict[str, str | None]):
    Parameters, as required by homer annotatepeaks's snakemake-wrapper
    """
    extra: str = config.get("params", {}).get("homer", {}).get("annotatepeaks", "")
    sample_data: dict[str, str | None] = get_sample_information(wildcards, samples)
    fragment_size: str | None = sample_data.get("fragment_size")
    downstream_file: str | None = sample_data.get("downstream_file")

    if fragment_size and ("fragLength" not in extra) and downstream_file:
        extra += " -fragLength {fragment_size} "

    mode: str | None = "tss"

    species: str = str(wildcards.species)
    release: str = str(wildcards.release)
    if species == "homo_sapiens":
        if release == "GRCh38":
            mode = "tss hg38"
        elif release == "GRCh37":
            mode = "tss hg19"
    elif species == "mus_musculus":
        if release == "GRCm38":
            mode = "tss mm10"

    return {
        "extra": extra,
        "mode": mode,
    }


def get_deeptools_plotcoverage_input(
    wildcards: snakemake.io.Wildcards, samples: pandas.DataFrame = samples, genomes: pandas.DataFrame = genomes
) -> dict[str, list[str]]:
    """
    Return expected input files for deeptools plotCoverage,
    according to user-input, and snakemake-wrapper requirements

    Parameters:
    wildcards (snakemake.io.Wildcards): Required for snakemake unpacking function
    samples   (pandas.DataFrame)      : Describe sample names and related paths/genome
    genomes   (pandas.DataFrame)      : Describe genome annotations

    Return (dict[str, list[str]]):
    Input files dict, as required by deeptools plotCoverage's snakemake-wrapper
    """
    results: dict[str, list[str]] = {"bams": [], "bais": []}
    species: str = str(wildcards.species)
    build: str = str(wildcards.build)
    release: str = str(wildcards.release)
    blacklist: str | None = get_blacklist(wildcards, genomes=genomes)
    if blacklist:
        results["blacklist"] = blacklist

    genome_restricted_samples: pandas.DataFrame = samples.loc[
        (samples.species == species)
        & (samples.build == build)
        & (samples.release == release)
    ]

    datatype: str = "dna"
    sample_iterator = zip(
        genome_restricted_samples.sample_id,
        genome_restricted_samples.species,
        genome_restricted_samples.build,
        genome_restricted_samples.release,
    )

    for sample, species, build, release in sample_iterator:
        results["bams"].append(
            f"results/{species}.{build}.{release}.{datatype}/Mapping/{sample}.bam"
        )
        results["bais"].append(
            f"results/{species}.{build}.{release}.{datatype}/Mapping/{sample}.bam.bai"
        )

    return results


def get_deeptools_fingerprint_input(
    wildcards: snakemake.io.Wildcards, samples: pandas.DataFrame = samples
) -> dict[str, list[str]]:
    """
    Return expected input files for deeptools plotfingerprint,
    according to user-input, and snakemake-wrapper requirements

    Parameters:
    wildcards (snakemake.io.Wildcards): Required for snakemake unpacking function
    samples   (pandas.DataFrame)      : Describe sample names and related paths/genome

    Return (dict[str, list[str]]):
    Input files dict, as required by deeptools plotfingerprint's snakemake-wrapper
    """
    results: dict[str, list[str]] = {"bam_files": [], "bam_idx": []}
    species: str = str(wildcards.species)
    build: str = str(wildcards.build)
    release: str = str(wildcards.release)

    genome_restricted_samples: pandas.DataFrame = samples.loc[
        (samples.species == species)
        & (samples.build == build)
        & (samples.release == release)
    ]

    datatype: str = "dna"
    sample_iterator = zip(
        genome_restricted_samples.sample_id,
        genome_restricted_samples.species,
        genome_restricted_samples.build,
        genome_restricted_samples.release,
    )

    for sample, species, build, release in sample_iterator:
        results["bam_files"].append(
            f"results/{species}.{build}.{release}.{datatype}/Mapping/{sample}.bam"
        )
        results["bam_idx"].append(
            f"results/{species}.{build}.{release}.{datatype}/Mapping/{sample}.bam.bai"
        )

    return results


def get_deeptools_multibigwig_summary_input(
    wildcards: snakemake.io.Wildcards,
    samples: pandas.DataFrame = samples,
    genomes: pandas.DataFrame = genomes,
    config: dict[str, Any] = config,
) -> dict[str, list[str] | str]:
    """
    Return expected input files for DeepTools MultiBigwigSummary,
    according to user-input, and snakemake-wrapper requirements

    Parameters:
    wildcards (snakemake.io.Wildcards): Required for snakemake unpacking function
    samples   (pandas.DataFrame)      : Describe sample names and related paths/genome
    genomes   (pandas.DataFrame)      : Describe genome annotations
    config    (dict[str, Any])        : User defined configuration

    Return (dict[str, Union[str, list[str]]]):
    Input files dict, as required by DeepTools MultiBigwigSummary's snakemake-script
    """
    species: str = str(wildcards.species)
    build: str = str(wildcards.build)
    release: str = str(wildcards.release)
    datatype: str = "dna"

    samples_list: list[str] = list(samples.loc[
        (samples.species == species)
        & (samples.build == build)
        & (samples.release == release)
    ].sample_id)

    results: dict[str, list[str] | str] = {
        "bed": f"tmp/bedtools/merge/{species}.{build}.{release}.{datatype}/{macs2_peak_type}.merged.bed",
        "bw": expand(
            "results/{species}.{build}.{release}.{datatype}/Coverage/{sample}.bw",
            sample=samples_list,
            species=[species],
            build=[build],
            release=[release],
            datatype=[datatype],
        ),
    }

    blacklist: str | None = get_blacklist(wildcards, genome)
    if blacklist:
        results["blacklist"] = blacklist

    return results


def get_multiqc_report_input(
    wildcards: snakemake.io.Wildcards,
    sample: pandas.DataFrame = samples,
    config: dict[str, Any] = config,
) -> dict[str, Union[str, list[str]]]:
    """
    Return expected input files for Multiqc peak calling report,
    according to user-input, and snakemake-wrapper requirements

    Parameters:
    wildcards (snakemake.io.Wildcards): Required for snakemake unpacking function
    samples   (pandas.DataFrame)      : Describe sample names and related paths/genome
    config    (dict[str, Any])        : User defined configuration

    Return (dict[str, Union[str, list[str]]]):
    Input files dict, as required by MultiQC's snakemake-wrapper
    """
    results: dict[str, Union[str, list[str]]] = {
        "picard_qc": [],
        "fastp": [],
        "fastqc": [],
        "macs2": [],
        "deeptools_coverage": [],
        "deeptools_fingerprint": [],
        "deeptools_pca": [],
        "deeptools_correlation": [],
        "bowtie2": [],
    }

    datatype: str = "dna"
    sample_iterator = zip(
        samples.sample_id,
        samples.species,
        samples.build,
        samples.release,
    )

    for sample, species, build, release in sample_iterator:
        results["picard_qc"] += multiext(
            f"tmp/picard/{species}.{build}.{release}.{datatype}/stats/{sample}",
            ".alignment_summary_metrics",
            ".insert_size_metrics",
            ".insert_size_histogram.pdf",
            ".base_distribution_by_cycle_metrics",
            ".base_distribution_by_cycle.pdf",
            ".gc_bias.detail_metrics",
            ".gc_bias.summary_metrics",
            ".gc_bias.pdf",
        )

        results["macs2"].append(
            f"tmp/macs2/{species}.{build}.{release}.dna/narrowPeak/{sample}_peaks.xls"
        )

        results["macs2"].append(
            f"tmp/macs2/{species}.{build}.{release}.dna/broadPeak/{sample}_peaks.xls"
        )

        results["deeptools_coverage"].append(
            f"tmp/deeptools/plot_coverage/{species}.{build}.{release}.{datatype}/Coverage.raw"
        )
        results["deeptools_coverage"].append(
            f"tmp/deeptools/plot_coverage/{species}.{build}.{release}.{datatype}/Coverage.metrics"
        )

        results["deeptools_fingerprint"].append(
            f"tmp/deeptools/plot_fingerprint/{species}.{build}.{release}.{datatype}/raw_counts.tab"
        )

        results["deeptools_fingerprint"].append(
            f"tmp/deeptools/plot_fingerprint/{species}.{build}.{release}.{datatype}/qc_metrics.txt"
        )

        results[
            "bowtie2"
        ] = f"logs/bowtie2/align/{species}.{build}.{release}.{datatype}/{sample}.log"

        sample_data: dict[str, str | None] = get_sample_information(
            snakemake.io.Wildcards(fromdict={"sample": sample}), samples
        )
        downstream_file: str | None = sample_data.get("downstream_file")

        if downstream_file:
            results["fastp"].append(f"tmp/fastp/report_pe/{sample}.fastp.json")
            results["fastqc"].append(f"results/QC/report_pe/{sample}.1_fastqc.zip")
            results["fastqc"].append(f"results/QC/report_pe/{sample}.2_fastqc.zip")
        else:
            results["fastp"].append(f"tmp/fastp/report_se/{sample}.fastp.json")
            results["fastqc"].append(f"results/QC/report_pe/{sample}_fastqc.zip")

    peak_types: list[str] = (
        config.get("params", {})
        .get("macs2", {})
        .get("modes", ["broadPeak", "narrowPeak"])
    )
    for macs2_peak_type in peak_types:
        results["deeptools_pca"].append(
            f"tmp/deeptools/plot_pca/{species}.{build}.{release}.{datatype}/{macs2_peak_type}.tab"
        )
        results["deeptools_correlation"].append(
            f"tmp/deeptools/plot_correlation/{species}.{build}.{release}.{datatype}/{macs2_peak_type}.tab"
        )

    return results


def get_merge_homer_summaries_input(
    wildcards: snakemake.io.Wildcards, samples: pandas.DataFrame = samples
) -> dict[str, list[str]]:
    """
    Return expected input files for Multiqc peak calling report,
    according to user-input, and snakemake-wrapper requirements

    Parameters:
    wildcards (snakemake.io.Wildcards): Required for snakemake unpacking function
    samples   (pandas.DataFrame)      : Describe sample names and related paths/genome

    Return (dict[str, list[str]]):
    Input files dict, as required by MultiQC's snakemake-wrapper
    """
    results: dict[str, Union[str, list[str]]] = {"summaries": []}
    species_req: str = str(wildcards.species)
    build_req: str = str(wildcards.build)
    release_req: str = str(wildcards.release)
    macs2_peak_type: str = str(wildcards.macs2_peak_type)
    datatype: str = "dna"
    filtered_samples = samples.loc[
        (samples.species == species_req)
        & (samples.build == build_req)
        & (samples.release == release_req)
    ]

    sample_iterator = zip(
        filtered_samples.sample_id,
        filtered_samples.species,
        filtered_samples.build,
        filtered_samples.release,
    )

    for sample, species, build, release in sample_iterator:
        results["summaries"].append(
            f"tmp/summarize_homer/{species}.{build}.{release}.{datatype}/{sample}.{macs2_peak_type}.json"
        )

    return results


def get_macs2_calling_pipeline_targets(
    wildcards: snakemake.io.Wildcards,
    samples: pandas.DataFrame = samples,
    config: dict[str, Any] = config,
) -> dict[str, list[str]]:
    """
    Return expected output files for this pipeline,
    according to user-input, and snakemake requirements

    Parameters:
    wildcards (snakemake.io.Wildcards): Required for snakemake unpacking function
    samples   (pandas.DataFrame)      : Describe sample names and related paths/genome
    config    (dict[str, Any])        : User defined configuration file content

    Return (dict[str, list[str]]):
    Output files dict
    """
    results: dict[str, list[str]] = {
        "coverage": [],
        "homer": [],
        "bedtools": [],
        "mapping": "tmp/targets/fair_bowtie2_mapping_target.flag",
        "multiqc": [
            "results/QC/MultiQC.html",
            "results/QC/MultiQC.PeakCalling.html",
        ],
        "datavzrd": [],
        "inhouse": [],
    }
    sample_iterator = zip(
        samples.sample_id,
        samples.species,
        samples.build,
        samples.release,
    )

    macs2_peak_types: list[str] = ["narrowPeak", "broadPeak"]

    for sample, species, build, release in sample_iterator:
        results["coverage"].append(
            f"results/{species}.{build}.{release}.dna/Coverage/{sample}.bw"
        )

        results["bedtools"].append(
            f"results/{species}.{build}.{release}.dna/Graphs/PlotCoverage.png"
        )

        results["bedtools"].append(
            f"results/{species}.{build}.{release}.dna/Graphs/PlotFingerprint.png"
        )

        for macs2_peak_type in macs2_peak_types:
            results["datavzrd"].append(
                f"results/{species}.{build}.{release}.dna/PeakCalling/{macs2_peak_type}/{sample}_reports"
            )

            results["homer"].append(
                f"results/{species}.{build}.{release}.dna/PeakCalling/{macs2_peak_type}/{sample}.{macs2_peak_type}.tsv"
            )

            for content in ["Annotations", "GeneTypes"]:
                results["inhouse"].append(
                    f"results/{species}.{build}.{release}.dna/Graphs/{macs2_peak_type}/{content}.catplot.png"
                )

                # results["inhouse"].append(
                #     f"results/{species}.{build}.{release}.dna/Graphs/{macs2_peak_type}/{content}.dotplot.png"
                # )

    return results
