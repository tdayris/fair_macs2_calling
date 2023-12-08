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


release_list: List[str] = list(set(genomes.release.tolist()))
build_list: List[str] = list(set(genomes.build.tolist()))
species_list: List[str] = list(set(genomes.species.tolist()))


wildcard_constraints:
    sample=r"|".join(samples.sample_id),
    release=r"|".join(release_list),
    build=r"|".join(build_list),
    species=r"|".join(species_list),


def get_sample_information(
    wildcards: snakemake.io.Wildcards, samples: pandas.DataFrame
) -> Dict[str, Optional[str]]:
    """
    Return sample information for a given {sample} wildcards

    Parameters:
    wildcards (snakemake.io.Wildcards): Required for snakemake unpacking function
    samples   (pandas.DataFrame)      : Describe samples and their input files

    Return (Dict[str, Optional[str]]):
    Sample information
    """
    result: Optional[str] = samples.loc[(samples["sample_id"] == str(wildcards.sample))]
    if len(result) > 0:
        return next(iter(result.to_dict(orient="index").values()))
    return defaultdict(lambda: None)


def get_reference_genome_data(
    wildcards: snakemake.io.Wildcards, genomes: pandas.DataFrame
) -> Dict[str, Optional[str]]:
    """
    Return genome information for a given set of {species, build, release} wildcards

    Parameters:
    wildcards (snakemake.io.Wildcards): Required for snakemake unpacking function
    genomes   (pandas.DataFrame)      : Describe genomes and reference file(s)

    Return (Dict[str, Optional[str]]):
    Genome information
    """
    result: Optional[str] = genomes.loc[
        (genomes["species"] == str(wildcards.species))
        & (genomes["build"] == str(wildcards.build))
        & (genomes["release"] == str(wildcards.release))
    ]
    if len(result) > 0:
        return next(iter(result.to_dict(orient="index").values()))
    return defaultdict(lambda: None)


def get_effective_genome_size(
    wildcards: snakemake.io.Wildcards, genomes: pandas.DataFrame = genomes
) -> Optional[str]:
    """
    Return effective genome size if available in genomes table

    Parameters:
    wildcards (snakemake.io.Wildcards): Required for snakemake unpacking function
    genomes   (pandas.DataFrame)      : Genome description and reference file paths

    Return    (Optional[str])
    effective genome size if available else None
    """
    genome_data: Dict[str, Optional[str]] = genomes.loc[
        (genomes.species == str(wildcards.species))
        & (genomes.build == str(wildcards.build))
        & (genomes.release == str(wildcards.release))
    ]
    if len(genome_data) > 0:
        return next(iter(genome_data.to_dict(orient="index").values())).get(
            "effective_genome_size"
        )
    return None


def get_deeptools_bamcoverage_input(
    wildcards: snakemake.io.Wildcards, genome: pandas.DataFrame = genomes
) -> Dict[str, str]:
    """
    Return expected input files for DeepTools bamCoverage, according to user-input,
    and snakemake-wrapper requirements

    Parameters:
    wildcards (snakemake.io.Wildcards): Required for snakemake unpacking function
    genomes   (pandas.DataFrame)      : Genome description and reference file paths

    Return (Dict[str, str]):
    Dictionnary of all input files as required by DeepTools bamCoverage's snakemake-wrapper
    """
    species: str = str(wildcards.species)
    build: str = str(wildcards.build)
    release: str = str(wildcards.release)
    sample_id: str = str(wildcards.sample)
    datatype: str = "dna"

    results: Dict[str, str] = {
        "bam": f"results/{species}.{build}.{release}.{datatype}/Mapping/{sample_id}.bam",
        "bai": f"results/{species}.{build}.{release}.{datatype}/Mapping/{sample_id}.bam.bai",
    }

    genome_data: Dict[str, Optional[str]] = genomes.loc[
        (genomes["species"] == str(wildcards.species))
        & (genomes["build"] == str(wildcards.build))
        & (genomes["release"] == str(wildcards.release))
    ]
    if len(genome_data) > 0:
        blacklist = next(iter(genome_data.to_dict(orient="index").values())).get(
            "blacklist"
        )

    if blacklist:
        results["blacklist"] = blacklist
    elif build in ["GRCh38", "GRCh37", "GRCm38", "NCBIM37"]:
        results[
            "blacklist"
        ] = f"reference/blacklist/{species}.{build}.{release}.merged.bed"

    return results


def get_macs2_callpeak_input(
    wildcards: snakemake.io.Wildcards, samples: pandas.DataFrame = samples
) -> Dict[str, str]:
    """
    Return expected input files for Macs2 peak calling, according to user-input,
    and snakemake-wrapper requirements

    Parameters:
    wildcards (snakemake.io.Wildcards): Required for snakemake unpacking function
    samples   (pandas.DataFrame)      : Describe sample names and related paths/genome

    Return (Dict[str, str]):
    Dictionnary of all input files as required by Macs2's snakemake-wrapper
    """
    sample_data = get_sample_information(wildcards, samples=samples)
    species: str = str(wildcards.species)
    build: str = str(wildcards.build)
    release: str = str(wildcards.release)
    sample_id: str = str(wildcards.sample)
    datatype: str = "dna"

    input_sample = sample_data.get("input")
    if input_sample:
        input_data = get_sample_information(wildcards, samples=samples)
        input_id: str = str(input_data["sample_id"])

        return {
            "treatment": f"results/{species}.{build}.{release}.{datatype}/Mapping/{sample_id}.bam",
            "treatment_index": f"results/{species}.{build}.{release}.{datatype}/Mapping/{sample_id}.bam.bai",
            "control": f"results/{species}.{build}.{release}.{datatype}/Mapping/{input_id}.bam",
            "control_index": f"results/{species}.{build}.{release}.{datatype}/Mapping/{input_id}.bam.bai",
        }

    return {
        "treatment": f"results/{species}.{build}.{release}.{datatype}/Mapping/{sample_id}.bam",
        "treatment_index": f"results/{species}.{build}.{release}.{datatype}/Mapping/{sample_id}.bam.bai",
    }


def get_macs2_callpeak_params(
    wildcards: snakemake.io.Wildcards,
    samples: pandas.DataFrame = samples,
    config: Dict[str, Any] = config,
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
) -> Dict[str, str]:
    """
    Return expected input files for Homer annotate peaks,
    according to user-input, and snakemake-wrapper requirements

    Parameters:
    wildcards (snakemake.io.Wildcards): Required for snakemake unpacking function
    samples   (pandas.DataFrame)      : Describe sample names and related paths/genome
    genomes   (pandas.DataFrame)      : Describe sample names and related paths/genome

    Return (Dict[str, str]):
    Input files dict, as required by homer annotatepeaks's snakemake-wrapper
    """
    datatype: str = "dna"
    species: str = str(wildcards.species)
    build: str = str(wildcards.build)
    release: str = str(wildcards.release)
    sample: str = str(wildcards.sample)
    macs2_peak_type: str = str(wildcards.macs2_peak_type)
    reference: Dict[str, str] = get_reference_genome_data(wildcards, genomes)

    result: Dict[str, str] = {
        "wig": f"results/{species}.{build}.{release}.{datatype}/Coverage/{sample}.bw",
        "peaks": f"results/{species}.{build}.{release}.{datatype}/PeakCalling/{macs2_peak_type}_bed/{sample}.{macs2_peak_type}.bed",
    }

    fasta = reference.get("fasta")
    result["genome"] = (
        fasta or f"reference/{species}.{build}.{release}.{datatype}.fasta"
    )

    fai = reference.get("fasta_index")
    result["genome_index"] = (
        fai or f"reference/{species}.{build}.{release}.{datatype}.fasta.fai"
    )

    gtf = reference.get("gtf")
    result["gtf"] = gtf or f"reference/{species}.{build}.{release}.gtf"

    return result


def get_homer_annotate_peaks_params(
    wildcards: snakemake.io.Wildcards,
    samples: pandas.DataFrame = samples,
    genomes: pandas.DataFrame = genomes,
    config: Dict[str, Any] = config,
) -> Dict[str, Optional[str]]:
    """
    Return expected parameters for Homer annotate peaks,
    according to user-input, and snakemake-wrapper requirements

    Parameters:
    wildcards (snakemake.io.Wildcards): Required for snakemake unpacking function
    samples   (pandas.DataFrame)      : Describe sample names and related paths/genome
    genomes   (pandas.DataFrame)      : Describe sample names and related paths/genome
    config    (Dict[str, Any])        : User defined configuration

    Return (Dict[str, Optional[str]]):
    Parameters, as required by homer annotatepeaks's snakemake-wrapper
    """
    extra: str = config.get("params", {}).get("homer", {}).get("annotatepeaks", "")
    sample_data: Dict[str, Optional[str]] = get_sample_information(wildcards, samples)
    fragment_size: Optional[str] = sample_data.get("fragment_size")
    downstream_file: Optional[str] = sample_data.get("downstream_file")

    if fragment_size and ("fragLength" not in extra) and downstream_file:
        extra += " -fragLength {fragment_size} "

    mode: Optional[str] = "tss"

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
    wildcards: snakemake.io.Wildcards, samples: pandas.DataFrame = samples
) -> Dict[str, List[str]]:
    """
    Return expected input files for deeptools plotCoverage,
    according to user-input, and snakemake-wrapper requirements

    Parameters:
    wildcards (snakemake.io.Wildcards): Required for snakemake unpacking function
    samples   (pandas.DataFrame)      : Describe sample names and related paths/genome

    Return (Dict[str, List[str]]):
    Input files dict, as required by deeptools plotCoverage's snakemake-wrapper
    """
    results: Dict[str, List[str]] = {"bams": [], "bais": []}
    species: str = str(wildcards.species)
    build: str = str(wildcards.build)
    release: str = str(wildcards.release)
    genome_data: Dict[str, Optional[str]] = genomes.loc[
        (genomes["species"] == str(wildcards.species))
        & (genomes["build"] == str(wildcards.build))
        & (genomes["release"] == str(wildcards.release))
    ]
    if len(genome_data) > 0:
        blacklist = next(iter(genome_data.to_dict(orient="index").values())).get(
            "blacklist"
        )

    if blacklist:
        results["blacklist"] = blacklist
    elif build in ["GRCh38", "GRCh37", "GRCm38", "NCBIM37"]:
        results[
            "blacklist"
        ] = f"reference/blacklist/{species}.{build}.{release}.merged.bed"

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
) -> Dict[str, List[str]]:
    """
    Return expected input files for deeptools plotfingerprint,
    according to user-input, and snakemake-wrapper requirements

    Parameters:
    wildcards (snakemake.io.Wildcards): Required for snakemake unpacking function
    samples   (pandas.DataFrame)      : Describe sample names and related paths/genome

    Return (Dict[str, List[str]]):
    Input files dict, as required by deeptools plotfingerprint's snakemake-wrapper
    """
    results: Dict[str, List[str]] = {"bam_files": [], "bam_idx": []}
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


def get_multiqc_report_input(
    wildcards: snakemake.io.Wildcards,
    sample: pandas.DataFrame = samples,
    config: Dict[str, Any] = config,
) -> Dict[str, Union[str, List[str]]]:
    """
    Return expected input files for Multiqc peak calling report,
    according to user-input, and snakemake-wrapper requirements

    Parameters:
    wildcards (snakemake.io.Wildcards): Required for snakemake unpacking function
    samples   (pandas.DataFrame)      : Describe sample names and related paths/genome
    config    (Dict[str, Any])        : User defined configuration

    Return (Dict[str, Union[str, List[str]]]):
    Input files dict, as required by MultiQC's snakemake-wrapper
    """
    results: Dict[str, Union[str, List[str]]] = {
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

        sample_data: Dict[str, Optional[str]] = get_sample_information(
            snakemake.io.Wildcards(fromdict={"sample": sample}), samples
        )
        downstream_file: Optional[str] = sample_data.get("downstream_file")

        if downstream_file:
            results["fastp"].append(f"tmp/fastp/report_pe/{sample}.fastp.json")
            results["fastqc"].append(f"results/QC/report_pe/{sample}.1_fastqc.zip")
            results["fastqc"].append(f"results/QC/report_pe/{sample}.2_fastqc.zip")
        else:
            results["fastp"].append(f"tmp/fastp/report_se/{sample}.fastp.json")
            results["fastqc"].append(f"results/QC/report_pe/{sample}_fastqc.zip")

    peak_types: List[str] = (
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


def get_macs2_calling_pipeline_targets(
    wildcards: snakemake.io.Wildcards,
    samples: pandas.DataFrame = samples,
    config: Dict[str, Any] = config,
) -> Dict[str, List[str]]:
    """
    Return expected output files for this pipeline,
    according to user-input, and snakemake requirements

    Parameters:
    wildcards (snakemake.io.Wildcards): Required for snakemake unpacking function
    samples   (pandas.DataFrame)      : Describe sample names and related paths/genome
    config    (Dict[str, Any])        : User defined configuration file content

    Return (Dict[str, List[str]]):
    Output files dict
    """
    results: Dict[str, List[str]] = {
        "coverage": [],
        "homer": [],
        "bedtools": [],
        "mapping": "tmp/targets/fair_bowtie2_mapping_target.flag",
        "multiqc": [
            "results/QC/MultiQC.html",
            "results/QC/MultiQC.PeakCalling.html",
        ],
    }
    sample_iterator = zip(
        samples.sample_id,
        samples.species,
        samples.build,
        samples.release,
    )
    for sample, species, build, release in sample_iterator:
        results["coverage"].append(
            f"results/{species}.{build}.{release}.dna/Coverage/{sample}.bw"
        )

        results["homer"].append(
            f"results/{species}.{build}.{release}.dna/PeakCalling/narrowPeak_annotation/{sample}.narrowPeak.tsv"
        )

        results["homer"].append(
            f"results/{species}.{build}.{release}.dna/PeakCalling/broadPeak_annotation/{sample}.broadPeak.tsv"
        )

        results["bedtools"].append(
            f"results/{species}.{build}.{release}.dna/Graphs/PlotCoverage.png"
        )

        results["bedtools"].append(
            f"results/{species}.{build}.{release}.dna/Graphs/PlotFingerprint.png"
        )

    return results
