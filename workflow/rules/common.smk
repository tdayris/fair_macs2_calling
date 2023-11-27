import csv
import pandas
import snakemake
import snakemake.utils

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
else:
    genomes: pandas.DataFrame = samples[
        ["species", "build", "release"]
    ].drop_duplicates(keep="first", ignore_index=True)
    genomes.to_csv("genomes.csv", sep=",", index=False, header=True)
    config["genomes"] = "genomes.csv"

snakemake.utils.validate(genomes, "../schemas/genomes.schema.yaml")

snakemake_wrappers_version: str = "v2.13.0"

report: "../report/workflows.rst"

release_list: List[str] = list(set(genomes.release.tolist()))
build_list: List[str] = list(set(genomes.build.tolist()))
species_list: List[str] = list(set(genomes.species.tolist()))

wildcard_constraints:
    sample=r"|".join(samples.sample_id),
    release=r"|".join(release_list),
    build=r"|".join(build_list),
    species=r"|".join(species_list),


def get_sample_data_as_dict(sample: str, samples: pandas.DataFrame = samples) -> Dict[str, Any]:
    """
    Return single-sample information as a dictionary
    """
    return samples[
        samples.sample_id == sample
    ].to_dict(orient="index")[0]


def get_macs2_callpeak_input(wildcards: snakemake.io.Wildcards, samples: pandas.DataFrame = samples) -> Dict[str, str]:
    """
    Return expected input files for Macs2 peak calling, according to user-input,
    and snakemake-wrapper requirements

    Parameters:
    wildcards (snakemake.io.Wildcards): Required for snakemake unpacking function
    samples   (pandas.DataFrame)      : Describe sample names and related paths/genome

    Return (Dict[str, str]):
    Dictionnary of all input files as required by Macs2's snakemake-wrapper
    """
    sample_data = get_sample_data_as_dict(sample=str(wildcards.sample), samples=samples)
    species: str = str(wildcards.species)
    build: str = str(wildcards.build)
    release: str = str(wildcards.release)
    sample_id: str = str(wildcards.sample)
    datatype: str = "dna"

    input_sample = sample_data.get("input")
    if input_sample:
        input_data = get_sample_data_as_dict(sample=input_sample, samples=samples)
        input_id: str = str(input_data["sample_id"])

        return {
            treatment: f"results/Mapping/{species}.{build}.{release}.{datatype}/{sample_id}.bam",
            treatment_index: f"results/Mapping/{species}.{build}.{release}.{datatype}/{sample_id}.bam.bai",
            control: f"results/Mapping/{species}.{build}.{release}.{datatype}/{input_id}.bam",
            control_index: f"results/Mapping/{species}.{build}.{release}.{datatype}/{input_id}.bam.bai",
        }
    
    return {
        treatment: f"results/Mapping/{species}.{build}.{release}.{datatype}/{sample_id}.bam",
        treatment_index: f"results/Mapping/{species}.{build}.{release}.{datatype}/{sample_id}.bam.bai",
    }



def get_macs2_callpeak_params(wildcards: snakemake.io.Wildcards, samples: pandas.DataFrame = samples, config: Dict[str, Any] = config) -> str:
    """
    Return expected parameters for Macs2 peak calling, according to user-input,
    and snakemake-wrapper requirements

    Parameters:
    wildcards (snakemake.io.Wildcards): Required for snakemake unpacking function
    samples   (pandas.DataFrame)      : Describe sample names and related paths/genome

    Return (str):
    Parameters as string, as required by Macs2's snakemake-wrapper
    """
    results: str = config.get("params" ,{}).get("macs2", {}).get("callpeak", "")

    sample_data = get_sample_data_as_dict(sample=str(wildcards.sample), samples=samples)
    if sample_data.get("downstream_file", False):
        if " BAMPE " not in results:
            results += " --format BAMPE "

        if "--nomodel" not in results:
            results += " --nomodel "


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


def get_multiqc_report_peak_input(wildcards: snakemake.io.Wildcards, sample: pandas.DataFrame = samples) -> Dict[str, Union[str, List[str]]]:
    """
    Return expected input files for Multiqc peak calling report, 
    according to user-input, and snakemake-wrapper requirements

    Parameters:
    wildcards (snakemake.io.Wildcards): Required for snakemake unpacking function
    samples   (pandas.DataFrame)      : Describe sample names and related paths/genome

    Return (Dict[str, Union[str, List[str]]]):
    Input files dict, as required by MultiQC's snakemake-wrapper
    """
    results: Dict[str, Union[str, List[str]]] = {
        "picard_qc": [], "fasp": [], "macs2": [], "deeptools_coverage": []
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
            f"tmp/deeptools/plotcoverage/{species}.{build}.{release}.{datatype}/{sample}_coverage.raw"
        )
        results["deeptools_coverage"].append(
            f"tmp/deeptools/plotcoverage/{species}.{build}.{release}.{datatype}/{sample}_coverage.metrics"
        )

        downstream_file: Optional[str] = (
            samples[samples.sample_id == sample]
            .to_dict(orient="index")[0]
            .get("downstream_file")
        )
        if downstream_file:
            results["fastp"].append(f"tmp/fastp/report_pe/{sample}.json")
            results["fastp"].append(f"results/QC/report_pe/{sample}.html")
        else:
            results["fastp"].append(f"tmp/fastp/report_se/{sample}.json")
            results["fastp"].append(f"results/QC/report_se/{sample}.html")

    return results
