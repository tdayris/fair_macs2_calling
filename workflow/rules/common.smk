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
) -> str:
    """
    Return expected parameters for Macs2 peak calling, according to user-input,
    and snakemake-wrapper requirements

    Parameters:
    wildcards (snakemake.io.Wildcards): Required for snakemake unpacking function
    samples   (pandas.DataFrame)      : Describe sample names and related paths/genome

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
    results: Dict[str, List[str]] = {"bam": [], "bai": []}
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
        results["bam"].append(
            f"results/{species}.{build}.{release}.{datatype}/Mapping/{sample}.bam"
        )
        results["bai"].append(
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


def get_multiqc_macs2_peakcalling_report_input(
    wildcards: snakemake.io.Wildcards, sample: pandas.DataFrame = samples
) -> Dict[str, Union[str, List[str]]]:
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
        "picard_qc": [],
        "fasp": [],
        "macs2": [],
        "deeptools_coverage": [],
        "deeptools_fingerprint": [],
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
            f"tmp/deeptools/plotcoverage/{species}.{build}.{release}.{datatype}/Coverage.raw"
        )
        results["deeptools_coverage"].append(
            f"tmp/deeptools/plotcoverage/{species}.{build}.{release}.{datatype}/Coverage.metrics"
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
            results["fastp"].append(f"tmp/fastp/report_pe/{sample}.json")
            results["fastp"].append(f"results/QC/report_pe/{sample}.html")
        else:
            results["fastp"].append(f"tmp/fastp/report_se/{sample}.json")
            results["fastp"].append(f"results/QC/report_se/{sample}.html")

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
            f"results/{species}.{build}.{release}.dna/PeakCalling/{sample}.narrowPeak.bed"
        )

        results["homer"].append(
            f"results/{species}.{build}.{release}.dna/PeakCalling/{sample}.broadPeak.bed"
        )

        results["bedtools"].append(
            f"results/{species}.{build}.{release}.dna/PlotCoverage.png"
        )

        results["bedtools"].append(
            f"results/{species}.{build}.{release}.dna/PlotFingerprint.png"
        )

    return results
