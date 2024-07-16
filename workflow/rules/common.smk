import csv
import os
import pandas
import snakemake
import snakemake.utils

from collections import defaultdict
from typing import Any, NamedTuple


snakemake_min_version: str = "8.14.0"
snakemake.utils.min_version(snakemake_min_version)

snakemake_docker_image: str = "docker://snakemake/snakemake:v8.14.0"


container: snakemake_docker_image


# Load and check configuration file
default_config_file: str = "config/config.yaml"


configfile: default_config_file


snakemake.utils.validate(config, "../schemas/config.schema.yaml")


# Load and check samples properties table
def load_table(path: str) -> pandas.DataFrame:
    """
    Load a table in memory, automatically inferring column separators

    Parameters:
    path (str): Path to the table to be loaded

    Return
    (pandas.DataFrame): The loaded table
    """
    with open(path, "r") as table_stream:
        dialect: csv.Dialect = csv.Sniffer().sniff(table_stream.readline())
        table_stream.seek(0)

    # Load table
    table: pandas.DataFrame = pandas.read_csv(
        path,
        sep=dialect.delimiter,
        header=0,
        index_col=None,
        comment="#",
        dtype=str,
    )

    # Remove empty lines
    table = table.where(table.notnull(), None)

    return table


def load_genomes(
    path: str | None = None, samples: pandas.DataFrame | None = None
) -> pandas.DataFrame:
    """
    Load genome file, build it if genome file is missing and samples is not None.

    Parameters:
    path    (str)               : Path to genome file
    samples (pandas.DataFrame)  : Loaded samples
    """
    if path is not None:
        genomes: pandas.DataFrame = load_table(path)

        if samples is not None:
            genomes = used_genomes(genomes, samples)
        return genomes

    elif samples is not None:
        return samples[["species", "build", "release"]].drop_duplicates(
            ignore_index=True
        )

    raise ValueError(
        "Provide either a path to a genome file, or a loaded samples table"
    )


def used_genomes(
    genomes: pandas.DataFrame, samples: pandas.DataFrame | None = None
) -> tuple[str]:
    """
    Reduce the number of genomes to download to the strict minimum
    """
    if samples is None:
        return genomes

    return genomes.loc[
        genomes.species.isin(samples.species.tolist())
        & genomes.build.isin(samples.build.tolist())
        & genomes.release.isin(samples.release.tolist())
    ]


# Load and check samples properties tables
try:
    if (samples is None) or samples.empty():
        sample_table_path: str = config.get("samples", "config/samples.csv")
        samples: pandas.DataFrame = load_table(sample_table_path)
except NameError:
    sample_table_path: str = config.get("samples", "config/samples.csv")
    samples: pandas.DataFrame = load_table(sample_table_path)

snakemake.utils.validate(samples, "../schemas/samples.schema.yaml")


# Load and check genomes properties table
genomes_table_path: str = config.get("genomes", "config/genomes.csv")
try:
    if (genomes is None) or genomes.empty:
        genomes: pandas.DataFrame = load_genomes(genomes_table_path, samples)
except NameError:
    genomes: pandas.DataFrame = load_genomes(genomes_table_path, samples)

snakemake.utils.validate(genomes, "../schemas/genomes.schema.yaml")


report: "../report/workflow.rst"


snakemake_wrappers_prefix: str = "v3.12.0"
release_tuple: tuple[str] = tuple(set(genomes.release.tolist()))
build_tuple: tuple[str] = tuple(set(genomes.build.tolist()))
species_tuple: tuple[str] = tuple(set(genomes.species.tolist()))
datatype_tuple: tuple[str] = ("dna", "cdna", "all", "transcripts")
gxf_tuple: tuple[str] = ("gtf", "gff3")
id2name_tuple: tuple[str] = ("t2g", "id_to_gene")
tmp: str = f"{os.getcwd()}/tmp"
samples_id_tuple: tuple[str] = tuple(samples.sample_id)
stream_tuple: tuple[str] = ("1", "2")
macs2_peak_types: tuple[str] = ("broadPeak", "narrowPeak")


wildcard_constraints:
    macs2_peak_type=r"|".join(macs2_peak_types),
    sample=r"|".join(samples_id_tuple),
    release=r"|".join(release_tuple),
    build=r"|".join(build_tuple),
    species=r"|".join(species_tuple),
    datatype=r"|".join(datatype_tuple),
    stream=r"|".join(stream_tuple),


def lookup_config(
    dpath: str, default: str | None = None, config: dict[str, Any] = config
) -> str:
    """
    Run lookup function with default parameters in order to search a key in configuration and return a default value
    """
    value: str | None = default

    try:
        value = lookup(dpath=dpath, within=config)
    except LookupError:
        value = default
    except WorkflowError:
        value = default

    return value


def lookup_genomes(
    wildcards: snakemake.io.Wildcards,
    key: str,
    default: str | list[str] | None = None,
    genomes: pandas.DataFrame = genomes,
) -> str:
    """
    Run lookup function with default parameters in order to search user-provided sequence/annotation files
    """
    query: str = (
        "species == '{wildcards.species}' & build == '{wildcards.build}' & release == '{wildcards.release}'".format(
            wildcards=wildcards
        )
    )

    query_result: str | float = getattr(
        lookup(query=query, within=genomes), key, default
    )
    if (query_result != query_result) or (query_result is None):
        # Then the result of the query is nan
        return default
    return query_result


def lookup_samples(
    wildcards: snakemake.io.Wildcards,
    key: str,
    default: str | list[str] | None = None,
    samples: pandas.DataFrame = samples,
) -> str | None:
    """
    Run lookup function with default parameters in order to search
    user-proveded sample information
    """
    query = str(
        "species == '{wildcards.species}' & build == '{wildcards.build}' & release == '{wildcards.release}' & sample_id == '{wildcards.sample}'".format(
            wildcards=wildcards
        )
    )
    query_result: str | float = getattr(
        lookup(query=query, within=samples), key, default
    )
    if (query_result != query_result) or (query_result is None):
        # Then the result of the query is nan
        return default
    return query_result


def get_dna_fasta(
    wildcards: snakemake.io.Wildcards, genomes: pandas.DataFrame = genomes
) -> str:
    """
    Return path to the final DNA fasta sequences
    """
    default = str(
        "reference/sequences/{wildcards.species}.{wildcards.build}.{wildcards.release}/{wildcards.species}.{wildcards.build}.{wildcards.release}.dna.fasta".format(
            wildcards=wildcards
        )
    )
    return lookup_genomes(wildcards, key="dna_fasta", default=default, genomes=genomes)


def get_cdna_fasta(
    wildcards: snakemake.io.Wildcards, genomes: pandas.DataFrame = genomes
) -> str:
    """
    Return path to the final cDNA fasta sequences
    """
    default = str(
        "reference/sequences/{wildcards.species}.{wildcards.build}.{wildcards.release}/{wildcards.species}.{wildcards.build}.{wildcards.release}.cdna.fasta".format(
            wildcards=wildcards
        )
    )
    return lookup_genomes(wildcards, key="cdna_fasta", default=default, genomes=genomes)


def get_transcripts_fasta(
    wildcards: snakemake.io.Wildcards, genomes: pandas.DataFrame = genomes
) -> str:
    """
    Return path to the final cDNA transcripts fasta sequences
    """
    default = str(
        "reference/sequences/{wildcards.species}.{wildcards.build}.{wildcards.release}/{wildcards.species}.{wildcards.build}.{wildcards.release}.transcripts.fasta".format(
            wildcards=wildcards
        )
    )
    return lookup_genomes(
        wildcards, key="transcripts_fasta", default=default, genomes=genomes
    )


def select_fasta(
    wildcards: snakemake.io.Wildcards, genomes: pandas.DataFrame = genomes
) -> str:
    """
    Evaluates the {datatype} wildcard, and return the right fasta file
    """
    return branch(
        condition=str(wildcards.datatype).lower(),
        cases={
            "dna": get_dna_fasta(wildcards),
            "cdna": get_cdna_fasta(wildcards),
            "transcripts": get_transcripts_fasta(wildcards),
        },
    )


def get_dna_fai(
    wildcards: snakemake.io.Wildcards, genomes: pandas.DataFrame = genomes
) -> str:
    """
    Return path to the final DNA fasta sequences index
    """
    default = str(
        "reference/sequences/{wildcards.species}.{wildcards.build}.{wildcards.release}/{wildcards.species}.{wildcards.build}.{wildcards.release}.dna.fasta.fai".format(
            wildcards=wildcards
        )
    )
    return lookup_genomes(wildcards, key="dna_fai", default=default, genomes=genomes)


def get_cdna_fai(
    wildcards: snakemake.io.Wildcards, genomes: pandas.DataFrame = genomes
) -> str:
    """
    Return path to the final cDNA fasta sequences index
    """
    default = str(
        "reference/sequences/{wildcards.species}.{wildcards.build}.{wildcards.release}/{wildcards.species}.{wildcards.build}.{wildcards.release}.cdna.fasta.fai".format(
            wildcards=wildcards
        )
    )
    return lookup_genomes(wildcards, key="cdna_fai", default=default, genomes=genomes)


def get_transcripts_fai(
    wildcards: snakemake.io.Wildcards, genomes: pandas.DataFrame = genomes
) -> str:
    """
    Return path to the final cDNA transcripts fasta sequences index
    """
    default = str(
        "reference/sequences/{wildcards.species}.{wildcards.build}.{wildcards.release}/{wildcards.species}.{wildcards.build}.{wildcards.release}.transcripts.fasta.fai".format(
            wildcards=wildcards
        )
    )
    return lookup_genomes(
        wildcards, key="transcripts_fai", default=default, genomes=genomes
    )


def get_dna_dict(
    wildcards: snakemake.io.Wildcards, genomes: pandas.DataFrame = genomes
) -> str:
    default = str(
        "reference/sequences/{wildcards.species}.{wildcards.build}.{wildcards.release}/{wildcards.species}.{wildcards.build}.{wildcards.release}.dna.dict".format(
            wildcards=wildcards
        )
    )
    return lookup_genomes(wildcards, key="dna_dict", default=default, genomes=genomes)


def get_cdna_dict(
    wildcards: snakemake.io.Wildcards, genomes: pandas.DataFrame = genomes
) -> str:
    default = str(
        "reference/sequences/{wildcards.species}.{wildcards.build}.{wildcards.release}/{wildcards.species}.{wildcards.build}.{wildcards.release}.cdna.dict".format(
            wildcards=wildcards
        )
    )
    return lookup_genomes(wildcards, key="cdna_dict", default=default, genomes=genomes)


def get_transcripts_dict(
    wildcards: snakemake.io.Wildcards, genomes: pandas.DataFrame = genomes
) -> str:
    default = str(
        "reference/sequences/{wildcards.species}.{wildcards.build}.{wildcards.release}/{wildcards.species}.{wildcards.build}.{wildcards.release}.transcripts.dict".format(
            wildcards=wildcards
        )
    )
    return lookup_genomes(
        wildcards, key="transcripts_dict", default=default, genomes=genomes
    )


def select_dict(
    wildcards: snakemake.io.Wildcards, genomes: pandas.DataFrame = genomes
) -> str:
    """
    Evaluates the {datatype} wildcards and return the right fasta dictionary
    """
    return branch(
        condition=str(wildcards.datatype).lower(),
        cases={
            "dna": get_dna_dict(wildcards),
            "cdna": get_cdna_dict(wildcards),
            "transcripts": get_transcripts_dict(wildcards),
        },
    )


def select_fai(
    wildcards: snakemake.io.Wildcards, genomes: pandas.DataFrame = genomes
) -> str:
    """
    Evaluates the {datatype} wildcard, and return the right fasta index file
    """
    return branch(
        condition=str(wildcards.datatype).lower(),
        cases={
            "dna": get_dna_fai(wildcards),
            "cdna": get_cdna_fai(wildcards),
            "transcripts": get_transcripts_fai(wildcards),
        },
    )


def get_gtf(
    wildcards: snakemake.io.Wildcards, genomes: pandas.DataFrame = genomes
) -> str:
    """
    Return path to the final genome annotation
    """
    default = str(
        "reference/annotation/{wildcards.species}.{wildcards.build}.{wildcards.release}/{wildcards.species}.{wildcards.build}.{wildcards.release}.gtf".format(
            wildcards=wildcards
        )
    )
    return lookup_genomes(wildcards, key="gtf", default=default, genomes=genomes)


def get_intervals(
    wildcards: snakemake.io.Wildcards, genomes: pandas.DataFrame = genomes
) -> str | None:
    """
    Return path to capturekit file
    """
    return lookup_genomes(wildcards, key="capture_kit", default=None, genomes=genomes)


def get_effective_genome_size(
    wildcards: snakemake.io.Wildcards, genomes: pandas.DataFrame = genomes
) -> str | int:
    """
    Return effective genome size if available in genomes table

    Parameters:
    wildcards (snakemake.io.Wildcards): Required for snakemake unpacking function
    genomes   (pandas.DataFrame)      : Genome description and reference file paths

    Return    (str | None)
    effective genome size if available else None
    """
    egs: int | None = lookup_genomes(
        wildcards, key="effective_genome_size", default=None
    )
    if not egs:
        return {
            "GRCh37": 2864785220,
            "GRCh38": 2913022398,
            "GRCm37": 2620345972,
            "GRCm38": 2652783500,
            "GRCm39": 2654621783,
        }.get(str(wildcards.release), 2652783500)
    return egs


def get_read_length(
    wildcards: snakemake.io.Wildcards, samples: pandas.DataFrame = samples
) -> int:
    """
    Return read length

    Parameters:

    wildcards (snakemake.io.Wildcards): Required for snakemake unpacking function
    samples   (pandas.DataFrame)      : Samples description and their file paths

    Return    (int): Read length (default 100)
    """
    return int(lookup_samples(wildcards, key="read_length", default=100))


def get_blacklist(
    wildcards: snakemake.io.Wildcards, genomes: pandas.DataFrame = genomes
) -> str | None:
    """
    Return blacklist file, if any.

    Parameters:
    wildcards (snakemake.io.Wildcards): Required for snakemake unpacking function
    genomes   (pandas.DataFrame)      : Genome description and reference file paths

    Return: str | None
    If blacklist exists, then return it. Else, return None
    """
    blacklist: str | None = lookup_genomes(wildcards, key="blacklist", default=[])
    if (not blacklist) and (
        str(wildcards.species).lower in ("homo_sapiens", "mus_musculus")
    ):
        return f"reference/blacklist/{wildcards.species}.{wildcards.build}.{wildcards.release}.merged.bed"
    return blacklist


def get_corresponding_input(
    wildcards: snakemake.io.Wildcards, samples: pandas.DataFrame = samples
) -> str | None:
    """
    Return the sample_id of the input file corresponding to the given sample's wildcards
    """
    return lookup_samples(wildcards, key="input_id", default=None)


def get_aln(
    wildcards: snakemake.io.Wildcards,
    index: bool = False,
    sample_id: str | None = None,
    samples: pandas.DataFrame = samples,
    config: dict[str, Any] = config,
) -> str:
    """
    Return the path to the correct bam file according to
    user configuration
    """
    use_sieve: bool = lookup_config(
        dpath="params/make_sieve", default=False, config=config
    )
    default: str = ""
    ext: str = ".bam"
    if index is True:
        ext += ".bai"

    if sample_id is None:
        sample_id = str(wildcards.sample)

    if use_sieve is True:
        return "tmp/fair_bowtie2_mapping_deeptools_alignment_sieve/{wildcards.species}.{wildcards.build}.{wildcards.release}.dna/{sample_id}{ext}".format(
            wildcards=wildcards, ext=ext, sample_id=sample_id
        )

    return "results/{wildcards.species}.{wildcards.build}.{wildcards.release}.dna/Mapping/{sample_id}{ext}".format(
        wildcards=wildcards, ext=ext, sample_id=sample_id
    )


def get_macs2_callpeak_input(
    wildcards: snakemake.io.Wildcards,
    samples: pandas.DataFrame = samples,
    config: dict[str, Any] = config,
) -> dict[str, str]:
    """
    Return expected input files for Macs2 peak calling, according to user-input,
    and snakemake-wrapper requirements

    Parameters:
    wildcards (snakemake.io.Wildcards): Required for snakemake unpacking function
    samples   (pandas.DataFrame)      : Describe sample names and related paths/genome
    config    (dict[str, Any])        : User configuration

    Return (dict[str, str]):
    Dictionnary of all input files as required by Macs2's snakemake-wrapper
    """
    macs2_callpeak_input: dict[str, str] = {
        "treatment": get_aln(wildcards, samples=samples, config=config),
        "treatment_index": get_aln(
            wildcards, index=True, samples=samples, config=config
        ),
    }

    input_id: str | None = get_corresponding_input(wildcards)
    if input_id:
        macs2_callpeak_input["control"] = get_aln(
            wildcards, sample_id=input_id, samples=samples, config=config
        )
        macs2_callpeak_input["control_index"] = get_aln(
            wildcards, sample_id=input_id, index=True, samples=samples, config=config
        )

    return macs2_callpeak_input


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
    macs2_parameters: str = lookup_config(
        dpath=f"params/fair_macs2_calling_macs2_callpeak", default=""
    )

    sample_data: str | bool = lookup_samples(
        wildcards, key="downstream_file", default=False
    )
    if sample_data:
        if " BAMPE " not in macs2_parameters:
            macs2_parameters += " --format BAMPE "

        if "--nomodel" not in macs2_parameters:
            macs2_parameters += " --nomodel "
    elif "--nolambda" not in macs2_parameters:
        macs2_parameters += " --nolambda "

    species: str = str(wildcards.species)
    if species.lower() == "homo_sapiens" and " hs " not in macs2_parameters:
        macs2_parameters += " --gsize hs "
    elif species.lower() == "mus_musculus" and " mm " not in macs2_parameters:
        macs2_parameters += " --gsize mm "
    elif (
        species.lower() == "drosophila_melanogaster" and " dm " not in macs2_parameters
    ):
        macs2_parameters += " --gsize dm "
    elif species.lower() == "caenorhabditis_elegans" and " ce " not in macs2_parameters:
        macs2_parameters += " --gsize ce "

    return macs2_parameters


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
    homer_params: str = lookup_config(
        dpath="params/fair_macs2_calling_homer_annotate_peaks", default=""
    )

    fragment_size: str | None = lookup_samples(wildcards, "fragment_size")
    downstream_path: str | None = lookup_samples(wildcards, "downstream_file")

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

    parameters: dict[str, str] = {
        "extra": homer_params,
        "mode": mode,
    }
    return parameters


def get_deeptools_plotcoverage_input(
    wildcards: snakemake.io.Wildcards,
    samples: pandas.DataFrame = samples,
    genomes: pandas.DataFrame = genomes,
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
    genome_restricted_samples = lookup(
        query=f"species == '{wildcards.species}' & release == '{wildcards.release}' & build == '{wildcards.build}'",
        within=samples,
    )
    results: dict[str, list[str]] = {
        "bams": collect(
            "results/{sample.species}.{sample.build}.{sample.release}.dna/Mapping/{sample.sample_id}.bam",
            sample=genome_restricted_samples,
        ),
        "bais": collect(
            "results/{sample.species}.{sample.build}.{sample.release}.dna/Mapping/{sample.sample_id}.bam.bai",
            sample=genome_restricted_samples,
        ),
    }

    blacklist: str | None = get_blacklist(wildcards, genomes=genomes)
    if blacklist:
        results["blacklist"] = blacklist

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
    macs2_peak_type: str = str(wildcards.macs2_peak_type)

    results: dict[str, list[str] | str] = {
        "bed": f"tmp/fair_macs2_calling_bedtools_merge_macs2_sorted_peaks/{species}.{build}.{release}/{macs2_peak_type}.merged.bed",
        "bw": expand(
            "results/{sample.species}.{sample.build}.{sample.release}.dna/Coverage/{sample.sample_id}.bw",
            sample=lookup(
                query=f"species == '{species}' & release == '{release}' & build == '{build}'",
                within=samples,
            ),
        ),
    }

    blacklist: str | None = get_blacklist(wildcards, genomes)
    if blacklist:
        results["blacklist"] = blacklist

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
            f"tmp/fair_macs2_calling_summarize_homer/{species}.{build}.{release}.{datatype}/{sample}.{macs2_peak_type}.json"
        )

    return results


def get_single_ended_samples(samples: pandas.DataFrame = samples) -> NamedTuple:
    """
    Return the list of single ended samples, as a NameTuple usable
    in expand/collect function

    Parameters:
    samples (pandas.DataFrame): Sample information

    Return: NamedTuple of single ended samples
    """
    return lookup(query="downstream_file != downstream_file", within=samples)


def get_pair_ended_samples(samples: pandas.DataFrame = samples) -> NamedTuple:
    """
    Return the list of pair ended samples, as a NameTuple usable
    in expand/collect function

    Parameters:
    samples (pandas.DataFrame): Sample information

    Return: NamedTuple of pair ended samples
    """
    return lookup(query="downstream_file == downstream_file", within=samples)


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
        "multiqc": [
            "results/QC/MultiQC_FastQC.html",
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

        results["multiqc"].append(
            f"results/{species}.{build}.{release}.dna/QC/MultiQC_PeakCalling.html",
        )

    results["multiqc"] = list(set(results["multiqc"]))

    return results
