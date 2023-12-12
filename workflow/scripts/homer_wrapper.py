# coding: utf-8

__author__ = "Thibault Dayris"
__copyright__ = "Copyright 2023, Thibault Dayris"
__email__ = "thibault.dayris@gustaveroussy.fr"
__license__ = "MIT"


from snakemake.shell import shell

# Optional parameters
log: str = snakemake.log_fmt_shell(stdout=False, stderr=True)
extra: str = snakemake.params.get("extra", "")


# Input optional files
gtf: str = snakemake.input.get("gtf", "")
if gtf:
    gtf = f"-gtf {gtf}"


# Only one coverage track is required here
wig: str = snakemake.input.get("wig", "")
bg: str = snakemake.input.get("bg", "")
if wig:
    wig = f"-wig {wig}"
elif bg:
    bg = f"-bedGraph {bg}"

genome: str = snakemake.input.get("genome", "")
if not genome:
    genome = "hg38"

motifs: str = snakemake.input.get("motifs", "")
if motifs:
    motifs = f"-m {motifs}"


# Build and launch command line
shell(
    "annotatePeaks.pl "
    "{snakemake.input.peaks} "
    "{genome} {wig} {bg} {gtf} {motifs} "
    "-cpu {snakemake.threads} {extra} "
    "> {snakemake.output.annotations} "
    "{log} "
)
