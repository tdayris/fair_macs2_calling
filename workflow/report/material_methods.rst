Matierial and methods
=====================

Genome information was download from Ensembl. Samtools_ [#samtoolspaper]_ 
and Picard_ [#picardpaper]_ were used to index genome sequences.
Agat_ [#agatpaper]_ was used to correct common issues found in Ensembl
genome annotation files.

Raw fastq file quality was assessed with FastQC_ [#fastqcpaper]_.
Raw fastq files were trimmed using Fastp_ [#fastppaper]_ . Cleaned reads were aligned 
over indexed Ensembl genome with Bowtie2_ [#bowtie2paper]_. Sambamba_ [#sambambapaper]_ 
was used to sort, filter, mark duplicates, and compress aligned reads. Quality 
controls were done using Picard_ [#picardpaper]_ and Samtools_ [#samtoolspaper]_.

Genome coverage and associated quality controls were produced by DeepTools_ [#deeptoolspaper]_,
peak calling was performed with Macs2_ [#macs2paper]_ and annotated with Homer_ [#homerpaper]_.
Quality repord produced during both trimming and mapping steps have been aggregated with 
MultiQC_ [#multiqcpaper]_. Additional graphs were done with in-house Python_ scripts, using
Seaborn_ [#seabornpaper]_ library. The whole pipeline was powered by Snakemake_ [#snakemakepaper]_.

This pipeline is freely available on Github_, details about installation
usage, and resutls can be found on the `Snakemake workflow`_ page.


.. [#samtoolspaper] Li, Heng, et al. "The sequence alignment/map format and SAMtools." bioinformatics 25.16 (2009): 2078-2079.
.. [#picardpaper] McKenna, Aaron, et al. "The Genome Analysis Toolkit: a MapReduce framework for analyzing next-generation DNA sequencing data." Genome research 20.9 (2010): 1297-1303.
.. [#agatpaper] Dainat J. AGAT: Another Gff Analysis Toolkit to handle annotations in any GTF/GFF format.  (Version v0.7.0). Zenodo. https://www.doi.org/10.5281/zenodo.3552717
.. [#fastqcpaper] Andrews, S. Fastqc. "A quality control tool for high throughput sequence data. Augen, J.(2004). Bioinformatics in the post-genomic era: Genome, transcriptome, proteome, and information-based medicine." (2010).
.. [#fastppaper] Chen, Shifu, et al. "fastp: an ultra-fast all-in-one FASTQ preprocessor." Bioinformatics 34.17 (2018): i884-i890.
.. [#bowtie2paper] Langmead, Ben, and Steven L. Salzberg. "Fast gapped-read alignment with Bowtie 2." Nature methods 9.4 (2012): 357-359.
.. [#sambambapaper] Tarasov, Artem, et al. "Sambamba: fast processing of NGS alignment formats." Bioinformatics 31.12 (2015): 2032-2034.
.. [#deeptoolspaper] Ramírez, Fidel, et al. "deepTools: a flexible platform for exploring deep-sequencing data." Nucleic acids research 42.W1 (2014): W187-W191.
.. [#macs2paper] Zhang, Yong, et al. "Model-based analysis of ChIP-Seq (MACS)." Genome biology 9.9 (2008): 1-9.
.. [#homerpaper] Heinz, Sven, et al. "Simple combinations of lineage-determining transcription factors prime cis-regulatory elements required for macrophage and B cell identities." Molecular cell 38.4 (2010): 576-589.
.. [#multiqcpaper] Ewels, Philip, et al. "MultiQC: summarize analysis results for multiple tools and samples in a single report." Bioinformatics 32.19 (2016): 3047-3048.
.. [#seabornpaper] Waskom, Michael L. "Seaborn: statistical data visualization." Journal of Open Source Software 6.60 (2021): 3021.
.. [#snakemakepaper] Köster, Johannes, and Sven Rahmann. "Snakemake—a scalable bioinformatics workflow engine." Bioinformatics 28.19 (2012): 2520-2522.

.. _Sambamba: https://snakemake-wrappers.readthedocs.io/en/v3.2.0/wrappers/sambamba.html
.. _Bowtie2: https://snakemake-wrappers.readthedocs.io/en/v3.2.0/wrappers/bowtie2.html
.. _Fastp: https://snakemake-wrappers.readthedocs.io/en/v3.2.0/wrappers/fastp.html
.. _Picard: https://snakemake-wrappers.readthedocs.io/en/v3.2.0/wrappers/picard/collectmultiplemetrics.html
.. _MultiQC: https://snakemake-wrappers.readthedocs.io/en/v3.2.0/wrappers/multiqc.html
.. _Snakemake: https://snakemake.readthedocs.io
.. _Github: https://github.com/tdayris/fair_bowtie2_mapping
.. _`Snakemake workflow`: https://snakemake.github.io/snakemake-workflow-catalog?usage=tdayris/fair_bowtie2_mapping
.. _Agat: https://agat.readthedocs.io/en/latest/index.html
.. _Samtools: https://snakemake-wrappers.readthedocs.io/en/v3.2.0/wrappers/samtools/faidx.html
.. _DeepTools: https://snakemake-wrappers.readthedocs.io/en/v3.2.0/wrappers/deeptools.html
.. _Macs2: https://snakemake-wrappers.readthedocs.io/en/v3.2.0/wrappers/macs2/callpeak.html
.. _Homer: https://snakemake-wrappers.readthedocs.io/en/v3.2.0/wrappers/homer/annotatePeaks.html
.. _FastQC: https://snakemake-wrappers.readthedocs.io/en/v3.2.0/wrappers/fastqc.html
.. _Python: docs.python.org
.. _Seaborn: https://seaborn.pydata.org/index.html

:Authors:
    Thibault Dayris

:Version: 1.0.0 of 12/12/2023
