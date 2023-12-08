Results
=======


Alongside with the report, you may find directories called `reference`,
and `results`.

Reference
---------

You shall find all genome-related files in it. Considering a genome named `XXX`,
the following files are present:

::

    reference/
    ├── XXX.all.vcf
    ├── XXX.cdna.fasta
    ├── XXX.cdna.fasta.fai
    ├── XXX.dna.dict
    ├── XXX.dna.fasta
    ├── XXX.dna.fasta.fai
    └── XXX.gtf


+---------------+-----------------------------+
| Extension     | Content                     |
+===============+=============================+
| `.gtf`        | Genome annotation           |
+---------------+-----------------------------+
| `.fasta`      | Genome sequences            |
+---------------+-----------------------------+
| `.fasta.fai`  | Genome sequences index      |
+---------------+-----------------------------+
| `.dict`       | Genome sequences dictionary |
+---------------+-----------------------------+
| `.vcf`        | Genome known variations     |
+---------------+-----------------------------+

These files are quite volumous and are not embeded in this HTML page. Please
find them directly on file system.


Results
-------

Given an samples called `YYY` and a genome called `XXX`,
the following files are present:


::

    results/
    ├── XXX
    │   ├── Coverage
    │   │   └── YYY.bw
    │   ├── Graphs
    │   │   ├── PlotCoverage.png
    │   │   └── PlotFingerprint.png
    │   ├── Mapping
    │   │   ├── YYY.bam
    │   │   └── YYY.bam.bai
    │   └── PeakCalling
    │       ├── broadPeak_annotation
    │       │   └── YYY.broadPeak.tsv
    │       ├── broadPeak_bed
    │       │   └── YYY.broadPeak.bed
    │       ├── narrowPeak_annotation
    │       │   └── YYY.narrowPeak.tsv
    │       └── narrowPeak_bed
    │           └── YYY.narrowPeak.bed
    └── QC
        ├── MultiQC_data.zip
        ├── MultiQC.html
        ├── MultiQC.PeakCalling_data.zip
        ├── MultiQC.PeakCalling.html
        └── report_pe
            ├── YYY.1_fastqc.zip
            ├── YYY.1.html
            ├── YYY.2_fastqc.zip
            ├── YYY.2.html
            └── YYY.html



+---------------------------------------+-----------------------------------+-----------------------------------------------+
| Directory                             | File Extension                    | Content                                       |
+=======================================+===================================+===============================================+
| XXX/Graphs                            | `PlotCoverage.png`                | Read coverage around TSS                      |
+                                       +-----------------------------------+-----------------------------------------------+
|                                       | `PlotFingerprint.png`             | QC of antigen specificity                     |
+---------------------------------------+-----------------------------------+-----------------------------------------------+
| XXX/Mapping                           | `YYY.bam`                         | Aligned reads                                 |
+                                       +-----------------------------------+-----------------------------------------------+
|                                       | `YYY.bam.bai`                     | Aligned reads index                           |
+---------------------------------------+-----------------------------------+-----------------------------------------------+
| XXX/Coverage                          | `YYY.bw`                          | Read coverage                                 |
+                                       +-----------------------------------+-----------------------------------------------+
| XXX/PeakCalling/broadPeak_annotation  | `YYY.boradPeak.tsv`               | Annotated peaks, from Macs2 Braod peaks       |
+---------------------------------------+-----------------------------------+-----------------------------------------------+
| XXX/PeakCalling/broadPeak_bed         | `YYY.boradPeak.tsv`               | Called peaks, from Macs2 Braod peaks (IGV)    |
+---------------------------------------+-----------------------------------+-----------------------------------------------+
| XXX/PeakCalling/narrowPeak_annotation | `YYY.boradPeak.tsv`               | Annotated peaks, from Macs2 Narrow peaks      |
+---------------------------------------+-----------------------------------+-----------------------------------------------+
| XXX/PeakCalling/narrowPeak_bed        | `YYY.boradPeak.tsv`               | Called peaks, from Macs2 Narrow peaks (IGV)   |
+---------------------------------------+-----------------------------------+-----------------------------------------------+
| QC                                    | `MultiQC_data.zip`                | Zipped figures and tables                     |
+                                       +-----------------------------------+-----------------------------------------------+
|                                       | `MultiQC.html`                    | Complete quality report, includes all samples |
+                                       +-----------------------------------+-----------------------------------------------+
|                                       | `MultiQC.PeakCalling.html`        | Complete quality report, includes all samples |
+                                       +-----------------------------------+-----------------------------------------------+
|                                       | `MultiQC.PeakCalling._data.zip`   | Complete quality report, includes all samples |
+---------------------------------------+-----------------------------------+-----------------------------------------------+
| QC/report_pe                          | `YYY.html`                        | Sequence quality report for PE sample `YYY`   |
+---------------------------------------+-----------------------------------+-----------------------------------------------+
| QC/report_se                          | `YYY.html`                        | Sequence quality report for SE sample `YYY`   |
+---------------------------------------+-----------------------------------+-----------------------------------------------+