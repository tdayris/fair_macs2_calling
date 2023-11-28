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
    ├── XXX.dna
    │   ├── Mapping
    │   |   ├── YYY.bam
    │   |   └── YYY.bam.bai
    |   ├── Coverage
    |   |   └── YYY.bw
    |   ├── PlotCoverage.png
    |   ├── PlotFingerprint.png
    |   └── PeakCalling
    |       └── YYY.bed
    └── QC
        ├── MultiQC_data.zip
        ├── MultiQC.html
        └── report_pe
            └── YYY.html


+------------------+------------------------+-----------------------------------------------+
| Directory        | File Extension         | Content                                       |
+==================+========================+===============================================+
| XXX              | `PlotCoverage.png`     | Read coverage around TSS                      |
+                  +------------------------+-----------------------------------------------+
|                  | `PlotFingerprint.png`  | QC of antigen specificity                     |
+------------------+------------------------+-----------------------------------------------+
| XXX/Mapping      | `YYY.bam`              | Aligned reads                                 |
+                  +------------------------+-----------------------------------------------+
|                  | `YYY.bam.bai`          | Aligned reads index                           |
+------------------+------------------------+-----------------------------------------------+
| XXX/Coverage     | `YYY.bw`               | Read coverage                                 |
+                  +------------------------+-----------------------------------------------+
| XXX/PeakCalling  | `YYY.bed`              | Annotated peaks                               |
+                  +------------------------+-----------------------------------------------+
| QC               | `MultiQC_data.zip`     | Zipped figures and tables                     |
+                  +------------------------+-----------------------------------------------+
|                  | `MultiQC.html`         | Complete quality report, includes all samples |
+------------------+------------------------+-----------------------------------------------+
| QC/report_pe     | `YYY.html`             | Sequence quality report for PE sample `YYY`   |
+------------------+------------------------+-----------------------------------------------+
| QC/report_se     | `YYY.html`             | Sequence quality report for SE sample `YYY`   |
+------------------+------------------------+-----------------------------------------------+