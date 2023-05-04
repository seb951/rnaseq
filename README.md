# Overview
 * This is a general purpose repository for rnaseq analyses tools/pipelines
 * It contains a script to do either **QC** or **counts**.
   * **QC** uses `fastqc` processes raw *fastq.gz* files and produces a .html general reports and associated QC files for each sample.
   * **counts** processes raw *.fastq.gz* files and produces a .tsv file of count data.
 * It also contains a pipeline to go from raw *.fastq.gz* files to count data as described below.


## Installation
  * This requires `R`, `fastqc`, `trim_galore`, `samtools`, `STAR`, `HTseq`, `PicardTools`, `GATK` (not yet)
  * Requires some specific `R` libraries as well (see in [R/](R/)).
  * Requires a genome (.gff, .fasta), if the genome is not indexed, `STAR` will do it.


## Basic usage
  * (Note that currently, some options are **QC** or **counts** specific)
``` bash
Usage: R/rna_pipeline.R [options] QC/counts


Options:
        --fqdir=CHARACTER
                fastq directory [default data/fastq]

        --trimdir=CHARACTER
                trimmed fastq directory [default out/fastq.trim]

        --genomeindexdir=CHARACTER
                reference genome index directory [default data/reference_genome/chr1_small_index]

        --annotationgtf=CHARACTER
                reference genome gtf file [default data/reference_genome/gencode.v43.primary_assembly.annotation_small.gtf]

        --genomefasta=CHARACTER
                reference genome fasta file [default data/reference_genome/chr1_small.fa]

        --outdir=CHARACTER
                output directory [default out/]

        --cutadapt=CHARACTER
                cutadapt directory [default /usr/bin/cutadapt]

        --threads=CHARACTER
                Nb of Threads to use [default 12]

        --nbfiles=CHARACTER
                specify how many files to process [default all] or an integer or two integers seperated by a comma (e.g. 2,4)

        --qcdir=CHARACTER
                QC argument: output directory [default data/fastqc_results]

        --metadata=CHARACTER
                QC argument: output directory [default data/librairies_1_a_70_et_RIN.xlsx]

        --fastqc=CHARACTER
                QC argument: where is fastqc installed [default data/fastqc]

        -h, --help
                Show this help message and exit

```



## Further information
  * sebastien.renaut.1@ulaval.ca
  * [www.yohanbosselab.com](https://www.yohanbosselab.com/)
  * May 2023

