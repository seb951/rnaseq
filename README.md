# Overview
 * This is a general purpose repository for rnaseq analyses tools/pipelines
 * It contains a script to do either **QC** or **counts**, **kallisto**, **cnv**, **fusion**
   * **QC** uses `fastqc` processes raw *fastq.gz* files and produces a .html general reports and associated QC files for each sample.
   * **star** processes raw *.fastq.gz* files and produces a .tsv file of count data.
 * It also contains a pipeline to go from raw *.fastq.gz* files to count data as described below.
 * Rmarkdown files for deconvolution, PCA, etc, survival curves,...


## Installation
  * This requires `R`, `fastqc`, `trim_galore`, `samtools`, `STAR`, `HTseq`, `PicardTools`,`Kallisto`,`GATK` (not yet)
  * Requires some specific `R` libraries as well (see in [R/](R/)).
  * Requires a genome (.gff, .fasta), if the genome is not indexed, `STAR` will do it.
  * Requires a transcriptome (.fa.gz), if the transcriptome is not indexed, `Kallisto` will do it.


## Basic usage
  * (Note that currently, some options are **QC** or **counts** specific)
``` bash
Usage: R/rna_pipeline.R [options] QC/counts/kallisto


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

## EXAMPLES
``` bash
Rscript R/rna_pipeline.R kallisto --reftranscriptome data/reference_transcriptome/chr1.fasta.gz --quantdir out/kallisto --trimdir out/fastq.trim

Rscript R/rna_pipeline.R star \
--cutadapt /mnt/sde/renseb01/miniconda3/envs/Rbase/bin/cutadapt \
--fqdir /mnt/sde/renseb01/Documents/rnaseq/data/fastq/ \
--trimdir /mnt/sde/renseb01/Documents/rnaseq/out/fastq.trim \
--outdir /mnt/sde/renseb01/Documents/rnaseq/out/ \
--threads 30 \
--genomeindexdir data/reference_genome/GRCh38/ \
--annotationgtf data/reference_genome/gencode.v43.basic.annotation.gtf \
--nbfiles all \
--genomefasta data/reference_genome/GRCh38.primary_assembly.genome.fa

```
## TO DO 
  * urgent:
    * verify that a minimal example exists for **QC** and **star**
    * integrate `gatk` and check the validity of `STAR` options.
  * not-so-urgent:
    * unit-test the functions, better handling of parameters

    

## Further information
  * sebastien.renaut.1@ulaval.ca
  * [www.yohanbosselab.com](https://www.yohanbosselab.com/)
  * May 2023

