# Overview
 * This is a general purpose repository for rnaseq analyses tools/pipelines
 * It contains a script to do QC (fastqc) on raw fastq files and produce a .html general reports and associated QC files for each sample.
 * It also contains a pipeline to go from raw *.fastq.gz* files to count data as described below.


## Installation
  * This requires`R`, `trim_galore`, `samtools`, `STAR`, `HTseq`, `PicardTools`, `GATK`


## Basic usage
``` bash
#This is for the fastQC report of raw sequences. I uses functions from `fastqcr`, which is itself a wrapper for [`fastqc`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
#You can knit this directly from `Rstudio` or you can create this document by running 
Rscript -e "rmarkdown::render('fastqc_reports.Rmd',params = list())"
```


``` r
#(rna_pipeline.R currently called directly in R)
params = list(fq.dir='data/fastq',
              trim.dir='out/fastq.trim',
              genomedir = 'data/reference_genome/chr1_small_index/',
              annotation.gtf ='data/reference_genome/gencode.v43.primary_assembly.annotation_small.gtf',
              genomefasta = 'data/reference_genome/chr1_small.fa',
              out.dir = 'out/')


rna_wrapper(fq.dir = params$fq.dir,
                         trim.dir= params$trim.dir,
                         genomedir = params$genomedir,
                         annotation.gtf = params$annotation.gtf,
                         genomefasta = params$genomefasta,
                         out.dir = params$out.dir)
```


## Further information
  * sebastien.renaut.1@ulaval.ca
  * [www.yohanbosselab.com](https://www.yohanbosselab.com/)

