suppressWarnings(library(optparse))
source('R/counts_functions.R')


#For info, in Windows, you can call Rscript as such: C:'\Program Files\'R\R-4.2.2\bin\x64\Rscript.exe script.R

#options (QC and counts)
option_list = list(
  make_option( "--fqdir",type="character", default='data/fastq', 
               help="fastq directory [default %default]", metavar="character"),
  
  make_option( "--trimdir",type="character", default='out/fastq.trim', 
               help="trimmed fastq directory [default %default]", metavar="character"),
  
  make_option( "--genomeindexdir",type="character", default='data/reference_genome/chr1_small_index', 
               help="reference genome index directory [default %default]", metavar="character"),
  
  make_option( "--annotationgtf",type="character", default='data/reference_genome/gencode.v43.primary_assembly.annotation_small.gtf', 
               help="reference genome gtf file [default %default]", metavar="character"),
  
  make_option( "--genomefasta",type="character", default='data/reference_genome/chr1_small.fa', 
               help="reference genome fasta file [default %default]", metavar="character"),
  
  make_option( "--outdir",type="character", default='out/', 
               help="output directory [default %default]", metavar="character"),
  
  make_option( "--cutadapt",type="character", default='/usr/bin/cutadapt', 
               help="cutadapt directory [default %default]", metavar="character"),
  
  make_option( "--threads",type="integer", default=12, 
               help="Nb of Threads to use [default %default]", metavar="character"),
  
  make_option( "--nbfiles",type="character", default='all', 
               help="specify how many files to process [default %default] or an integer or two integers seperated by a comma (e.g. 2,4)", metavar="character"),
  
  ###QC specific inputs
  make_option( "--qcdir",type="character", default='fastqc_results', 
               help="QC argument: output directory [default %default]", metavar="character"),
  
  make_option( "--metadata",type="character", default='data/librairies_1_a_70_et_RIN.xlsx', 
               help="QC argument: output directory [default %default]", metavar="character"),
  
  make_option( "--fastqc",type="character", default='fastqc', 
               help="QC argument: where is fastqc installed [default %default]", metavar="character")
  ) 

#
parser <- OptionParser(usage = "%prog [options] QC/counts", option_list=option_list)
arguments <- parse_args(parser,positional_arguments = 1)

if(arguments$args == 'QC' | arguments$args == 'counts') {
  sprintf("Running command ( %s )", arguments$args)
  } else {
         stop(sprintf("Specified command ( %s ) does not exist", arguments$args))
  }



#running rnaseq wrapper
if(arguments$args == 'counts') {
rna_wrapper(fq.dir = arguments$options$fqdir,
                         trim.dir= arguments$options$trimdir,
                         genomedir = arguments$options$genomeindexdir,
                         annotation.gtf = arguments$options$annotationgtf,
                         genomefasta = arguments$options$genomefasta,
                         out.dir = arguments$options$outdir,
                         cutadapt=arguments$options$cutadapt,
                         threads= arguments$options$threads,
                         nbfiles= arguments$options$nbfiles
            )
}

#running QC notebook
if(arguments$args == 'QC') {
dir.create(arguments$options$outdir)
library(rmarkdown)
rmarkdown::render('./Rmarkdown/fastqc_reports.Rmd',params = list(fq.dir = arguments$options$fqdir,
                                                                 qc.dir = file.path(arguments$options$outdir,arguments$options$qcdir),
                                                                 threads =  arguments$options$threads,
                                                                 metadata = arguments$options$metadata,
                                                                 fastqc.path = arguments$options$fastqc,
                                                                 adapters.dir = NULL),
                  output_dir = arguments$options$outdir
                  )
}


#running another module here:
