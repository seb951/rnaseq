library(optparse)
source('R/rnaseq_functions.R')


#For info, in Windows, you can call Rscript as such: C:'\Program Files\'R\R-4.2.2\bin\x64\Rscript.exe script.R

#make options 
option_list = list(
  make_option( "--fqdir",type="character", default='data/fastq', 
               help="fastq directory [default %default]", metavar="character"),
  
  make_option( "--trimdir",type="character", default='out/fastq.trim', 
               help="trimmed fastq directory [default %default]", metavar="character"),
  
  make_option( "--genomedir",type="character", default='data/reference_genome/chr1_small_index', 
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
               help="Nb of Threads to use [default %default]", metavar="character")

) 

opt_parser = OptionParser(option_list=option_list)

params = parse_args(opt_parser)


#running function
rna_wrapper(fq.dir = params$fqdir,
                         trim.dir= params$trimdir,
                         genomedir = params$genomedir,
                         annotation.gtf = params$annotationgtf,
                         genomefasta = params$genomefasta,
                         out.dir = params$outdir,
                         cutadapt=params$cutadapt,
                        threads= params$threads)


