source('R/rnaseq_functions.R')


params = list(fq.dir='/mnt/c/Users/renseb01/Documents/rnaseq/data/fastq',
              trim.dir='/mnt/c/Users/renseb01/Documents/rnaseq/out/fastq.trim',
              genomedir = '/mnt/c/Users/renseb01/Documents/rnaseq/data/reference_genome/chr1_small_index/',
              annotation.gtf ='/mnt/c/Users/renseb01/Documents/rnaseq/data/reference_genome/gencode.v43.primary_assembly.annotation_small.gtf',
              genomefasta = '/mnt/c/Users/renseb01/Documents/rnaseq/data/reference_genome/chr1_small.fa',
              out.dir = 'out/')


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



