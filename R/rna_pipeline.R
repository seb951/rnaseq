source('R/rnaseq_functions.R')

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



