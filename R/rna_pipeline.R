source('rnaseq/R/rnaseq_functions.R')


params = list(fq.dir='/mnt/c/Users/renseb01/Documents/rnaseq/data/fastq',
              trim.dir='/mnt/c/Users/renseb01/Documents/rnaseq/out/fastq.trim',
              genomedir = '/mnt/c/Users/renseb01/Documents/rnaseq/data/reference_genome/chr1_index/',
              annotation.gtf ='/mnt/c/Users/renseb01/Documents/rnaseq/data/reference_genome/gencode.v43.primary_assembly.annotation_small.gtf',
              genomefasta = '/mnt/c/Users/renseb01/Documents/rnaseq/data/reference_genome/chr1.fa',
              out.dir = 'rnaseq/out/')


#fastq files
fastq_files = sequences(fq.dir = params$fq.dir,
                        trim.dir = params$trim.dir)

#out_prefix
out_prefix = paste0(params$out.dir,fastq_files[[5]][1],'_')

#trimming
trimming(trim.dir = params$trim.dir,
         R1 = fastq_files[[1]][1],
         R2 = fastq_files[[2]][1])


#mapping
mapping(genomedir = params$genomedir,
        genomefasta = params$genomefasta,
        threads = 12,
        R1_trim = fastq_files[[3]][1],
        R2_trim = fastq_files[[4]][1],
        annotation.gtf = params$annotation.gtf,
        out_prefix = out_prefix)
        
     
#bamtosam
bamtosam(out_prefix = out_prefix)
         
#picardtools         
picardtools(out_prefix = out_prefix,
            out_prefix_nopath = fastq_files[[5]][1])

#sortbam
sort_bam_by_name(out_prefix = out_prefix)

#htseq
htseq_count(out_prefix = out_prefix,
            annotation.gtf = params$annotation.gtf)

#cleanup
cleanup(out_prefix = out_prefix)
