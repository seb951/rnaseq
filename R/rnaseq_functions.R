

#=====================
###get sequences ready
#=====================
sequences = function(fq.dir='data/fastq',
                     trim.dir='out/fastq.trim',
                     out.dir = 'out/') {
  R12 = list.files(gsub("/mnt/c",ifelse(Sys.info()['sysname'] == 'Windows','C:',''),fq.dir),pattern = 'R[12].fastq.gz$')
  sample_names = unique(gsub('_R[12].fastq.gz','',R12))
  
  R1 = paste0(file.path(fq.dir,sample_names),'_R1.fastq.gz')
  R2 = paste0(file.path(fq.dir,sample_names),'_R2.fastq.gz')

  
  R1_trim = paste0(file.path(trim.dir,sample_names),'_R1_val_1.fq.gz')
  R2_trim = paste0(file.path(trim.dir,sample_names),'_R2_val_2.fq.gz')
  
  
  out.dir.length = length(list.files(gsub("/mnt/c",ifelse(Sys.info()['sysname'] == 'Windows','C:',''),out.dir))) 
  
  if(length(R12)==0) warning('It appears like your directory contains no .fastq.gz file',immediate. = TRUE)
  if(out.dir.length>1) warning(paste0(out.dir,' already contains ',out.dir.length,' files. I will add files to this directory'),immediate. = TRUE)
  message(paste0('Done listing sequences, Time is: ',Sys.time()))
  
  return(list(R1,R2,R1_trim,R2_trim,sample_names))
}

#=====================
###trimming
#=====================
trimming = function(trim.dir = 'out/fastq.trim',
                    R1 = sequences()[[1]][1],
                    R2 = sequences()[[2]][1],
                    out_seq = sequences()[[4]][1],
                    adaptor1='AGATCGGAAGAGCACACGTCTGAACTCCAGTCA',
                    adaptor2='AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT',
                    cutadapt = '/usr/bin/cutadapt',
                    i=1,
                    out.dir=out.dir) {
  
  if(file.exists(gsub("/mnt/c",ifelse(Sys.info()['sysname'] == 'Windows','C:',''),trim.dir))==F) dir.create(gsub("/mnt/c",ifelse(Sys.info()['sysname'] == 'Windows','C:',''),trim.dir))

 
  cmd = paste0(ifelse(Sys.info()['sysname'] == 'Windows','wsl.exe ',''),
               'trim_galore -q 20 -o ',
               trim.dir,
               ' --phred33 --adapter=',adaptor1,' --adapter2=',adaptor2,' --paired ',
               R1,
               ' ',
               R2,
               ' --path_to_cutadapt ',
               cutadapt, ' 1>',ifelse(i>1,'>',''),file.path(out.dir,'logs'),'/trim_galore.out 2>',ifelse(i>1,'>',''),file.path(out.dir,'logs'),'/trim_galore.err')
  
  if(file.exists(gsub("/mnt/c",ifelse(Sys.info()['sysname'] == 'Windows','C:',''),out_seq))==F) {system(cmd);message(paste0('Done Trimming, Time is: ',Sys.time()))} else {message(paste0('Skipping Trimming, Time is: ',Sys.time()))}
  
  
  
  return('')
}


#=====================
#Prepare genome (fasta, gff downloaded from here: https://www.gencodegenes.org/human/)
#=====================
.generate_genome = function(genomefasta = 'data/reference_genome/chr1.fa',
                           annotation.gtf = 'data/reference_genome/gencode.v43.primary_assembly.annotation_small.gtf',
                           genomedir = 'data/reference_genome/chr1_index',
                           threads = 12){
  cmd = paste0(ifelse(Sys.info()['sysname'] == 'Windows','wsl.exe ',''),
                'STAR --genomeSAindexNbases 12 --runThreadN ',threads,' --runMode genomeGenerate --genomeDir ',
               genomedir,
               ' --genomeFastaFiles ',
               genomefasta,
               ' --sjdbGTFfile ',
               annotation.gtf,
               ' --sjdbOverhang 99')
  
  system(cmd)
  
  #message(cmd)
  
  message(paste0('Done generate_genome, Time is: ',Sys.time()))
  
  return('')
}

#=====================
#Mapping step
#=====================
mapping = function(genomedir = 'data/reference_genome/chr1_index/',
                   genomefasta = 'data/reference_genome/chr1.fa',
                   threads = 12,
                   R1_trim = sequences()[[3]][1],
                   R2_trim = sequences()[[4]][1],
                   annotation.gtf = 'data/reference_genome/gencode.v43.primary_assembly.annotation_small.gtf',
                   out_prefix = paste0('rnaseq/out/',sequences()[[5]][1],'_'),
                   i=1,
                   out.dir=out.dir) {
  
  #generate genome if its not there.
  if(length(list.files(gsub("/mnt/c",ifelse(Sys.info()['sysname'] == 'Windows','C:',''),genomedir)))==0) .generate_genome(genomedir = genomedir,
                                                                                                                          genomefasta = genomefasta,
                                                                                                                          annotation.gtf = annotation.gtf,
                                                                                                                          threads = threads)
  
  #STAR cmd
  cmd = paste0(ifelse(Sys.info()['sysname'] == 'Windows','wsl.exe ',''),'STAR --runMode alignReads --genomeDir ',
               genomedir,
               ' --readFilesCommand zcat --limitBAMsortRAM 100000000000 --readFilesIn ',
               R1_trim,
               ' ',
               R2_trim,
               ' --outSAMtype BAM SortedByCoordinate --sjdbOverhang 99 --outFilterMultimapNmax 1 --outReadsUnmapped None --quantMode GeneCounts --sjdbGTFfile ',
               annotation.gtf,
               ' --runThreadN ',
               threads ,
               ' --outFileNamePrefix ', out_prefix,' 1>',ifelse(i>1,'>',''),file.path(out.dir,'logs'),'/star_aligner.out 2>',ifelse(i>1,'>',''),file.path(out.dir,'logs'),'/star_aligner.err')
  
  system(cmd)

 # message(cmd)
  
  message(paste0('Done STAR, Time is: ',Sys.time()))
  
  return('')
}


#======================
#  Bam to sam
#=====================
bamtosam = function(out_prefix = 'rnaseq/out/toto_',
                    bam = paste0(out_prefix,'Aligned.sortedByCoord.out.bam'),
                    bam_out = paste0(out_prefix,'trimmed_Aligned_PP_UM.bam')){

  #cmd1=paste0(ifelse(Sys.info()['sysname'] == 'Windows','wsl.exe ',''),"samtools view -h -f 0x0002 ",bam," |  grep -P '^@|NH:i:1\t' | samtools view -h -b -  >",bam_out)
  cmd2=paste0(ifelse(Sys.info()['sysname'] == 'Windows','wsl.exe ',''),"samtools view -h -o ",gsub(".bam",".sam",bam_out)," ",bam)

  #system(cmd1)
  system(cmd2)
  
  #message(cmd1)
  #message(cmd2)
  
  message(paste0('Done samtools bamtosam, Time is: ',Sys.time()))
  
  return('')
}

#======================
# Add RG, sort, mark 
# duplicates & create 
# index
#======================
picardtools = function(out_prefix = 'rnaseq/out/toto_',
                       bam = paste0(out_prefix,'Aligned.sortedByCoord.out.bam'),
                       bam_added = paste0(out_prefix,'trimmed_Aligned_PP_UM_rgAdded.bam'),
                       bam_rmdup = paste0(out_prefix,'trimmed_Aligned_PP_UM_rgAdded_dup.bam'),
                       metrics = paste0(out_prefix,'trimmed_Aligned_PP_UM_rgAdded_metrics.txt'),
                       i=i,
                       out.dir=out.dir){
  
  out_prefix_nopath = tail(strsplit(out_prefix,"/")[[1]],1)
  
  cmd1 = paste0(ifelse(Sys.info()['sysname'] == 'Windows','wsl.exe ',''),
                'PicardCommandLine AddOrReplaceReadGroups -I ',
                bam,
                ' -O ',
                bam_added,
                ' -SO coordinate -RGID ',out_prefix_nopath,' -RGLB ',out_prefix_nopath,' -RGPL illumina -RGPU 1 -RGSM ',out_prefix_nopath,' 1>',ifelse(i>1,'>',''),file.path(out.dir,'logs'),'/picardtools.out 2>',ifelse(i>1,'>',''),file.path(out.dir,'logs'),'/picardtools.err')
  
  
  cmd2 = paste0(ifelse(Sys.info()['sysname'] == 'Windows','wsl.exe ',''),
                'PicardCommandLine MarkDuplicates -I ',
                bam_added,
                ' -O ',
                bam_rmdup,
                ' -CREATE_INDEX TRUE -VALIDATION_STRINGENCY SILENT -M ',
                metrics,' 1>>',file.path(out.dir,'logs'),'/picardtools.out 2>>',file.path(out.dir,'logs'),'/picardtools.err')
  
  
  system(cmd1)
  system(cmd2)
  
  #message(cmd1)
  #message(cmd2)

  message(paste0('Done picard tools MarkDuplicates, Time is: ',Sys.time()))
  
  return('')
}

#======================
# splitNtrim & reassign
# mapping quality
#======================
#/lustre03/project/6032391/GROUP/bin/gatk-4.1.2.0/gatk SplitNCigarReads -R /lustre03/project/6032391/GROUP/References/GRCh38/Homo_sapiens_assembly38.fasta -I /Users/jerry/Documents/IUCPQ/rnaseq/src/star/806rcbc2881/806rcbc2881_trimmed_Aligned_PP_UM_rgAdded_dup.bam -O /Users/jerry/Documents/IUCPQ/rnaseq/src/star/806rcbc2881/806rcbc2881_trimmed_Aligned_PP_UM_rgAdded_dup_split.bam




#======================
# sort cleaned bam remove unmapped reads   sort bam by name 
#======================
sort_bam_by_name = function (out_prefix = 'rnaseq/out/toto_',
                             bam_rmdup = paste0(out_prefix,'trimmed_Aligned_PP_UM_rgAdded_dup.bam'),
                             bam_sort =  paste0(out_prefix,'trimmed_Aligned_PP_UM_rgAdded_dup_split_sortP.bam'),
                             bam_nounmapped = paste0(out_prefix,     'trimmed_Aligned_PP_UM_rgAdded_dup_split_noUnmapped_sortP.bam'),
                             bam_nounmapped_sort = paste0(out_prefix,'trimmed_Aligned_PP_UM_rgAdded_dup_split_noUnmapped_sortN.bam'))
                             {
  sort_noGATK = paste0(ifelse(Sys.info()['sysname'] == 'Windows','wsl.exe ',''),"samtools sort -o ",bam_sort,' ',bam_rmdup)
  view = paste0(ifelse(Sys.info()['sysname'] == 'Windows','wsl.exe ',''),"samtools view -b -F 4 ",bam_sort,' >',bam_nounmapped)
  index =  paste0(ifelse(Sys.info()['sysname'] == 'Windows','wsl.exe ',''),"samtools index ",bam_nounmapped)
  resort = paste0(ifelse(Sys.info()['sysname'] == 'Windows','wsl.exe ',''),"samtools sort -n -o ",bam_nounmapped_sort," ",bam_nounmapped)
  
  system(sort_noGATK)
  system(view)
  system(index)
  system(resort)
  
  #message(sort_noGATK)
  #message(view)
  #message(index)
  #message(resort)
  
  message(paste0('Done sort_bam_by_name, Time is: ',Sys.time()))
  
  return('')
}
  


#======================
# htseq-count step 
#======================
htseq_count = function (
    annotation.gtf = 'data/reference_genome/gencode.v43.primary_assembly.annotation_small.gtf',
    #htseq_fwd = paste0(out_prefix,'htseq_fwd.txt'),
    #htseq_reverse = paste0(out_prefix,'htseq_reverse.txt'),
    #htseq_onstranded = paste0(out_prefix,'htseq_unstranded.txt'),
    #out_prefix = 'out/toto_',
    out.dir=out.dir,
    threads=12)
                        {
  
     htseq1 = paste0(ifelse(Sys.info()['sysname'] == 'Windows','wsl.exe ',''),'htseq-count --format=bam --mode=intersection-nonempty --stranded=yes -n ',threads,' --idattr=gene_id ',file.path(out.dir,'*sortN.bam'),' ',annotation.gtf, ' >',file.path(out.dir,'htseq_counts.tsv'),' 2>',file.path(out.dir,'logs','htseq.err'))
 #   htseq2 = paste0(ifelse(Sys.info()['sysname'] == 'Windows','wsl.exe ',''),'htseq-count --format=bam --mode=intersection-nonempty --stranded=reverse --idattr=gene_id ',bam_nounmapped_sort,' ',annotation.gtf, ' >',htseq_reverse,' 2>>',file.path(out.dir,'logs'),'/htseq.err')
 #   htseq3 = paste0(ifelse(Sys.info()['sysname'] == 'Windows','wsl.exe ',''),'htseq-count --format=bam --mode=intersection-nonempty --stranded=no --idattr=gene_id ',bam_nounmapped_sort,' ',annotation.gtf, ' >',htseq_onstranded,' 2>>',file.path(out.dir,'logs'),'/htseq.err')

     system(htseq1)
  #  system(htseq2)
  #  system(htseq3)
  
  #  message(htseq1)
  #  message(htseq2)
  #  message(htseq3)
  
  message(paste0('Done htseq-count for all samples, Time is: ',Sys.time()))  
    
  return('')
}


#======================
# cleanup 
#======================
cleanup = function (out_prefix = 'out/toto_')
  {
  message(paste0('Done cleanup (nothing yet), Time is: ',Sys.time()))
  
  return('')
}


#======================
# wrapper 
#======================
rna_wrapper = function(fq.dir = params$fq.dir,
                       trim.dir= params$trim.dir,
                       genomedir = params$genomedir,
                       annotation.gtf = params$annotation.gtf,
                       genomefasta = params$genomefasta,
                       out.dir = params$out.dir,
                       cutadapt = params$cutadapt,
                       threads = params$threads,
                       nbfiles = params$nbfiles){
  
  dir.create(out.dir)
  dir.create(file.path(out.dir,'logs'))
  
  #fastq files
  fastq_files = sequences(fq.dir = fq.dir,
                          trim.dir = trim.dir,
                          out.dir = out.dir)
  
  
  #files to process depends on nbfiles
  files = seq_along(fastq_files[[5]])
  if(nbfiles=='all') files = files else {files = files[1:min(length(files),as.numeric(nbfiles))]}
  
  for(i in files) {
    
    #out_prefix
    out_prefix = paste0(out.dir,fastq_files[[5]][i],'_')
    
    #trimming
    trimming(trim.dir = trim.dir,
             R1 = fastq_files[[1]][i],
             R2 = fastq_files[[2]][i],
             out_seq = fastq_files[[4]][i],
             cutadapt = cutadapt,
             i=i,
             out.dir = out.dir)
    
    
    #mapping
    mapping(genomedir = genomedir,
            genomefasta = genomefasta,
            threads = threads,
            R1_trim = fastq_files[[3]][i],
            R2_trim = fastq_files[[4]][i],
            annotation.gtf = annotation.gtf,
            out_prefix = out_prefix,
            i = i,
            out.dir = out.dir)
    
    
    #bamtosam
    #bamtosam(out_prefix = out_prefix)
    
    #picardtools         
    picardtools(out_prefix = out_prefix,
                i = i,
                out.dir = out.dir)
    
    #sortbam
    sort_bam_by_name(out_prefix = out_prefix)
    
    #message
    message(paste0('--- Done sample, ',fastq_files[[5]][i],', Time is: ',Sys.time(),' ---'))
  }
  
  #htseq for all
  htseq_count(annotation.gtf = annotation.gtf,
              out.dir = out.dir,
              threads = threads)
  
  #final cleanup
  cleanup(out_prefix = out_prefix)
  
}


