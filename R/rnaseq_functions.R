

#=====================
###get sequences ready
#=====================
sequences = function(fq.dir='/mnt/c/Users/renseb01/Documents/rnaseq/data/fastq',
                     trim.dir='/mnt/c/Users/renseb01/Documents/rnaseq/out/fastq.trim') {
  R12 = list.files(gsub("/mnt/c",ifelse(Sys.info()['sysname'] == 'Windows','C:',''),fq.dir),pattern = 'R[12].fastq.gz')
  sample_names = unique(gsub('_R[12].fastq.gz','',R12))
  
  R1 = paste0(file.path(fq.dir,sample_names),'_R1.fastq.gz')
  R2 = paste0(file.path(fq.dir,sample_names),'_R2.fastq.gz')

  
  R1_trim = paste0(file.path(trim.dir,sample_names),'_R1_val_1.fq.gz')
  R2_trim = paste0(file.path(trim.dir,sample_names),'_R2_val_2.fq.gz')
  
  return(list(R1,R2,R1_trim,R2_trim,sample_names))
}

#=====================
###trimming
#=====================
trimming = function(trim.dir = '/mnt/c/Users/renseb01/Documents/rnaseq/out/fastq.trim',
                    R1 = sequences()[[1]][1],
                    R2 = sequences()[[2]][1],
                    adaptor1='AGATCGGAAGAGCACACGTCTGAACTCCAGTCA',
                    adaptor2='AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT',
                    cutadapt = '/usr/bin/cutadapt') {
  
  if(file.exists(trim.dir)==F) dir.create(gsub("/mnt/c",ifelse(Sys.info()['sysname'] == 'Windows','C:',''),trim.dir))

 
  cmd = paste0(ifelse(Sys.info()['sysname'] == 'Windows','wsl.exe ',''),
               'trim_galore -q 20 -o ',
               trim.dir,
               ' --phred33 --adapter=',adaptor1,' --adapter2=',adaptor2,' --paired ',
               R1,
               ' ',
               R2,
               ' --path_to_cutadapt ',
               cutadapt, ' 1>>trim_galore.out 2>>trim_galore.err')
  
  message(cmd)
  system(cmd)

return('Done Trimming')
}


#=====================
#Prepare genome (fasta, gff downloaded from here: https://www.gencodegenes.org/human/)
#=====================
.generate_genome = function(genomefasta = '/mnt/c/Users/renseb01/Documents/rnaseq/data/reference_genome/chr1.fa',
                           annotation.gtf = '/mnt/c/Users/renseb01/Documents/rnaseq/data/reference_genome/gencode.v43.primary_assembly.annotation_small.gtf',
                           genomedir = '/mnt/c/Users/renseb01/Documents/rnaseq/data/reference_genome/chr1_index'){
  cmd = paste0(ifelse(Sys.info()['sysname'] == 'Windows','wsl.exe ',''),
                'STAR --runThreadN 6 --runMode genomeGenerate --genomeDir ',
               genomedir,
               ' --genomeFastaFiles ',
               genomefasta,
               ' --sjdbGTFfile ',
               annotation.gtf,
               ' --sjdbOverhang 99 1>>star_aligner.out 2>>star_aligner.err')
  
  message(cmd)
  system(cmd)
  return('Done generate_genome')
}

#=====================
#Mapping step
#=====================
mapping = function(genomedir = '/mnt/c/Users/renseb01/Documents/rnaseq/data/reference_genome/chr1_index/',
                   genomefasta = '/mnt/c/Users/renseb01/Documents/rnaseq/data/reference_genome/chr1.fa',
                   threads = 12,
                   R1_trim = sequences()[[3]][1],
                   R2_trim = sequences()[[4]][1],
                   annotation.gtf = '/mnt/c/Users/renseb01/Documents/rnaseq/data/reference_genome/gencode.v43.primary_assembly.annotation_small.gtf',
                   out_prefix = paste0('rnaseq/out/',sequences()[[5]][1],'_')) {
  
  #generate genome if its not there.
  if(length(list.files(gsub("/mnt/c",ifelse(Sys.info()['sysname'] == 'Windows','C:',''),genomedir)))==0) .generate_genome(genomedir = genomedir,genomefasta=genomefasta,annotation.gtf=annotation.gtf)
  
  
  cmd = paste0(ifelse(Sys.info()['sysname'] == 'Windows','wsl.exe ',''),'STAR --genomeDir ',
               genomedir,
               ' --readFilesCommand zcat --readFilesIn ',
               R1_trim,
               ' ',
               R2_trim,
               ' --outSAMtype BAM SortedByCoordinate --outReadsUnmapped Fastx --twopassMode Basic --twopass1readsN -1 --quantMode GeneCounts --sjdbGTFfile ',
               annotation.gtf,
               ' --runThreadN ',
               threads ,
               ' --outFileNamePrefix ', out_prefix)
  
  message(cmd)
  system(cmd)
  return('STAR done')
}


#======================
#  Bam to sam
#=====================
bamtosam = function(out_prefix = 'rnaseq/out/toto_',
                    bam = paste0(out_prefix,'Aligned.sortedByCoord.out.bam'),
                    bam_out = paste0(out_prefix,'trimmed_Aligned_PP_UM.bam')){

  cmd1=paste0(ifelse(Sys.info()['sysname'] == 'Windows','wsl.exe ',''),"samtools view -h -f 0x0002 ",bam," |  grep -P '^@|NH:i:1\t' | samtools view -h -b -  >",bam_out)
  cmd2=paste0(ifelse(Sys.info()['sysname'] == 'Windows','wsl.exe ',''),"samtools view -h -o",gsub(".bam",".sam",bam_out)," ",bam_out)

  system(cmd1)
  system(cmd2)
  
  message(cmd1)
  message(cmd2)
  
  return('samtools bamtosam done')
}

#======================
# Add RG, sort, mark 
# duplicates & create 
# index
#======================
picardtools = function(out_prefix = 'rnaseq/out/toto_',
                       out_prefix_nopath = 'toto_',
                       bam_out = paste0(out_prefix,'trimmed_Aligned_PP_UM.bam'),
                      bam_added = paste0(out_prefix,'trimmed_Aligned_PP_UM_rgAdded.bam'),
                      bam_rmdup = paste0(out_prefix,'trimmed_Aligned_PP_UM_rgAdded_dup.bam'),
                     metrics = paste0(out_prefix,'trimmed_Aligned_PP_UM_rgAdded_metrics.txt')){
  
  cmd1 = paste0(ifelse(Sys.info()['sysname'] == 'Windows','wsl.exe ',''),
                'PicardCommandLine AddOrReplaceReadGroups -I ',
                gsub(".bam",".sam",bam_out),
                ' -O ',
                bam_added,
                ' -SO coordinate -RGID ',out_prefix_nopath,' -RGLB ',out_prefix_nopath,' -RGPL illumina -RGPU 1 -RGSM ',out_prefix_nopath)
  
  
  cmd2 = paste0(ifelse(Sys.info()['sysname'] == 'Windows','wsl.exe ',''),
                'PicardCommandLine MarkDuplicates -I ',
                bam_added,
                ' -O ',
                bam_rmdup,
                ' -CREATE_INDEX TRUE -VALIDATION_STRINGENCY SILENT -M ',
                metrics)
  
  
  message(cmd1)
  message(cmd2)
  
  system(cmd1)
  system(cmd2)
  return('picard tools MarkDuplicates done')
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
                             bam_nounmapped = paste0(out_prefix,'trimmed_Aligned_PP_UM_rgAdded_dup_split_sortP_noUnmapped.bam'),
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
  
  message(sort_noGATK)
  message(view)
  message(index)
  message(resort)
  
  return('Done sort_bam_by_name')
}
  


#======================
# htseq-count step 
#======================
htseq_count = function (
    out_prefix = 'rnaseq/out/toto_',
    annotation.gtf = '/mnt/c/Users/renseb01/Documents/rnaseq/data/reference_genome/gencode.v43.primary_assembly.annotation_small.gtf',
    bam_nounmapped_sort = paste0(out_prefix,'trimmed_Aligned_PP_UM_rgAdded_dup_split_noUnmapped_sortN.bam'),
    htseq_fwd = paste0(out_prefix,'trimmed_Aligned_PP_UM_rgAdded_dup_split_noUnmapped_sortN_htseq_fwd.txt'),
    htseq_reverse = paste0(out_prefix,'trimmed_Aligned_PP_UM_rgAdded_dup_split_noUnmapped_sortN_htseq_reverse.txt'),
    htseq_onstranded = paste0(out_prefix,'trimmed_Aligned_PP_UM_rgAdded_dup_split_noUnmapped_sortN_htseq_unstranded.txt'))
                        {
  
  
  htseq1 = paste0(ifelse(Sys.info()['sysname'] == 'Windows','wsl.exe ',''),'htseq-count --format=bam --mode=intersection-nonempty --stranded=yes --idattr=gene_id ',bam_nounmapped_sort,' ',annotation.gtf, ' >',htseq_fwd)
  htseq2 = paste0(ifelse(Sys.info()['sysname'] == 'Windows','wsl.exe ',''),'htseq-count --format=bam --mode=intersection-nonempty --stranded=reverse --idattr=gene_id ',bam_nounmapped_sort,' ',annotation.gtf, ' >',htseq_reverse)
  htseq3 = paste0(ifelse(Sys.info()['sysname'] == 'Windows','wsl.exe ',''),'htseq-count --format=bam --mode=intersection-nonempty --stranded=no --idattr=gene_id ',bam_nounmapped_sort,' ',annotation.gtf, ' >',htseq_onstranded)

  
  system(htseq1)
  system(htseq2)
  system(htseq3)
  
  message(htseq1)
  message(htseq2)
  message(htseq3)
  
return('Done htseq')
}


#======================
# cleanup 
#======================
cleanup = function (out_prefix = 'rnaseq/out/toto_')
  {
  message('Nothing so far, but will need to empty the garbage eventually')
  }


