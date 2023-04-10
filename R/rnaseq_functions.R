


params = list(fq.dir = 'C:/Users/renseb01/Documents/rnaseq/data/fastq',
              trim.dir = 'C:/Users/renseb01/Documents/rnaseq/data/fastq.trim')


###get sequences ready
sequences = function(fq.dir=params$fq.dir,
                     trim.dir=params$trim.dir) {
  R12 = list.files(fq.dir,pattern = 'R[12].fastq.gz')
  sample_names = unique(gsub('_R[12].fastq.gz','',R12))
  
  R1 = paste0(file.path(sub('C:','/mnt/c',fq.dir),sample_names),'_R1.fastq.gz')
  R2 = paste0(file.path(sub('C:','/mnt/c',fq.dir),sample_names),'_R2.fastq.gz')

  
  R1_trim = paste0(file.path(sub('C:','/mnt/c',trim.dir),sample_names),'_R1_val_1.fq.gz')
  R2_trim = paste0(file.path(sub('C:','/mnt/c',trim.dir),sample_names),'_R1_val_1.fq.gz')
  
  return(list(R1,R2,R1_trim,R2_trim))
}


###trimming (note that trim_galore is Unix-based)
trimming = function(trim.dir = params$trim.dir,
                    R1 = sequences()[[1]][1],
                    R2 = sequences()[[2]][1],
                    adaptor1='AGATCGGAAGAGCACACGTCTGAACTCCAGTCA',
                    adaptor2='AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT',
                    cutadapt = '/usr/bin/cutadapt') {
  
  if(file.exists(trim.dir)==F) dir.create(trim.dir)

 
  cmd = paste0(ifelse(Sys.info()['sysname'] == 'Windows','wsl.exe ',''),
               'trim_galore -q 20 -o ',
               sub('C:','/mnt/c',trim.dir),
               ' --phred33 --adapter=',adaptor1,' --adapter2=',adaptor2,' --paired ',
               R1,
               ' ',
               R2,
               ' --path_to_cutadapt ',
               cutadapt, '1>>trim_galore.out 2>>trim_galore.err')
  
  message(cmd)
  system(cmd)

return('Done Trimming')
}



#Prepare genome: note that STAR is Unix-based)
#genome (fasta, gff downloaded from here: https://www.gencodegenes.org/human/)
generate_genome = function(genomefasta = 'C:/Users/renseb01/Documents/rnaseq/data/reference_genome/chr1.fa',
                           genomegtf = 'C:/Users/renseb01/Documents/rnaseq/data/reference_genome/gencode.v43.primary_assembly.annotation_small.gtf',
                           genomedir = 'C:/Users/renseb01/Documents/rnaseq/data/reference_genome/chr1_index'){
  cmd = paste0(ifelse(Sys.info()['sysname'] == 'Windows','wsl.exe ',''),
                'STAR --runThreadN 6 --runMode genomeGenerate --genomeDir ',
               sub('C:','/mnt/c',genomedir),
               ' --genomeFastaFiles ',
               sub('C:','/mnt/c',genomefasta),
               ' --sjdbGTFfile ',
               sub('C:','/mnt/c',genomegtf),
               ' --sjdbOverhang 99 1>>star_aligner.out 2>>star_aligner.err')
  
  message(cmd)
  system(cmd)
  return('Done generate_genome')
}


#Mapping step (note that star is Unix-based)
mapping = function(genomedir = 'C:/Users/renseb01/Documents/rnaseq/data/reference_genome/chr1_index/',
                   threads = 12,
                   R1_trim = sequences()[[1]][3],
                   R2_trim = sequences()[[1]][4],
                   annotation.gtf = 'gencode.v38.annotation.gtf') {
  
  #generate genome if its not there.
  if(length(list.files(genomedir,''))==0) generate_genome(genomedir = genomedir)
  
  
  cmd = paste0(ifelse(Sys.info()['sysname'] == 'Windows','wsl.exe ',''),'STAR --genomeDir ',
               sub('C:','/mnt/c',genomedir),
               '--readFilesCommand zcat --readFilesIn ',
               R1_trim,
               ' ',
               R2_trim,
               '--outSAMtype BAM SortedByCoordinate --outReadsUnmapped Fastx --twopassMode Basic --twopass1readsN -1 --quantMode GeneCounts --sjdbGTFfile ',
               annotation.gtf,
               ' --runThreadN ',
               threads ,
               '--outFileNamePrefix toto')
  
  message(cmd)
  system(cmd)
  return('STAR done')
}


#======================
#  Bam to sam
#=====================
bamtosam = function(){

  cmd1=paste0('samtools view -h -f 0x0002 /Users/jerry/Documents/IUCPQ/rnaseq/src/star/806rcbc2881/806rcbc2881_Aligned.sortedByCoord.out.bam |  grep -P "^@|NH:i:1\t" | samtools view -h -b -  > /Users/jerry/Documents/IUCPQ/rnaseq/src/star/806rcbc2881/806rcbc2881_trimmed_Aligned_PP_UM.bam')
  cmd2=paste0('samtools view -h -o /Users/jerry/Documents/IUCPQ/rnaseq/src/star/806rcbc2881/806rcbc2881_trimmed_Aligned_PP_UM.sam /Users/jerry/Documents/IUCPQ/rnaseq/src/star/806rcbc2881/806rcbc2881_trimmed_Aligned_PP_UM.bam')

  system(cmd1)
  system(cmd2)
  
  return('samtools bamtosam done')
}

#======================
# Add RG, sort, mark 
# duplicates & create 
# index
#======================
java -Xmx20000m -jar 
/lustre03/project/6032391/GROUP/bin/picard-2.20.0/picard/build/libs/picard.jar \
AddOrReplaceReadGroups \
I=/Users/jerry/Documents/IUCPQ/rnaseq/src/star/806rcbc2881/806rcbc2881_trimmed_Aligned_PP_UM.sam O=/Users/jerry/Documents/IUCPQ/rnaseq/src/star/806rcbc2881/806rcbc2881_trimmed_Aligned_PP_UM_rgAdded.bam \
SO=coordinate \
RGID=806rcbc2881 \
RGLB=806rcbc2881 \
RGPL=illumina \
RGPU=1 \
RGSM=806rcbc2881 \
java -Xmx20000m -jar /lustre03/project/6032391/GROUP/bin/picard-2.20.0/picard/build/libs/picard.jar MarkDuplicates I=/Users/jerry/Documents/IUCPQ/rnaseq/src/star/806rcbc2881/806rcbc2881_trimmed_Aligned_PP_UM_rgAdded.bam O=/Users/jerry/Documents/IUCPQ/rnaseq/src/star/806rcbc2881/806rcbc2881_trimmed_Aligned_PP_UM_rgAdded_dup.bam CREATE_INDEX=TRUE VALIDATION_STRINGENCY=SILENT M=/Users/jerry/Documents/IUCPQ/rnaseq/src/star/806rcbc2881/806rcbc2881_trimmed_Aligned_PP_UM_rgAdded_metrics.txt 


#======================
# splitNtrim & reassign
# mapping quality
#======================
#/lustre03/project/6032391/GROUP/bin/gatk-4.1.2.0/gatk SplitNCigarReads -R /lustre03/project/6032391/GROUP/References/GRCh38/Homo_sapiens_assembly38.fasta -I /Users/jerry/Documents/IUCPQ/rnaseq/src/star/806rcbc2881/806rcbc2881_trimmed_Aligned_PP_UM_rgAdded_dup.bam -O /Users/jerry/Documents/IUCPQ/rnaseq/src/star/806rcbc2881/806rcbc2881_trimmed_Aligned_PP_UM_rgAdded_dup_split.bam




#======================
# sort cleaned bam remove unmapped reads   sort bam by name 
#======================
samtools sort -o /Users/jerry/Documents/IUCPQ/rnaseq/src/star/806rcbc2881/806rcbc2881_trimmed_Aligned_PP_UM_rgAdded_dup_split_sortP.bam /Users/jerry/Documents/IUCPQ/rnaseq/src/star/806rcbc2881/806rcbc2881_trimmed_Aligned_PP_UM_rgAdded_dup_split.bam
samtools view -b -F 4 /Users/jerry/Documents/IUCPQ/rnaseq/src/star/806rcbc2881/806rcbc2881_trimmed_Aligned_PP_UM_rgAdded_dup_split_sortP.bam > /Users/jerry/Documents/IUCPQ/rnaseq/src/star/806rcbc2881/806rcbc2881_trimmed_Aligned_PP_UM_rgAdded_dup_split_sortP_noUnmapped.bam
samtools index /Users/jerry/Documents/IUCPQ/rnaseq/src/star/806rcbc2881/806rcbc2881_trimmed_Aligned_PP_UM_rgAdded_dup_split_sortP_noUnmapped.bam
samtools sort -n -o /Users/jerry/Documents/IUCPQ/rnaseq/src/star/806rcbc2881/806rcbc2881_trimmed_Aligned_PP_UM_rgAdded_dup_split_noUnmapped_sortN.bam /Users/jerry/Documents/IUCPQ/rnaseq/src/star/806rcbc2881/806rcbc2881_trimmed_Aligned_PP_UM_rgAdded_dup_split_sortP_noUnmapped.bam





#======================
# htseq-count step 
#======================
cd /Users/jerry/Documents/IUCPQ/rnaseq/src/htseq_count

/lustre03/project/6032391/GROUP/Programs/python_virtualenv/bin/htseq-count --format=bam --mode=intersection-nonempty --stranded=yes --idattr=gene_id /Users/jerry/Documents/IUCPQ/rnaseq/src/star/806rcbc2881/806rcbc2881_trimmed_Aligned_PP_UM_rgAdded_dup_split_noUnmapped_sortN.bam /lustre03/project/6032391/GROUP/References/GRCh38/gencode.v38.annotation.gtf > /Users/jerry/Documents/IUCPQ/rnaseq/src/htseq_count/806rcbc2881_trimmed_Aligned_PP_UM_rgAdded_dup_split_noUnmapped_sortN_htseq_fwd.txt

/lustre03/project/6032391/GROUP/Programs/python_virtualenv/bin/htseq-count --format=bam --mode=intersection-nonempty --stranded=reverse --idattr=gene_id /Users/jerry/Documents/IUCPQ/rnaseq/src/star/806rcbc2881/806rcbc2881_trimmed_Aligned_PP_UM_rgAdded_dup_split_noUnmapped_sortN.bam /lustre03/project/6032391/GROUP/References/GRCh38/gencode.v38.annotation.gtf > /Users/jerry/Documents/IUCPQ/rnaseq/src/htseq_count/806rcbc2881_trimmed_Aligned_PP_UM_rgAdded_dup_split_noUnmapped_sortN_htseq_reverse.txt

/lustre03/project/6032391/GROUP/Programs/python_virtualenv/bin/htseq-count --format=bam --mode=intersection-nonempty --stranded=no --idattr=gene_id /Users/jerry/Documents/IUCPQ/rnaseq/src/star/806rcbc2881/806rcbc2881_trimmed_Aligned_PP_UM_rgAdded_dup_split_noUnmapped_sortN.bam /lustre03/project/6032391/GROUP/References/GRCh38/gencode.v38.annotation.gtf > /Users/jerry/Documents/IUCPQ/rnaseq/src/htseq_count/806rcbc2881_trimmed_Aligned_PP_UM_rgAdded_dup_split_noUnmapped_sortN_htseq_unstranded.txt

