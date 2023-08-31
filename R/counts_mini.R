

#=====================
###get sequences ready
#=====================
sequences = function(fq.dir='data/fastq',
                     trim.dir='out/fastq.trim',
                     mini.dir='out/fastq.mini',
                     out.dir = 'out/') {
    R12 = list.files(gsub("/mnt/c",ifelse(Sys.info()['sysname'] == 'Windows','C:',''),fq.dir),pattern = 'R[12].fastq.gz$')
    sample_names = unique(gsub('_R[12].fastq.gz','',R12))
    
    R1 = paste0(file.path(fq.dir,sample_names),'_R1.fastq.gz')
    R2 = paste0(file.path(fq.dir,sample_names),'_R2.fastq.gz')
    
    R1_trim = paste0(file.path(trim.dir,sample_names),'_R1_val_1.fq.gz')
    R2_trim = paste0(file.path(trim.dir,sample_names),'_R2_val_2.fq.gz')
    
    R1_mini = paste0(file.path(trim.dir,sample_names),'_R1_mini.fq')
    R2_mini = paste0(file.path(trim.dir,sample_names),'_R2_mini.fq')
    
    if(length(list.files(trim.dir)>0)) {
    message('It appears that your trimdir is not empty, I will list what is in there in output[[2]],[[3]] & [[4]]')
    R1_trim = file.path(trim.dir,list.files(trim.dir,pattern = '_R1_val_1.fq.gz$'))
    R2_trim = file.path(trim.dir,list.files(trim.dir,pattern = '_R2_val_2.fq.gz$'))
    sample_names = unique(gsub('_R1_val_1.fq.gz','',list.files(trim.dir,pattern = '_R1_val_1.fq.gz$')))
    }
    
    out.dir.length = length(list.files(gsub("/mnt/c",ifelse(Sys.info()['sysname'] == 'Windows','C:',''),out.dir))) 
    
    if(length(R12)==0) warning('It appears like your directory contains no .fastq.gz file',immediate. = TRUE)
    #if(out.dir.length>1) warning(paste0(out.dir,' already contains ',out.dir.length,' files. I will add files to this directory'),immediate. = TRUE)
    message(paste0('Done listing sequences, ',length(sample_names),' files to process, Time is: ',Sys.time()))
    
    return(list(R1,R2,R1_trim,R2_trim,sample_names,R1_mini,R2_mini))
}


#=====================
###zcat, random (1%)
#=====================
minisample = function(r1 = sequences()[[1]][1],
                      r2 = sequences()[[2]][1],
                      r1_mini = sequences()[[6]][1],
                      r2_mini = sequences()[[7]][1],
                      head = 40000000){

  print('make sure you have this installed (https://github.com/rieseberglab/fastq-examples) in a place your system can find it')
  
  
    cmd1 = paste0('zcat ', r1,' | head -40000000 >10millions_r1.fq') #10million reads (roughly 10%)
    cmd2 = paste0('zcat ', r2,' | head -40000000 >10millions_r2.fq') 
    cmd3 = paste0('generateFastqSample.py 10 10millions_r1.fq 10millions_r2.fq 3>',r1_mini, ' 4>',r2_mini) # 10 % of 10million (roughly 1% of the data)

    system(cmd1)
    system(cmd2)
    system(cmd3)
    
    message(paste0('Done mini sample (4 million PE reads), Time is: ',Sys.time()))
    
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
                 ' --limitBAMsortRAM 100000000000 --readFilesIn ',
                 R1_trim,
                 ' ',
                 R2_trim,
                 ' --outSAMtype BAM SortedByCoordinate --sjdbOverhang 99 --outFilterMultimapNmax 100 --outReadsUnmapped None --quantMode GeneCounts --sjdbGTFfile ',
                 annotation.gtf,
                 ' --runThreadN ',
                 threads ,
                 ' --outFileNamePrefix ', out_prefix,' 1>',ifelse(i>1,'>',''),file.path(out.dir,'logs'),'/star_aligner.out 2>',ifelse(i>1,'>',''),file.path(out.dir,'logs'),'/star_aligner.err')
    
    if(file.exists(paste0(out_prefix,'Aligned.sortedByCoord.out.bam'))==F) {system(cmd);message(paste0('Done STAR, Time is: ',Sys.time()))} else {message(paste0('Skipping STAR (already done), Time is: ',Sys.time()))}
    
    return('')
}



#======================
# cleanup 
#======================
cleanup = function (out.dir = 'out/alignments')
{

    rm_cmd1 = paste0('rm ',out.dir,'/alignments/*Aligned.sortedByCoord.out.bam*')
    system(rm_cmd1)

    message(paste0('Done cleanup of temporary (bam) because they are really big, Time is: ',Sys.time()))
    return('')
}

#======================
# merge genetabs from STAR 
#======================
mergeSTARGeneTabs = function (out.dir = out.dir,
                samples = fastq_files[[5]][files])
{

    for(i in seq_along(samples)){
      expression_temp = read.delim(paste0(out.dir,'/alignments/',samples[i],'_ReadsPerGene.out.tab'), header = F)
      if(i==1) {
        expression_df = as.data.frame(matrix(0,nrow =nrow(expression_temp), ncol = length(fastq[[5]])+1));
        colnames(expression_df)[1] = 'gene_names';
        colnames(expression_df)[-1] = fastq[[5]];
        expression_df[,1]  =expression_temp[,1]
        }
      expression_df[,i+1] = expression_temp[,4]
        
    }
  
    
    write.csv(expression_df,file.path(out.dir,'expression_df.csv'),row.names =F)
    
    message(paste0('Done meerging all STAR gene tabs, Time is: ',Sys.time()))
    return('')
}



#======================
# wrapper : counts_rnaseq
#======================
counts_test_rnaseq = function(fq.dir = params$fq.dir,
                         trim.dir= params$trim.dir,
                         genomedir = params$genomedir,
                         annotation.gtf = params$annotation.gtf,
                         genomefasta = params$genomefasta,
                         out.dir = params$out.dir,
                         cutadapt = params$cutadapt,
                         threads = params$threads,
                         nbfiles = params$nbfiles,
                         head = params$head){
    
    dir.create(out.dir,showWarnings =F)
    dir.create(file.path(out.dir,'logs'),showWarnings =F)
    dir.create(file.path(out.dir,'alignments'),showWarnings =F)
    
    #fastq files
    fastq_files = sequences(fq.dir = fq.dir,
                            trim.dir = trim.dir,
                            out.dir = out.dir)
    
    
    #files to process depends on nbfiles parameter
    files = seq_along(fastq_files[[5]])
    nbfiles =  strsplit(nbfiles,',')[[1]]
    if(nbfiles[1] == 'all') files = files
    if(length(nbfiles) == 1 & nbfiles[1] != 'all') {files = files[1:min(length(files),as.numeric(nbfiles))]}
    if(length(nbfiles) == 2) {files = files[as.numeric(nbfiles[1]):as.numeric(nbfiles[2])]}
    
    
    for(i in files) {
        
        #out_prefix
        out_prefix = file.path(out.dir,'alignments',paste0(fastq_files[[5]][i],'_'))

        #mini sample (1%)
        minisample(r1 = sequences()[[1]][i],
                      r2 = sequences()[[2]][i],
                      r1_mini = sequences()[[6]][i],
                      r2_mini = sequences()[[7]][i],
                      head = head)
        
        #mapping
        mapping(genomedir = genomedir,
                genomefasta = genomefasta,
                threads = threads,
                R1_trim = fastq_files[[6]][i],
                R2_trim = fastq_files[[7]][i],
                annotation.gtf = annotation.gtf,
                out_prefix = out_prefix,
                i = i,
                out.dir = out.dir)
              
        #message
        message(paste0('--- Done sample, ',fastq_files[[5]][i],', Time is: ',Sys.time(),' ---'))
    }
    
    #merge
    mergeSTARGeneTabs(out.dir = out.dir,
                samples = fastq_files[[5]][files])
  
    
    #final cleanup
    cleanup(out.dir = out.dir)
    
}
