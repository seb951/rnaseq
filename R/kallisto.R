#kallisto function

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
                    out.dir= 'out/') {
    
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
#Pseudo-mapping step
#=====================
kallisto = function(R1_trim = sequences()[[3]][1],
                   R2_trim = sequences()[[4]][1],
                   in.dir = file.path('data/index'),
                   out.dir = file.path('data/kallisto'),
                   out_prefix = paste0('rnaseq/out/',sequences()[[5]][1],'_')
                   ) {
    
    #index cmd
    # index downloaded from here: https://www.gencodegenes.org/human/
    #Crééer fichier contenant l'index
    index = paste0(ifelse(Sys.info()['sysname'] == 'Windows','wsl.exe ',''),'kallisto index -i', in.dir, 'gencode.v43.transcripts.idx', in.dir, 'gencode.v43.transcripts.fa.gz')
    system(index)
    
    #kallisto cmd
    cmd = paste0(ifelse(Sys.info()['sysname'] == 'Windows','wsl.exe ',''),'kallisto quant -i', in.dir, 'gencode.v43.transcripts.idx -o',
                 out.dir,
                 '-b 100', 
                 in.dir,
                 R1_trim,
                 ' ',
                 in.dir,
                 R2_trim, 
                 out_prefix,' 1>',ifelse(i>1,'>',''),file.path(out.dir,'logs'),'/kallisto.out 2>',ifelse(i>1,'>',''),file.path(out.dir,'logs'),'/kallisto.err')
    
    system(cmd)
    
    #view results abundance.h5 with sleuth 
    
    # message(cmd)
    
    message(paste0('Done Kallisto, Time is: ',Sys.time()))
    
    return('')
}

#Observe results in R