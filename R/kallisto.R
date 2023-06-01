#kallisto function
source('R/counts_functions.R')


#=====================
#Pseudo-mapping step
#=====================
index = function(kallindex.dir = 'data/index')
{
    
    #index cmd, index downloaded from here: https://www.gencodegenes.org/human/ (gencode.v43.transcripts.idx)
    cmd = paste0(ifelse(Sys.info()['sysname'] == 'Windows','wsl.exe ',''),'kallisto index --index=', kallindex.dir, 'chr1.idx',' --make-unique ', kallindex.dir, 'chr1.fasta.gz')
    system(cmd)
    
    # message(cmd)
    message(paste0('Done index, Time is: ', Sys.time()))
    
}


kallisto = function(trim.dir = 'out/fastq.trim',
                    kallindex.dir = 'data/index',
                    trim1 = list.files(path = "out/fastq.trim", pattern = "val_1.fq.gz"),
                    trim2 = list.files(path = "out/fastq.trim", pattern = "val_2.fq.gz"),
                    quant.dir = 'out/kallisto',
                    i = 1)
{
    #kallisto cmd
    cmd = paste0(ifelse(Sys.info()['sysname'] == 'Windows','wsl.exe ',''),'kallisto quant --index=', kallindex.dir,'chr1.idx',
                 ' -o ',
                 quant.dir,
                 ' -b 100 ',
                 trim.dir,
                 trim1,
                 ' ',
                 trim.dir,
                 trim2,
                 ' 1>',ifelse(i>1,'>',''), 'out/logs/kallisto.out 2>',ifelse(i>1,'>',''),'out/logs/kallisto.err')
    
    system(cmd)
    
    # message(cmd)
    message(paste0('Done Kallisto, Time is: ',Sys.time()))
    
}


#=====================
#Wrapper: kallisto_rnaseq
#=====================

kallisto_rnaseq = function(kallindex.dir = params$kallindex.dir,
                           trim.dir = params$trim.dir,
                           quant.dir = params$quant.dir,
                           out.dir = params$out.dir)
{
    if(length(list.files(gsub("/mnt/c",ifelse(Sys.info()['sysname'] == 'Windows','C:',''),trim.dir))==0))
        
    {
        #fastq files
        fastq_files = sequences(fq.dir='data/fastq',
                                trim.dir='out/fastq.trim',
                                out.dir = 'out/')
        
        #trimming
        trimming(trim.dir = trim.dir,
                 R1 = fastq_files[[1]][i],
                 R2 = fastq_files[[2]][i],
                 out_seq = fastq_files[[4]][i],
                 cutadapt = '/usr/bin/cutadapt',
                 i=i,
                 out.dir = 'out/')
    }
    
    #generate index if its not there.
    if(length(list.files(gsub("/mnt/c",ifelse(Sys.info()['sysname'] == 'Windows','C:',''),kallindex.dir)))!=2) 
        
    {
        index(kallindex.dir = kallindex.dir)
    }
    
    #name of files
    trim1 = list.files(path = "out/fastq.trim", pattern = "val_1.fq.gz")
    trim2 = list.files(path = "out/fastq.trim", pattern = "val_2.fq.gz")
    name_samples = sapply(strsplit(trim1,"_val_"), `[`, 1)
    
    for(i in seq_along(name_samples)) 
    {
        #outpout directory
        dir.create(file.path("out/kallisto", name_samples[i]), recursive = TRUE)
        
        #kallisto cmd
        kallisto(trim.dir = trim.dir,
                 kallindex.dir = kallindex.dir,
                 trim1 = trim1[i],
                 trim2 = trim2[i],
                 quant.dir = file.path("out/kallisto", name_samples[i]),
                 i = i)
        
        #zip files
        #files_zip = file.path(out_prefix, name_samples[i])
        #zip(zipfile = files_zip, files = files_zip)
        
        #message
        message(paste0('--- Done sample,', name_samples[i] ,' Time is: ',Sys.time(),' ---'))
        
    }
    
}