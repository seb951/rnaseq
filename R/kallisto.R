#kallisto function
source('R/counts_functions.R')


#=====================
#Pseudo-mapping step
#=====================
index = function(in.dir = file.path('data/index/'))
{
    #index cmd, index downloaded from here: https://www.gencodegenes.org/human/ (gencode.v43.transcripts.idx)
    index = paste0(ifelse(Sys.info()['sysname'] == 'Windows','wsl.exe ',''),'kallisto index -i ', in.dir, 'chr1.idx', ' --make-unique ', in.dir, 'chr1.fa.gz')
    system(index)
    
    # message(cmd)
    message(paste0('Done index, Time is: ',Sys.time()))
    
}


kallisto = function(fastq.dir = file.path('out/fastq.trim/'),
                    in.dir = file.path('data/index/'),
                    trim1 = list.files(path = "out/fastq.trim/", pattern = "val_1.fq.gz"),
                    trim2 = list.files(path = "out/fastq.trim/", pattern = "val_2.fq.gz"),
                    out.dir = file.path('out/kallisto/'),
                    i = 1)
    {
    #kallisto cmd
    cmd = paste0(ifelse(Sys.info()['sysname'] == 'Windows','wsl.exe ',''),'kallisto quant -i ', in.dir, 'chr1.idx ',
                     '-o ',
                     out.dir,
                     ' -b 100 ',
                     fastq.dir,
                     trim1,
                     ' ',
                     fastq.dir,
                     trim2,
                     ' ', '1>', ifelse(i>1, '>', ''), file.path('out/logs'), '/kallisto.out 2>', ifelse(i>1, '>', ''), file.path('out/logs'), '/kallisto.err')
        
        system(cmd)
        
        # message(cmd)
        message(paste0('Done Kallisto, Time is: ',Sys.time()))

    }




#=====================
#Wrapper: kallisto_rnaseq
#=====================

kallisto_rnaseq = function(fastq.dir = params$fastq.dir,
                           in.dir = params$in.dir,
                           out.dir = params$out.dir,
                           out_prefix = params$out_prefix)
{
    index(in.dir = file.path('data/index/'))
    
    trim1 = list.files(path = "out/fastq.trim/", pattern = "val_1.fq.gz")
    trim2 = list.files(path = "out/fastq.trim/", pattern = "val_2.fq.gz")
    
    name_samples = sapply(strsplit(trim1,"_val_"), `[`, 1)
    
    for(i in seq_along(name_samples)) 
    {
        #create outpour
        out_prefix = "out/kallisto"
        dir.create(file.path(out_prefix, name_samples[i]), recursive = TRUE)
        
        #kallisto
        kallisto(fastq.dir = file.path('out/fastq.trim/'),
                 in.dir = file.path('data/index/'),
                 trim1 = trim1[i],
                 trim2 = trim2[i],
                 out.dir = file.path(out_prefix, name_samples[i]),
                 i = 1)
        #zip files
        #files_zip = file.path(out_prefix, name_samples[i])
        #zip(zipfile = files_zip, files = files_zip)
        
        #message
        message(paste0('--- Done sample, Time is: ',Sys.time(),' ---'))
    }
    
}