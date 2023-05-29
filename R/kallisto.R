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
    
    files = list.files(path = "out/fastq.trim/", pattern = ".fq.gz")
    trim1 = list.files(path = "out/fastq.trim/", pattern = "val_1.fq.gz")
    trim2 = list.files(path = "out/fastq.trim/", pattern = "val_2.fq.gz")
    
    name_samples = unique(substr(files, 1, 20))
    
    
    #output directory
    path = 'out/kallisto/'
    out_prefix = ifelse(!dir.exists(file.path(path, paste0(name_samples[1]))), dir.create(file.path(path, paste0(name_samples[1]))), FALSE)
    #out_prefix = dir.create(path, paste0(name_samples[1]))
    
    for(i in files)
    {
        #kallisto
        kallisto(fastq.dir = file.path('out/fastq.trim/'),
                 in.dir = file.path('data/index/'),
                 trim1 = list.files(path = "out/fastq.trim/", pattern = "val_1.fq.gz"),
                 trim2 = list.files(path = "out/fastq.trim/", pattern = "val_2.fq.gz"),
                 out.dir = out_prefix,
                 i = 1)
        
        #message
        message(paste0('--- Done sample, Time is: ',Sys.time(),' ---'))
    }
    
}

