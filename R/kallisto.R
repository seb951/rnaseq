#kallisto function
source('R/counts_functions.R')


#=====================
#Pseudo-mapping step
#=====================
kallisto = function(fastq.dir = file.path('out/fastq.trim/'),
                    out_prefix = 'out/kallisto/toto_',
                    R1_trim = sequences()[[3]][1],
                    R2_trim = sequences()[[4]][1],
                    in.dir = file.path('data/index/'),
                    out.dir = file.path('out/kallisto'),
                    i = 1)
    {
    
    #index cmd, index downloaded from here: https://www.gencodegenes.org/human/ (gencode.v43.transcripts.idx)
    index = paste0(ifelse(Sys.info()['sysname'] == 'Windows','wsl.exe ',''),'kallisto index -i ', in.dir, 'chr1.idx ', ' --make-unique ', in.dir, 'chr1.fa.gz')
    system(index)
    
    # message(cmd)
    message(paste0('Done index, Time is: ',Sys.time()))
    
    #files
    files = list.files(path = "out/fastq.trim/", pattern = ".fq.gz")
    
    for(i in files)
        {
        out_prefix_nopath = paste0('out/kallisto/',sequences()[[5]][1])
        
        #kallisto cmd
        cmd = paste0(ifelse(Sys.info()['sysname'] == 'Windows','wsl.exe ',''),'kallisto quant -i ', in.dir, 'chr1.idx -o ',
                     out_prefix_nopath,
                     ' -b 100 ', 
                     fastq.dir,
                     R1_trim,
                     ' ',
                     fastq.dir,
                     R2_trim,
                     ' ', '1>', ifelse(i>1, '>', ''), file.path('out/logs'), '/kallisto.out 2>', ifelse(i>1, '>', ''), file.path('out/logs'), '/kallisto.err')
        
        system(cmd)
        
        # message(cmd)
        message(paste0('Done Kallisto, Time is: ',Sys.time()))
        }
    }




#=====================
#Wrapper: kallisto_rnaseq
#=====================

kallisto_rnaseq = function(fastq.dir = params$fastq.dir,
                           out_prefix = params$out_prefix,
                           in.dir = params$in.dir,
                           out.dir = params$out.dir)
{  
        
        #kallisto
        kallisto(fastq.dir = file.path('out/fastq.trim/'),
                 out_prefix = 'out/kallisto/toto_',
                 R1_trim = sequences()[[3]][1],
                 R2_trim = sequences()[[4]][1],
                 in.dir = file.path('data/index/'),
                 out.dir = file.path('out/kallisto'),
                 i = 1)
        
        #message
        message(paste0('--- Done sample, Time is: ',Sys.time(),' ---'))
    
}
