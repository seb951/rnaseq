#kallisto function
source('R/counts_functions.R')

#=====================
#Index command
#=====================
index <- function(idx.dir = 'data/index/') {
    
    if (!file.exists(file.path(paste0(getwd(), idx.dir)))) dir.create(file.path(paste0(getwd(), idx.dir)))
    
    # Index command, index downloaded from here: https://www.gencodegenes.org/human/ (gencode.v43.transcripts.idx)
    cmd <- paste0('kallisto index -i ', getwd(),idx.dir, 'chr1.idx --make-unique ', getwd(), idx.dir, 'chr1.fasta.gz')
    system(cmd)
    
    # Print message
    message(paste0('Done indexing, Time is: ', Sys.time()))
}


#=====================
#Kallisto command
#=====================
kallisto <- function(trim.dir = 'out/fastq.trim/',
                     idx.dir = 'data/index/',
                     trim1 = list.files(path = "out/fastq.trim", pattern = "val_1.fq.gz"),
                     trim2 = list.files(path = "out/fastq.trim", pattern = "val_2.fq.gz"),
                     quant.dir = 'out/kallisto/',
                     i = 1) {
    
    # Kallisto command
    cmd <- paste0('kallisto quant -i ', file.path(idx.dir, 'chr1.idx'),
                  ' -o ',
                  file.path(quant.dir),
                  ' -b 100 ',
                  file.path(trim.dir, trim1),
                  ' ',
                  file.path(trim.dir, trim2),
                  ' > ', if (i > 1) file.path('out/logs', 'kallisto.out') else file.path('out/logs', 'kallisto.out'),
                  ' 2> ', if (i > 1) file.path('out/logs', 'kallisto.err') else file.path('out/logs', 'kallisto.err'))
    
    system(cmd)
    
    # Print message
    message(paste0('Done Kallisto, Time is: ', Sys.time()))
}



#=====================
#Wrapper: kallisto_rnaseq
#=====================

kallisto_rnaseq = function(idx.dir = params$idx.dir,
                           trim.dir = params$trim.dir,
                           quant.dir = params$quant.dir,
                           out.dir = params$out.dir)
{
    #create out directory
    if (!file.exists(file.path(paste0(getwd(), quant.dir)))) dir.create(file.path(paste0(getwd(), quant.dir)))
    
    #index
    if(length(list.files(gsub("/mnt/c",ifelse(Sys.info()['sysname'] == 'Windows','C:',''),idx.dir)))!=2) index(idx.dir = idx.dir)
    
    
    #name of files
    trim1 = list.files(path = "out/fastq.trim", pattern = "val_1.fq.gz")
    trim2 = list.files(path = "out/fastq.trim", pattern = "val_2.fq.gz")
    name_samples = sapply(strsplit(trim1,"_val_"), `[`, 1)
    
    
    for(i in seq_along(name_samples)) 
    {
        
        #outpout directory
        dir.create(file.path(paste0(getwd(),"/out/kallisto/", name_samples[i])), recursive = TRUE)
        out_prefix = file.path('out/kallisto', name_samples[i])
        
        #kallisto cmd
        kallisto(trim.dir = trim.dir,
                 idx.dir = idx.dir,
                 trim1 = trim1[i],
                 trim2 = trim2[i],
                 quant.dir = out_prefix,
                 i = i)
        
        
        #zip files
        #files_zip = file.path(out_prefix, name_samples[i])
        #zip(zipfile = files_zip, files = files_zip)
        
        #message
        message(paste0('--- Done sample,', name_samples[i] ,' Time is: ',Sys.time(),' ---'))
        
    }
    
    #}
}
