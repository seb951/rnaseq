#kallisto function
source('R/counts_functions.R')

#=====================
#Pseudo-mapping step
#=====================
index = function(idx.dir = '/data/index/')
{
    #working path
    idx_path = sub('C:','/mnt/c', getwd())
    
    #index cmd, index downloaded from here: https://www.gencodegenes.org/human/ (gencode.v43.transcripts.idx)
    cmd <- paste0('kallisto index -i ', idx_path, idx.dir, 'chr1.idx --make-unique ', idx_path, idx.dir, 'chr1.fasta.gz')
    system(cmd)
    
    # message(cmd)
    message(paste0('Done index, Time is: ', Sys.time()))
    
}


kallisto = function(trim.dir = '/out/fastq.trim/',
                    idx.dir = '/data/index/',
                    trim1 = list.files(path = "out/fastq.trim", pattern = "val_1.fq.gz"),
                    trim2 = list.files(path = "out/fastq.trim", pattern = "val_2.fq.gz"),
                    quant.dir = '/out/kallisto',
                    i = 1)
{
    #working path
    kal_path = sub('C:','/mnt/c', getwd())
    
    #kallisto cmd
    cmd <- paste0('kallisto quant -i ', kal_path, idx.dir,'chr1.idx',
                 ' -o ',
                 kal_path,
                 quant.dir,
                 ' -b 100 ',
                 kal_path,
                 trim.dir,
                 trim1,
                 ' ',
                 kal_path,
                 trim.dir,
                 trim2,
                 ' 1 > ', if(i > 1) file.path(kal_path, 'out/logs', 'kallisto.out') else file.path(kal_path,'out/logs', 'kallisto.out'),
                 ' 2 > ', if (i > 1) file.path(kal_path, 'out/logs', 'kallisto.err') else file.path(kal_path, 'out/logs', 'kallisto.err'))
    
    system(cmd)
    
    # message(cmd)
    message(paste0('Done Kallisto, Time is: ',Sys.time()))
    
}


#=====================
#Wrapper: kallisto_rnaseq
#=====================

kallisto_rnaseq = function(idx.dir = params$idx.dir,
                           trim.dir = params$trim.dir,
                           quant.dir = params$quant.dir,
                           out.dir = params$out.dir)
{
    #working path
#    trim_path = paste0(sub('C:','/mnt/c', getwd()), trim.dir)
    
#    if(length(list.files(gsub(trim_path)))==0)
        
#    {
        #fastq files
#        fastq_files = sequences(fq.dir = paste0(getwd(), '/data/fastq'),
#                                trim.dir =paste0(getwd(),'/out/fastq.trim'),
#                                out.dir = paste0(getwd(),'/out'))
        
        #trimming
#        trimming(trim.dir = trim_path,
#                 R1 = fastq_files[[1]][i],
#                 R2 = fastq_files[[2]][i],
#                 out_seq = fastq_files[[4]][i],
#                 cutadapt = '/usr/bin/cutadapt',
#                 i=i,
#                 out.dir = paste0(trim_path,'/out'))
#    }
    
    #working path
    ref_path = paste0(sub('C:','/mnt/c', getwd()), idx.dir)
    
    #generate index if its not there.
    if(length(list.files(ref_path)!=2)) 
        
    {
        index(idx.dir = '/data/index/')
    }
    
    #name of files
    trim1 = list.files(path = "out/fastq.trim", pattern = "val_1.fq.gz")
    trim2 = list.files(path = "out/fastq.trim", pattern = "val_2.fq.gz")
    name_samples = sapply(strsplit(trim1,"_val_"), `[`, 1)
    
    for(i in seq_along(name_samples)) 
    {
        #outpout directory
        dir.create(file.path(getwd(),"out/kallisto", name_samples[i]), recursive = TRUE)
        out_prefix = paste0('/out/kallisto/', name_samples[i])
        
        #kallisto cmd
        kallisto(trim.dir = '/out/fastq.trim/',
                 idx.dir = '/data/index/',
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
    
}