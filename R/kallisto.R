#kallisto function
source('R/counts_functions.R')


#=====================
#Pseudo-mapping step
#=====================
kallisto = function(R1_trim = sequences()[[3]][1],
                    R2_trim = sequences()[[4]][1],
                    in.dir = file.path('data/index/'),
                    out.dir = file.path('data/kallisto/'),
                    out_prefix = paste0('rnaseq/out/',sequences()[[5]][1],'_')) 
    {
    #if sequences not trimmed
    if(length(list.files(gsub("/mnt/c",ifelse(Sys.info()['sysname'] == 'Windows','C:',''),file.path('out/fastq.trim'))))==0)
        {
        #=====================
        ###get sequences ready
        #=====================
        sequences(fq.dir='data/fastq',
                  trim.dir='out/fastq.trim',
                  out.dir = 'out/')
        
        #=====================
        ###trimming
        #=====================
        trimming(trim.dir = 'out/fastq.trim',
                 R1 = sequences()[[1]][1],
                 R2 = sequences()[[2]][1],
                 out_seq = sequences()[[4]][1],
                 adaptor1='AGATCGGAAGAGCACACGTCTGAACTCCAGTCA',
                 adaptor2='AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT',
                 cutadapt = '/usr/bin/cutadapt',
                 i=1,
                 out.dir=out.dir)
        }
    
    #index cmd, index downloaded from here: https://www.gencodegenes.org/human/ (gencode.v43.transcripts.idx)
    index = paste0(ifelse(Sys.info()['sysname'] == 'Windows','wsl.exe ',''),'kallisto index -i ', in.dir, 'chr1.idx ', in.dir, 'chr1.fa.gz')
    system(index)
    
    #kallisto cmd
    cmd = paste0(ifelse(Sys.info()['sysname'] == 'Windows','wsl.exe ',''),'kallisto quant -i ', in.dir, 'chr1.idx -o ',
                 out.dir,
                 ' -b 100 ', 
                 in.dir ,
                 R1_trim,
                 ' ',
                 in.dir,
                 R2_trim,
                 ' ',
                 out_prefix,' 1>',ifelse(i>1,'>',''),file.path(out.dir,'logs'),'/kallisto.out 2>',ifelse(i>1,'>',''),file.path(out.dir,'logs'),'/kallisto.err')
    
    system(cmd)
    
    # message(cmd)
    message(paste0('Done Kallisto, Time is: ',Sys.time()))
}


#=====================
#Load kallisto data
#=====================
load_data = function(in.dir = "/data/",
                     sample_id = dir(file.path("/data/kallisto")),
                     kal_dirs = file.path("/data/kallisto", sample_id))
    
{
    #Grep info from run file csv & covert to txt
    #data_file = readxl::read_xlsx("/data/*.xlsx", sheet = 1, col_names = TRUE, col_types = NULL, na = "", skip = 0)
    #file.dir = file.path("/data/")
    #data_frame = write.csv(data_file, file.path(file.dir, "test.csv"))
    #df = read.csv("/data/test.csv")
    #col = df[, c("Name...2", "Sample...5", "Type")] > "/data/sampleTable.txt.txt
    
    #Text file of run info
    sampleTable = read.delim(file.path(in.dir, "sampleTable.txt"))
    rownames(sampleTable) = sampleTable$Name
    
    #load TxDB from gencode (Release 43 (GRCh38.p13)
    txdb_v31_hg19_chr1 = GenomicState::gencode_txdb("43", "hg19", chrs = "chr1") #chrs = paste0("chr", c(seq_len(22), "X", "Y", "M"))
    k = keys(txdb_v31_hg19_chr1, keytype = "TXNAME")
    tx2gene = select(txdb_v31_hg19_chr1, k, "GENEID", "TXNAME")
    tx2gene$TXNAME = sapply(strsplit(tx2gene$TXNAME,".",fixed = T), `[`, 1)
    
    #Load kallisto data
    files = file.path(paste0(in.dir, "kallisto"), list.files(paste0(in.dir, "kallisto")), "abundance.h5")
    names(files) = list.files(paste0(in.dir, "kallisto"))
    
    #import abundance.h5 files
    txi = tximport::tximport(files, type = "kallisto", tx2gene = tx2gene, 
                             txIn = TRUE, txOut = FALSE, countsFromAbundance = "lengthScaledTPM",  ignoreTxVersion = TRUE, ignoreAfterBar=TRUE)
    
    # Sleuth
    #SampleTable (info of samples associated with kallisto quant)
    s2c = read.table(file.path(in.dir, "sampleTable.txt"), header = TRUE, stringsAsFactors=FALSE)
    s2c = dplyr::select(s2c, sample = Name, Type)
    s2c = dplyr::mutate(s2c, path = kal_dirs)
    
    #Initialize sleuth object 
    so = sleuth::sleuth_prep(s2c, read_bootstrap_tpm = TRUE)
    
    #Including gene names for transcript-level analysis
    mart = biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL", 
                            dataset = "hsapiens_gene_ensembl",
                            host = 'https:///www.ensembl.org')
    
    t2g = biomaRt::getBM(attributes = c("ensembl_transcript_id","ensembl_transcript_id_version", "ensembl_gene_id", 
                                        "ensembl_gene_id_version","external_gene_name","description",
                                        "chromosome_name","start_position",
                                        "end_position","strand",
                                        "entrezgene_id"), mart = mart)
    
    t2g = dplyr::rename(t2g, target_id = ensembl_transcript_id, 
                        ens_gene = ensembl_gene_id, ext_gene = external_gene_name)
    
    t2g = dplyr::select(t2g, c('target_id', 'ens_gene', 'ext_gene'))
    
    #Adding gene names to sleuth table
    so = sleuth::sleuth_prep(s2c, target_mapping = t2g)
    
    #fitting measurement error models
    so = sleuth::sleuth_fit(so, ~1, 'reduced')
    so = sleuth::sleuth_lrt(so, 'reduced', 'reduced')
    
    #PCA plot
    sleuth::plot_pca(so, color_by = 'Type')
    
    #Group density plot for count distribution/sample type
    
    sleuth::plot_group_density(so, use_filtered = TRUE, units = "tpm", 
                               trans = "log", grouping = setdiff(colnames(so$sample_to_covariates),
                                                                 "sample"), offset = 1)
    #View results with sleuth_live(so)
    
    
    # message(cmd)
    message(paste0('Done data_load, Time is: ',Sys.time()))
}


#=====================
#Wrapper: kallisto_rnaseq
#=====================

kallisto_rnaseq = function(fq.dir = params$fq.dir,
                           trim.dir= params$trim.dir,
                           sample_id = params$sample_id,
                           kal_dirs = params$kal_dirs,
                           in.dir = params$in.dir,
                           out.dir = params$out.dir,
                           cutadapt = params$cutadapt,
                           nbfiles = params$nbfiles)
{    
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
    
    for(i in files) 
        {
        
        #out_prefix
        out_prefix = file.path(out.dir,'alignments',paste0(fastq_files[[5]][i],'_'))
        
        
        #kallisto
        kallisto(R1_trim = sequences()[[3]][1],
                 R2_trim = sequences()[[4]][1],
                 in.dir = file.path('data/index/'),
                 out.dir = file.path('data/kallisto/'),
                 out_prefix = paste0('rnaseq/out/',sequences()[[5]][1],'_'))
        
        #load data
        load_data(in.dir = "/data/",
                  sample_id = dir(file.path("/data/kallisto")),
                  kal_dirs = file.path("/data/kallisto", sample_id))
        
        #message
        message(paste0('--- Done sample, ',fastq_files[[5]][i],', Time is: ',Sys.time(),' ---'))
        }
    
    
    #final cleanup
    cleanup(out_prefix = out_prefix)
    
}
