# kallisto function
source('R/counts_functions.R')

#libraries needed: GenomicFeatures, AnnotationDbi, tximport & DESeq2

#=====================
# Index command
#=====================
index <- function(ref.transcriptome = 'data/reference_transcriptome/chr1.fasta.gz'){
    # Reference transcriptome: full transcriptome can be found here: https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_43/gencode.v43.transcripts.fa.gz

    # Kallisto index command 
    cmd <- paste0(ifelse(Sys.info()['sysname'] == 'Windows','wsl.exe ',''),
                  'kallisto index -i ',
                  gsub('fasta.gz','idx',ref.transcriptome),
                  ' --make-unique ',
                  ref.transcriptome)
    system(cmd)
    
    # Print message
    message(paste0('Done Kallisto indexing, Time is: ', Sys.time()))
    return('')
}


#=====================
# kallisto command
#=====================
kallisto <- function(trim1 = 'out/fastq.trim/RNA_0006_5598_Tumeur_R1_val_1.fq.gz',
                     trim2 = 'out/fastq.trim/RNA_0006_5598_Tumeur_R2_val_2.fq.gz',
                     quant.dir = 'out/kallisto',
                     ref.transcriptome = 'data/reference_transcriptome/chr1.fasta.gz',
                     threads = 12,
                     i = 1) {
    # Create output directory
    if (!file.exists(quant.dir)) {
        dir.create(quant.dir, recursive = TRUE)
    }
    
    
    # Kallisto quant command
    cmd <- paste0(ifelse(Sys.info()['sysname'] == 'Windows','wsl.exe ',''),
                  'kallisto quant -i ',
                  gsub('fasta.gz','idx',ref.transcriptome),
                  ' -o ',
                  file.path(quant.dir),
                  ' -b 100 ',
                  trim1,
                  ' ',
                  trim2,
                  ' --threads ', threads,
                  ' 1>',ifelse(i>1,'>',''),
                  file.path('out/logs', 'kallisto.out'),
                  ' 2>',ifelse(i>1,'>',''),
                  file.path('out/logs', 'kallisto.err'))
    
    system(cmd)
    
    # Print message
    message(paste0('Done Kallisto count, Time is: ', Sys.time()))
    return('')
}


#==========================================
# TxDb object step
#==========================================
ref_txb <- function(txdb.dir = "data/reference_transcriptome/gencode.v43.annotation.gtf",
                    txdb.file = "data/reference_transcriptome/gencode.v43.annotation.sqlite"){
    
    # Create TxDb database and save to SQLite for later use
    gtf = file.path(txdb.dir)
    txdb.filename = file.path("data/reference_transcriptome/gencode.v43.annotation.sqlite")
    txdb = GenomicFeatures::makeTxDbFromGFF(gtf)
    AnnotationDbi::saveDb(txdb, txdb.filename)
}


#==========================================
# tximport counts step
#==========================================
# gtf downloaded from here: https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_43/gencode.v43.annotation.gtf.gz
tximp_counts = function(txdb.dir = "data/reference_transcriptome/gencode.v43.annotation.gtf",
                        quant.dir = "out/kallisto",
                        out.dir = "out/gene_counts",
                        txdb.file = "data/reference_transcriptome/gencode.v43.annotation.sqlite"){
    
    # Load TxDB & create tx2gene
    txdb = AnnotationDbi::loadDb(txdb.file)
    k = AnnotationDbi::keys(txdb, keytype = "TXNAME")
    tx2gene = AnnotationDbi::select(txdb, k, "GENEID", "TXNAME")
    #tx2gene$TXNAME = sapply(strsplit(tx2gene$TXNAME,".",fixed = T), `[`, 1)
    
    # Import kallisto data
    files = file.path(paste0(quant.dir),list.files(quant.dir, pattern = "_R1"), "abundance.tsv")
    names(files) = list.files(quant.dir, pattern = "_R1")
    txi_tsv = tximport::tximport(files, type = "kallisto", tx2gene = tx2gene, ignoreAfterBar = TRUE)
    count_df = round(txi_tsv$counts)
    write.table(count_df,file.path(out.dir,'kallisto_counts.tsv'),row.names =T, quote = F, sep = '\t')
    df_count = read.csv(file.path(out.dir,'kallisto_counts.tsv'),sep = '\t',header = F)
    colnames(df_count) = c('GENE_ID', names(files))
    write.table(df_count,file.path(out.dir,'kallisto_counts.tsv'),row.names =F, quote = F,sep = '\t')
    df_c = read.csv(file.path(out.dir,'kallisto_counts.tsv'),sep = '\t',header = T)
    df_c[-1,]
    write.table(df_c[-1,],file.path(out.dir,'kallisto_counts.tsv'),row.names =F, quote = F,sep = '\t')
    

    # Print message
    message(paste0('Done tximport, Time is: ', Sys.time()))
}


#==========================================
# Wrapper: kallisto_rnaseq
#==========================================
kallisto_rnaseq <- function(idx.dir ='data/index',
                            ref.transcriptome = 'data/reference_transcriptome/chr1.fasta.gz',
                            trim.dir = 'out/fastq.trim',
                            quant.dir = 'out/kallisto',
                            threads = 12,
                            txdb.dir = "data/reference_transcriptome/gencode.v43.annotation.gtf",
                            txdb.file = "data/reference_transcriptome/gencode.v43.annotation.sqlite"){
 
    # Index, if not present yet, but reference is there, otherwise, kaboom!...
    if(file.exists(ref.transcriptome) & !file.exists(gsub('fasta.gz','idx',ref.transcriptome))) {
        index(ref.transcriptome = ref.transcriptome)
    }
    
    # Name of sequencing files
    sequencing_files <- sequences(trim.dir=trim.dir)
    trim1 <- sequencing_files[[3]]
    trim2 <- sequencing_files[[4]]

    
    for (i in seq_along(sequencing_files[[5]])) {
        # kallisto quant command
        kallisto(trim1 = trim1[i],
                 trim2 = trim2[i],
                 quant.dir =  file.path(quant.dir, sequencing_files[[5]][i]),
                 ref.transcriptome = ref.transcriptome,
                 threads = threads,
                 i = i)
        
        # Message
        message(paste0('--- Done sample, ', sequencing_files[[5]][i], ' Time is: ', Sys.time(), ' ---'))
    }
    
    # Create TxDb SQL if not present
    if(!file.exists(txdb.file)) {
        ref_txb(txdb.dir = txdb.dir,
                txdb.file = txdb.file)
    }
    
    # tximport for all
    tximp_counts(txdb.dir = txdb.dir,
                 quant.dir = quant.dir,
                 out.dir = out.dir)
    
}