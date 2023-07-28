# kallisto function
source('R/counts_functions.R')

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
}


#=====================
# kallisto command
#=====================
kallisto <- function(trim1 = 'out/fastq.trim/RNA_0006_5598_Tumeur_R1_val_1.fq.gz',
                     trim2 = 'out/fastq.trim/RNA_0006_5598_Tumeur_R1_val_1.fq.gz',
                     quant.dir = 'out/kallisto',
                     ref.transcriptome = 'data/reference_transcriptome/chr1.fasta.gz',
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
                  ' 1>',ifelse(i>1,'>',''),
                  file.path('out/logs', 'kallisto.out'),
                  ' 2>',ifelse(i>1,'>',''),
                  file.path('out/logs', 'kallisto.err'))
    
    system(cmd)
    
    # Print message
    message(paste0('Done Kallisto count, Time is: ', Sys.time()))
}


#==========================================
# tximport counts step
#==========================================
# gtf downloaded from here: https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_43/gencode.v43.annotation.gtf.gz
tximp_counts = function(gtf.dir = "data/reference_genome/gencode.v33.annotation.gtf",
                        txdb.dir = "data/gencode.v43.annotation.sqlite",
                        quant.dir = "out/kallisto",
                        out.dir = "out/gene_counts"){
    
    # Create TxDb database and save to SQLite for later uses
    txdb = GenomicFeatures::makeTxDbFromGFF(gtf.dir)
    saveDb(txdb.dir, txdb.filename)
    
    # Load TxDB & create tx2gene
    txdb = AnnotationDbi::loadDb(txdb.filename)
    k = AnnotationDbi::keys(txdb, keytype = "TXNAME")
    tx2gene = AnnotationDbi::select(txdb, k, "GENEID", "TXNAME")
    #tx2gene$TXNAME = sapply(strsplit(tx2gene$TXNAME,".",fixed = T), `[`, 1)
    
    # Import kallisto data
    files = file.path(paste0(kal.dir),list.files(test.dir, pattern = "_R1"), "abundance.tsv")
    names(files) = list.files(kal.dir, pattern = "_R1")
    txi_tsv = tximport::tximport(files, type = "kallisto", tx2gene = tx2gene, ignoreAfterBar = TRUE)
    #head(txi_tsv$counts)
    
    # Samples information text file (Contains: Sample name, type(Tumeur/Sain) & sample number)
    sampleTable = read.delim(file.path("data/sampleTable.txt"))
    rownames(sampleTable) = sampleTable$Name
    sampleTable$Type = as.factor(sampleTable$Type)
    sampleTable$Sample = as.factor(sampleTable$Sample)
    
    # Generate kallisto gene counts intto dataframe
    dds_kal = DESeq2::DESeqDataSetFromTximport(txi_tsv, colData = sampleTable, ~ Type)
    dds_deseq = DESeq2::DESeq(dds_kal)
    dds_counts = counts(dds_deseq, normalized = TRUE)
    dds_counts = round(dds_counts)
    write.table(dds_counts, file.path(out.dir,'kallisto_counts.tsv'), sep="\t")
    kal_count = read.csv(file.path(out.dir,'kallisto_counts.tsv'),sep = '\t',header = F)
    kal_count = kal_count[- 1, ] 
    colnames(kal_count) = c('GENE_ID',names(files))
    write.table(kal_count,file.path(out.dir,'kallisto_counts.tsv'),row.names =F, quote = F,sep = '\t')
}


#==========================================
# Wrapper: kallisto_rnaseq
#==========================================
kallisto_rnaseq <- function(idx.dir ='data/index',
                            ref.transcriptome = 'data/reference_transcriptome/chr1.fasta.gz',
                            trim.dir = 'out/fastq.trim',
                            quant.dir = 'out/kallisto') {
  
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
                 i = i)
        
        # Message
        message(paste0('--- Done sample, ', sequencing_files[[5]][i], ' Time is: ', Sys.time(), ' ---'))
    }
    
    # tximport for all
    tximp_counts(gtf.dir = gtf.dir,
                 txdb.dir = txdb.dir,
                 quant.dir = quant.dir,
                 out.dir = out.dir)
    
}
