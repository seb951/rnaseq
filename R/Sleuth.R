#Sleuth
source('R/kallisto.R')

gene_level = function(in.dir = "out/kallisto")
    
{
    #=====================
    #Import kallisto data
    #=====================
    
    #load TxDB from gencode (https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_43/gencode.v43.annotation.gtf.gz)
    gtf = file.path("test/reference_genome/gencode.v33.annotation.gtf")
    txdb.filename = file.path("test/gencode.v43.annotation.sqlite")
    txdb = GenomicFeatures::makeTxDbFromGFF(gtf)
    
    #Save TxDb database (SQLite database) for later uses
    saveDb(txdb, txdb.filename)
    
    #load to use the TxDb database
    txdb = AnnotationDbi::loadDb(txdb.filename)
    
    
    #load TxDB & create tx2gene
    k = AnnotationDbi::keys(txdb, keytype = "TXNAME")
    tx2gene = AnnotationDbi::select(txdb, k, "GENEID", "TXNAME")
    #tx2gene$TXNAME = sapply(strsplit(tx2gene$TXNAME,".",fixed = T), `[`, 1)
    
    #Load kallisto data
    test.dir = "test"
    files = file.path(paste0(test.dir),list.files(test.dir, pattern = "_R1"), "abundance.tsv")
    names(files) = list.files(test.dir, pattern = "_R1")
    
    #with tsv
    txi_tsv = tximport::tximport(files, type = "kallisto", tx2gene = tx2gene, ignoreAfterBar = TRUE)
    #head(txi_tsv$counts)
    
    #Text file of run info
    sampleTable = read.delim(file.path("data/sampleTable.txt"))
    rownames(sampleTable) = sampleTable$Name
    sampleTable$Type = as.factor(sampleTable$Type)
    sampleTable$Sample = as.factor(sampleTable$Sample)
    
    #kallisto gene_level table
    dds_kal = DESeq2::DESeqDataSetFromTximport(txi_tsv, colData = sampleTable, ~ Type)
    dds_deseq = DESeq2::DESeq(dds_kal)
    est_dds = estimateSizeFactors(dds_deseq)
    dds_counts = counts(est_dds, normalized = TRUE)
    write.table(dds_counts, file = "test/kallisto_counts.tsv", sep="\t")

    #PCA plot
    log_kal = DESeq2::rlog(dds_deseq)
    log_kal = DESeq2::rlog(dds_deseq)
    pcaData = plotPCA(r_log, intgroup=c("Type"), 
                      returnData=TRUE)
    percentVar = round(100 * attr(pcaData, "percentVar"))
    
    pdf("DGE_PCA-rlog.kallistio.pdf", width=7, height=7, units = "in", res = 300)
    ggplot2::ggplot(pcaData, ggplot2::aes(PC1, PC2, colour = Type)) + 
        ggplot2::geom_point(size = 2) + ggplot2::theme_bw() + 
        ggplot2::scale_color_manual(values = c("blue", "red")) +
        ggrepel::geom_text_repel(ggplot2::aes(label = Type), nudge_x = -1, nudge_y = 0.2, size = 3) +
        ggplot2::ggtitle("Principal Component Analysis (PCA)", subtitle = "rlog transformation") +
        ggplot2::xlab(paste0("PC1: ",percentVar[1],"% variance")) +
        ggplot2::ylab(paste0("PC2: ",percentVar[2],"% variance"))
    dev.off()
    
    
}
    
    
#=====================
#Data visualisation
#=====================  
sleuth = function()
{     
    #abundance
    kal_files = file.path(paste0(test.dir),list.files(test.dir, pattern = "_R1"), "*.tsv")
    names(kal_files) = list.files(test.dir, pattern = "_R1")
    
    #creat csv file
    readr::read_tsv(kal_files)
    colnames(htseq_count) = c('features',samples)
    write.table(htseq_count,file.path(out.dir,'htseq_counts.tsv'),row.names =F, quote = F,sep = '\t')
    
    
    
    #Grep info from run file csv & convert to txt
    #data_file = readxl::read_xlsx("/data/*.xlsx", sheet = 1, col_names = TRUE, col_types = NULL, na = "", skip = 0)
    #file.dir = file.path("/data/")
    #data_frame = write.csv(data_file, file.path(file.dir, "test.csv"))
    #df = read.csv("/data/test.csv")
    #col = df[, c("Name...2", "Sample...5", "Type")] > "/data/sampleTable.txt.txt
    

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