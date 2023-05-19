source('R/counts_functions.R')
#kallisto function
library(data.table)
library(GenomicFeatures)
library(tximport)
library(DESeq2)
library(apeglm)
library(ggplot2)
library(ggrepel)
library(rhdf5)
library(GenomicState)
library(EnhancedVolcano)
#=====================
###get sequences ready
#=====================
sequences()

#=====================
###trimming
#=====================
trimming()
#=====================
#Pseudo-mapping step
#=====================
kallisto = function(R1_trim = sequences()[[3]][1],
                   R2_trim = sequences()[[4]][1],
                   in.dir = file.path('data/index/'),
                   out.dir = file.path('data/kallisto/'),
                   out_prefix = paste0('rnaseq/out/',sequences()[[5]][1],'_')
                   ) {
    
    #index cmd
    # index downloaded from here: https://www.gencodegenes.org/human/
    #Crééer fichier contenant l'index
    index = paste0(ifelse(Sys.info()['sysname'] == 'Windows','wsl.exe ',''),'kallisto index -i ', in.dir, 'gencode.v43.transcripts.idx ', in.dir, 'gencode.v43.transcripts.fa.gz')
    system(index)
    
    #kallisto cmd
    cmd = paste0(ifelse(Sys.info()['sysname'] == 'Windows','wsl.exe ',''),'kallisto quant -i ', in.dir, 'gencode.v43.transcripts.idx -o ',
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
    
    return('')
}

#=====================
#Data infos
#=====================
#Text file of run infos
sampleTable = read.delim("data/sampleTable.txt")
rownames(sampleTable) = sampleTable$Name

#=====================
#load TxDB from gencode (Release 43 (GRCh38.p13))
#=====================
txdb_v31_hg19_chr1 = gencode_txdb("43", "hg19", chrs = "chr1") 
#chrs = paste0("chr", c(seq_len(22), "X", "Y", "M"))
k = keys(txdb_v31_hg19_chr1, keytype = "TXNAME")
tx2gene = select(txdb_v31_hg19_chr1, k, "GENEID", "TXNAME")
tx2gene$TXNAME = sapply(strsplit(tx2gene$TXNAME,".",fixed = T), `[`, 1)

#=====================
#Load expression data
#=====================
dir = "data/"
files = file.path(paste0(dir, "kallisto"), list.files(paste0(dir, "kallisto")), "abundance.h5")
names(files) = list.files(paste0(dir, "kallisto"))

#import abundance.h5 files
txi = tximport(files, type = "kallisto", tx2gene = tx2gene, 
               txIn = TRUE, txOut = FALSE, countsFromAbundance = "lengthScaledTPM",  ignoreTxVersion = TRUE, ignoreAfterBar=TRUE)

#=====================
#Create DESeqDataSet object
#=====================
dds = DESeq2::DESeqDataSetFromTximport(txi, colData = sampleTable)
#ERROR if (all(assay(se) == 0)) {
#stop("all samples have 0 counts for all genes. check the counting script.")

#Differential expression analysis
dds = DESeq(dds)
res = results(dds, name = "Sain vs tumeur", alpha = 0.05)

#=====================
#Generating plots
#=====================
#count data transformation with variance stabilizing transformations (vst)
variance = vst(dds)
# PCA_plot (Principal Componenet Analysis)
pcaData = plotPCA(vst, intgroup=c("Name","Sain vs tumeur"), 
                  returnData=TRUE)
percentVar = round(100 * attr(pcaData, "percentVar"))

png("DGE_PCA-vst.kallisto.png", width=7, height=7, units = "in", res = 300)
ggplot2::ggplot(pcaData, aes(PC1, PC2, colour = paris_classification)) + 
    geom_point(size = 2) +theme_bw() + scale_color_manual(values = c("blue", "red")) +
    geom_text_repel(aes(label = individual), nudge_x = -1, nudge_y = 0.2, size = 3) +
    ggtitle("Principal Component Analysis (PCA)", subtitle = "vst transformation") +
    xlab(paste0("PC1: ",percentVar[1],"% variance")) +
    ylab(paste0("PC2: ",percentVar[2],"% variance"))
dev.off()

#Volcano_plot for DGE
pCutoff = 0.05
FCcutoff = 1.0

p = EnhancedVolcano(data.frame(res), lab = NA, x = 'log2FoldChange', y = 'padj',
                    xlab = bquote(~Log[2]~ 'fold change'), ylab = bquote(~-Log[10]~adjusted~italic(P)),
                    pCutoff = pCutoff, FCcutoff = FCcutoff, pointSize = 1.0, labSize = 2.0,
                    title = "Volcano plot", subtitle = "SSA/P vs. Normal",
                    caption = paste0('log2 FC cutoff: ', FCcutoff, '; p-value cutoff: ', pCutoff, '\nTotal = ', nrow(res), ' variables'),
                    legend=c('NS','Log2 FC','Adjusted p-value', 'Adjusted p-value & Log2 FC'),
                    legendPosition = 'bottom', legendLabSize = 14, legendIconSize = 5.0)

png("DGE_VolcanoPlots.kallisto.png", width=7, height=7, units = "in", res = 300)
print(p)
dev.off()




#==============================================
# OTHER ALTERNATIVE
#==============================================


# Script from : https://pachterlab.github.io/sleuth_walkthroughs/trapnell/analysis.html
#=====================
#Path kallisto results
#=====================
dir = "/data/"
sample_id <- dir(file.path("/data/kallisto"))
kal_dirs <- file.path("/data/kallisto", sample_id)
kal_dirs

#=====================
#SampleTable (infos samples assiciated with kallisto quant)
#=====================
s2c = read.table(file.path("/data/", "sampleTable.txt"), header = TRUE, stringsAsFactors=FALSE)
s2c = dplyr::select(s2c, sample = Name, Type)
s2c

s2c = dplyr::mutate(s2c, path = kal_dirs)

#=====================
#Initialize sleuth object
#=====================
so = sleuth_prep(s2c, extra_bootstrap_summary = TRUE)

#=====================
#Including gene names for transcript-level analysis
#=====================
#Collect gene names
mart = biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                        dataset = "hsapiens_gene_ensembl",
                        host = 'ensembl.org')
t2g = biomaRt::getBM(attributes = c("ensembl_transcript_id","ensembl_transcript_id_version", "ensembl_gene_id",
                                    "ensembl_gene_id_version","external_gene_name","description",
                                    "chromosome_name","start_position",
                                    "end_position","strand",
                                    "entrezgene_id"), mart = mart)
t2g = dplyr::rename(t2g, target_id = ensembl_transcript_id,
                    ens_gene = ensembl_gene_id, ext_gene = external_gene_name)

t2g = dplyr::select(t2g, c('target_id', 'ens_gene', 'ext_gene'))

#Adding gene names to sleuth table
so = sleuth_prep(s2c, target_mapping = t2g)

#fitting measurement error models
so = sleuth_fit(so, ~1, 'reduced')

so = sleuth_lrt(so, 'reduced', 'reduced')

#View results
sleuth_live(so)

#PCA plot
plot_pca(so, color_by = 'condition')

#Plot group density for count distribution for each sample
plot_group_density(so, use_filtered = TRUE, units = "est_counts",
                   trans = "log", grouping = setdiff(colnames(so$sample_to_covariates),
                                                     "sample"), offset = 1)