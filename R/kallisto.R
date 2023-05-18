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
#load TxDB from https://www.gencodegenes.org/human/
#=====================
txdb = gencode_txdb("43", "hg19", chrs = paste0("chr", c(seq_len(22), "X", "Y", "M")))
k = keys(txdb, keytype = "TXNAME")
tx2gene = select(txdb, k, "GENEID", "TXNAME")


#=====================
#Load expression data
#=====================
dir = ('data/')
files = file.path(paste0(dir, "kallisto"), list.files(paste0(dir, "kallisto")), "abundance.h5")
names(files) = list.files(paste0(dir, "kallisto"))

#import abundance.h5 files
txi = tximport(files, type = "kallisto", tx2gene = tx2gene, 
               txIn = TRUE, txOut = FALSE, countsFromAbundance = "no",  ignoreTxVersion = TRUE, ignoreAfterBar = TRUE)

#=====================
#Create DESeqDataSet object
#=====================
dds = DESeqDataSetFromTximport(txi, 
                               colData = sampleData, ~ individual + paris_classification)

#Differential expression analysis
dds = DESeq(dds)
res <- results(dds, name = "paris_classification_0.IIa_vs_normal", alpha = 0.05)

#=====================
#Generating plots
#=====================
#count data transformation with variance stabilizing transformations (vst)
variance = vst(dds)
# PCA_plot (Principal Componenet Analysis)
pcaData = plotPCA(vsd, intgroup=c("individual","paris_classification"), 
                  returnData=TRUE)
percentVar = round(100 * attr(pcaData, "percentVar"))

png("DGE_PCA-vst.Salmon.png", width=7, height=7, units = "in", res = 300)
ggplot(pcaData, aes(PC1, PC2, colour = paris_classification)) + 
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

png("DGE_VolcanoPlots.Salmon.png", width=7, height=7, units = "in", res = 300)
print(p)
dev.off()

#=====================
#Exporting results
#=====================
annoData = "data/refanno/gencode.v33.annotation_genes.txt"
annoData = data.frame(fread(annoData))

normCounts = as.data.frame(counts(dds, normalized = TRUE))
baseMeans = as.data.frame(sapply( levels(dds$paris_classification), 
                                  function(lvl) rowMeans( counts(dds, normalized = TRUE)[, dds$paris_classification == lvl, drop = FALSE] ) ))

normData = merge(annoData, merge(baseMeans, normCounts, by.x = 'row.names', by.y = 'row.names'), by.x = 'GeneID', by.y = 'Row.names')
normData = normData[order(normData$Chromosome, normData$Start, normData$End),]

deData = data.frame(res[,c(1,2,5,6)])
colnames(deData) = c("baseMean","log2fc","pvalue","padj")
deData = merge(annoData, deData, by.x = 'GeneID', by.y = 'row.names')
deData = deData[order(deData$Chromosome, deData$Start, deData$End),]

write.table(normData, file="DGE_DESeq2_Means_and_NormalisedCount.kallisto.txt", 
            sep = "\t", quote = F, row.names = F, col.names = T)

write.table(deData, file="DGE_DESeq2_DE_results.kallisto.txt", 
            sep = "\t", quote = F, row.names = F, col.names = T)