#CNV functions

#libraries needed: CaSpER, GenomicRanges, S4Vectors & IRanges

#==================================
###Cytoband information
#==================================
# Cytoband downloaded from:"https://hgdownload.cse.ucsc.edu/goldenPath/hg38/database/cytoBand.txt.gz"
# USE UNZIP FILE!
cytoband <-  function(cytoband.dir = 'data/cnv/cytoband/cytoBand.txt'){
    
    #Create cytoband information file
    cytoband = read.delim(cytoband.dir, header=F)
    cytoband = data.frame(V1=gsub("chr", "", cytoband[,1]), V2=cytoband[,2], V3=cytoband[,3], V4=substring(cytoband$V4, 1, 1), stringsAsFactors=F)
    start = do.call(rbind, lapply(split(cytoband$V2, paste0(cytoband$V1, cytoband$V4)), min))
    end = do.call(rbind, lapply(split(cytoband$V3, paste0(cytoband$V1, cytoband$V4)), max))
    cytoband = data.frame(V1=gsub("p", "", gsub("q", "", rownames(start))), V2=start, V3=end, V4=rownames(start), stringsAsFactors=F)
    cytoband = cytoband [as.vector(unlist(sapply(c(1:22, "X"), function(x) which(cytoband$V1 %in% x)))), ]
    cytoband$V4[grep("q", cytoband$V4)] <- "q"
    cytoband$V4[grep("p", cytoband$V4)] <- "p"
    rownames(cytoband) = NULL
    write.table(cytoband,"data/cnv/cytoband/cytoband.txt",sep="\t",row.names=FALSE, col.names=TRUE)
    
    # Print message
    message(paste0('Done cytoband, Time is: ', Sys.time()))
    return('')
}


#==================================
###Centromere information
#==================================
# Centromere information downloaded from: https://hgdownload.cse.ucsc.edu/goldenPath/hg38/database/cytoBand.txt.gz
# USE UNZIP FILE!
centromere <-  function(centromere.dir = 'data/cnv/centromere/cytoBand.txt'){
    
    # Create centromere.txt cmd
    cmd <- paste0('grep acen ', centromere.dir, ' >> data/cnv/centromere/centromere.txt')
    system(cmd)
    
    # Print message
    message(paste0('Done centromere, Time is: ', Sys.time()))
    return('')
}


#==================================
###CaSpER step
#==================================
casper <- function(cnv.dir = 'out/cnv',
                   cyto.dir = 'data/cnv/cytoband/cytoband.txt',
                   centro.dir = 'data/cnv/centromere/centromere.txt',
                   counts.dir = 'out/gene_counts',
                   sample.dir ='out/gene_counts/htseq_counts.tsv'){
    
    # Create output directory
    if (!file.exists(cnv.dir)) {
        dir.create(cnv.dir, recursive = TRUE)
    }
    
    #Import Gene expression counts
    counts_df = read.csv(sample.dir, sep = '\t',header = TRUE)
    
    #Modify gene IDs and rownames for rannotation file
    counts_df$GENE_ID = sapply(strsplit(counts_df$GENE_ID,".",fixed = T), `[`, 1)
    write.table(counts_df, file.path(counts.dir, 'transfo_counts.tsv'),row.names =F, quote = F,sep = '\t')
    row_ids = read.delim(file.path(counts.dir,'transfo_counts.tsv'), header = TRUE)
    df = data.matrix(row_ids)
    row.names(df) = row_ids$GENE_ID
    df = df[,2:9]
    
    #Generate annotation & match to counts dataframe
    centromere = read.delim(centro.dir, header=FALSE)
    
    rannotation = CaSpER::generateAnnotation(id_type="ensembl_gene_id", genes=rownames(df),
                                             ishg19=T, centromere=centromere, host="https://www.ensembl.org/")
    counts_df = df[match(rannotation$Gene,rownames(df)), ]

    #Samples IDs (IMPORTANT!!!colnames(control.sample.ids) == colnames(counts_df))
    control.sample.ids = colnames(counts_df)
    
    #Cytoband file modif
    cytoband = read.delim(cyto.dir, header=TRUE)
    
    #Create CaSpER object
    object = CaSpER::CreateCasperObject(raw.data=counts_df, loh.name.mapping=NULL, sequencing.type="bulk",
                                        cnv.scale=3, loh.scale=3, matrix.type="raw", expr.cutoff=4.5,
                                        annotation=rannotation, method="iterative", loh=NULL, filter="median",
                                        control.sample.ids=control.sample.ids, cytoband=cytoband, genomeVersion="hg38")
    #Run CaSpER object
    final.objects <- list()
    loh.list <- list()
    cnv.list <- list()
    
    message("Performing HMM segmentation...")
    
    for (i in 1:object@cnv.scale) {
        cnv.list[[i]] <- CaSpER::PerformSegmentationWithHMM(object, cnv.scale = i, removeCentromere = F, cytoband = cytoband)
    }
    
    for (i in 1:3) {
        
        object <- cnv.list[[i]]
        object@segments$states2 <- rep("neut", length(object@segments$state))
        
        object@segments$states2[as.numeric(as.character(object@segments$state)) == 1] <- "del"
        object@segments$states2[as.numeric(as.character(object@segments$state)) == 5] <- "amp"
        final.objects[[i]] <- CaSpER::generateLargeScaleEvents(object)
    }
    
    #Large-Scale CNV Summarization.
    finalChrMat <- CaSpER::extractLargeScaleEvents(final.objects, thr=0.75)
    write.table(finalChrMat, file.path(cnv.dir, "cnv_sum.tsv"),row.names=TRUE)
    
    #Segment based CNV summarization
    all.segments <- do.call(rbind, lapply(final.objects, function(x) x@segments))
    write.table(all.segments, file.path(cnv.dir, "all_segments.tsv"), row.names=FALSE)
    gamma <- 1 #Set gamma parameter for CNV calls
    suppressWarnings(segment.summary <- CaSpER::extractSegmentSummary(final.objects))
    loss <- segment.summary$all.summary.loss
    gain <- segment.summary$all.summary.gain
    loss.final <- loss[loss$count>=gamma, ]
    gain.final <- gain[gain$count>=gamma, ]
    
    #Gene based CNV Summarization
    all.summary<- rbind(loss.final, gain.final)
    colnames(all.summary) [2:4] <- c("Chromosome", "Start",   "End")
    rna <- GenomicRanges::GRanges(seqnames = S4Vectors::Rle(gsub("q", "", gsub("p", "", all.summary$Chromosome))),
                                  IRanges::IRanges(all.summary$Start, all.summary$End))   
    ann.gr <- GenomicRanges::makeGRangesFromDataFrame(final.objects[[1]]@annotation.filt, keep.extra.columns = TRUE, seqnames.field="Chr")
    hits <- IRanges::findOverlaps(rna, ann.gr)
    genes <- CaSpER::splitByOverlap(ann.gr, rna, "GeneSymbol")
    genes.ann <- lapply(genes, function(x) x[!(x=="")])
    all.genes <- unique(final.objects[[1]]@annotation.filt[,2])
    all.samples <- unique(as.character(final.objects[[1]]@segments$ID))
    rna.matrix <- CaSpER::gene.matrix(seg=all.summary, all.genes=all.genes, all.samples=all.samples, genes.ann=genes.ann)
    write.table(rna.matrix, file=file.path(cnv.dir, "rna_matrix.tsv"))
    
    # Print message
    message(paste0('Done CaSpER, Time is: ', Sys.time()))
}


#======================
# cleanup 
#======================
cleanup <-  function(counts.dir = 'out/gene_counts')
{
    # rm command
    rm_cmd = paste0('rm ',file.path(counts.dir,'*transfo_counts.tsv*'))
    
    system(rm_cmd)
    
    # Print message
    message(paste0('Done cleanup of temporary transfo_counts (tsv), Time is: ',Sys.time()))
    return('')
}
 
   
#=====================
#Wrapper: cnv_rnaseq
#=====================
cnv_rnaseq <- function(cyto.dir = 'data/cnv/cytoband/cytoband.txt',
                       cytoband.dir = 'data/cnv/cytoband/cytoBand.txt',
                       centro.dir = 'data/cnv/centromere/centromere.txt',
                       centromere.dir = 'data/cnv/cytoband/cytoBand.txt',
                       counts.dir = 'out/gene_counts',
                       cnv.dir = 'out/cnv')
    suppressWarnings(
    {
    
    #Create cytoband info if not present
    if (!file.exists(cyto.dir)){
        cytoband(cytoband.dir = cytoband.dir)}
    
    #Create centromere info if not present
    if (!file.exists(centro.dir)){
        centromere(centromere.dir = centromere.dir)}
    
    #Name of count files
    list_files <- list.files(counts.dir)
    
    for (i in seq_along(list_files)){
        #CaSpER command
        casper(cnv.dir = cnv.dir,
               cyto.dir = cyto.dir,
               centro.dir = centro.dir,
               counts.dir = counts.dir,
               sample.dir = file.path(counts.dir,list_files[[1]]))
        
        # Message
        message(paste0('--- Done sample, ', list_files[[1]], ' Time is: ', Sys.time(), ' ---'))
        
    }
    
    # Print message
    message(paste0('Done cnv, Time is: ', Sys.time(),
                   '\nNote: In result matrix, 0 = neutral, 1 = amplification & -1 = deletion'))
    # Cleanup
    cleanup(counts.dir = counts.dir)
}
)