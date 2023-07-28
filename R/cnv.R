#CNV functions

#==================================
###Cytoband information
#==================================
#Cytiband downloaded from:"https://hgdownload.cse.ucsc.edu/goldenPath/hg38/database/cytoBand.txt.gz"
cytoband <-  function(cytoband.dir = 'data/cnv/cytoband/cytoBand.txt',
                      cyto.ref= "https://hgdownload.cse.ucsc.edu/goldenPath/hg38/database/cytoBand.txt.gz"){
    
    # Create output directory
    if (!file.exists(cytoband.dir)) {
        dir.create(cytoband.dir, recursive = TRUE)
    }
    
    #Cytoband cmd
    cmd <- paste0('wget ', cyto.ref, ' | gunzip cytoBand.txt.gz' )
    system(cmd)
    
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
    
    # Print message
    message(paste0('Done cytoband, Time is: ', Sys.time()))
}


#==================================
###Centromere information
#==================================
#Centromere information downloaded from: https://hgdownload.cse.ucsc.edu/goldenPath/hg38/database/cytoBand.txt.gz
centromere <-  function(centromere.dir = 'data/cnv/centromere/centromere.txt',
                        centro.ref= "https://hgdownload.cse.ucsc.edu/goldenPath/hg38/database/cytoBand.txt.gz"){
    
    # Create output directory
    if (!file.exists(centromere.dir)) {
        dir.create(centromere.dir, recursive = TRUE)
    }
    
    #Centromere cmd
    cmd <- paste0('curl -s ', centro.ref,' | gunzip cytoBand.txt.gz | grep acen >> ', centromere.dir)
    system(cmd)
    
    # Print message
    message(paste0('Done centromere, Time is: ', Sys.time()))
}


#==================================
###CaSpER step
#==================================
casper <- function(out.dir = 'out/cnv',
                   cytoband.dir = 'data/cnv/cytoband/cytoBand.txt',
                   centromere.dir = 'data/cnv/centromere/centromere.txt',
                   counts.dir = 'out/gene_counts/htseq_counts.tsv'){
    
    # Create output directory
    if (!file.exists(out.dir)) {
        dir.create(out.dir, recursive = TRUE)
    }
    
    #Import Gene expression counts
    counts_df = read.delim('data/gene_counts/htseq_counts.tsv', header = TRUE)
    
    #Modify gene IDs and rownames for rannotation file
    counts_df$GENE_ID = sapply(strsplit(counts_df$GENE_ID,".",fixed = T), `[`, 1)
    write.table(counts_df,file.path('data/cnv/transfo_counts_gene/transfo_counts.tsv'),row.names =F, quote = F,sep = '\t')
    row_ids = read.delim('data/cnv/transfo_counts_gene/transfo_counts.tsv', header = TRUE)
    df = data.matrix(row_ids)
    row.names(df) = row_ids$GENE_ID
    
    #Generate annotation & match to counts dataframe
    centromere = read.delim(centromere.dir, header=FALSE)
    
    rannotation = CaSpER::generateAnnotation(id_type="ensembl_gene_id", genes=rownames(df),
                                             ishg19=T, centromere=centromere, host="https://www.ensembl.org/")
    counts_df = df[match( rannotation$Gene,rownames(df)), ]

    #Samples IDs (IMPORTANT!!!colnames(control.sample.ids) == colnames(data_fr))
    control.sample.ids = colnames(counts_df)
    
    #Create CaSpER object
    counts_df = read.delim('data/counts_gene/htseq_counts.tsv', header=TRUE)
    
    object = CaSpER::CreateCasperObject(raw.data=counts_df, loh.name.mapping=NULL, sequencing.type="bulk",
                                        cnv.scale=3, loh.scale=3, matrix.type="normalized", expr.cutoff=4.5,
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
    finalChrMat <- CaSpER::extractLargeScaleEvents (final.objects, thr=0.75)
    
    #Segment based CNV summarization
    gamma <- 6
    all.segments <- do.call(rbind, lapply(final.objects, function(x) x@segments))
    segment.summary <- CaSpER::extractSegmentSummary (final.objects)
    loss <- segment.summary$all.summary.loss
    gain <- segment.summary$all.summary.gain
    loh <- segment.summary$all.summary.loh
    loss.final <- loss[loss$count>gamma, ]
    gain.final <- gain[gain$count>gamma, ]
    loh.final <- loh[loh$count>gamma, ]
    
    #Gene based CNV Summarization
    all.summary<- rbind(loss.final, gain.final)
    colnames(all.summary) [2:4] <- c("Chromosome", "Start",   "End")
    rna <- GenomicRanges::GRanges(seqnames = Rle(gsub("q", "", gsub("p", "", all.summary$Chromosome))),
                   IRanges::IRanges(all.summary$Start, all.summary$End))   
    ann.gr <- GenomicRanges::makeGRangesFromDataFrame(final.objects[[1]]@annotation.filt, keep.extra.columns = TRUE, seqnames.field="Chr")
    hits <- IRanges::findOverlaps(rna, ann.gr)
    genes <- CaSpER::splitByOverlap(ann.gr, rna, "GeneSymbol")
    genes.ann <- lapply(genes, function(x) x[!(x=="")])
    all.genes <- unique(final.objects[[1]]@annotation.filt[,2])
    all.samples <- unique(as.character(final.objects[[1]]@segments$ID))
    rna.matrix <- CaSpER::gene.matrix(seg=all.summary, all.genes=all.genes, all.samples=all.samples, genes.ann=genes.ann)
    write.table(rna.matrix,file.path("out/cnv/rna_matrix.tsv"))
    
    # Print message
    message(paste0('Done CaSpER, Time is: ', Sys.time(),
                   '\nNote: In result matrix, 0 = alteration, 1 = amplification & -1 = deletion'))
}
 
   
#=====================
#Wrapper: cnv_rnaseq
#=====================
cnv_rnaseq <- function(cytoband.dir = params$cytoband.dir,
                       centromere.dir = params$centromere.dir,
                       counts.dir = params$counts.dir,
                       out.dir = params$out.dir){
    
    #Create cytoband info if not present
    if (!file.exists(cytoband.dir)){
        cytoband(cytoband.dir = cytoband.dir,
                 cyto.ref= cyto.ref)}
    
    #Create centromere info if not present
    if (!file.exists(centromere.dir)){
        centromere(centromere.dir = centromere.dir,
                   centro.ref= centro.ref)}
    
    #Name of count files
    counts_file <- list.files(counts.dir)
    
    for (i in seq_along(counts_file)){
        #CaSpER command
        casper(out.dir = out.dir,
               cytoband.dir = cytoband.dir,
               centromere.dir = centromere.dir,
               centromere = centromere,
               counts_df = counts_file[i])
        
        # Message
        message(paste0('--- Done sample, ', counts_file[i], ' Time is: ', Sys.time(), ' ---'))
        
    }
}