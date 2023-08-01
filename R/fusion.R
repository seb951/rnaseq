#fusion functions
source('R/counts_functions.R')

### USE LATEST VERSION OF STAR-Fusion (v1.12.0) ###

#=====================
#CTAT genome library
#=====================
genome_lib <- function(ctat.dir = 'data/ctat',
                       ctat.lib = 'https://data.broadinstitute.org/Trinity/CTAT_RESOURCE_LIB/GRCh38_gencode_v37_CTAT_lib_Mar012021.plug-n-play.tar.gz'){
    
    # Create output directory
    if (!file.exists(ctat.dir)) {
        dir.create(ctat.dir, recursive = TRUE)
    }
    
    cmd <- paste0('cd ', ctat.dir,' | wget ', ctat.lib,' | tar -xvf GRCh38_gencode_v37_CTAT_lib_Mar012021.plug-n-play.tar.gz')
    
    system(cmd)
    
    #message(cmd)
    message(paste0('Done creating CTAT Genome library, Time is: ',Sys.time()))
    
}

#=====================
#STAR-Fusion
#=====================

fusion <- function(ctat.dir = 'data/ctat/GRCh38_gencode_v37_CTAT_lib_Mar012021.plug-n-play/ctat_genome_lib_build_dir',
                   trim1 = 'out/fastq.trim/RNA_0006_5598_Tumeur_R1_val_1.fq.gz',
                   trim2 = 'out/fastq.trim/RNA_0006_5598_Tumeur_R2_val_2.fq.gz',
                   out.dir = 'out/fusion',
                   CPU = 10,
                   i=1){
    
    # Create output directory
    if (!file.exists(out.dir)) {
        dir.create(out.dir, recursive = TRUE)
    }
    
    cmd <- paste0('STAR-Fusion --left_fq ', trim1, ' --right_fq ', trim2,
                 ' --genome_lib_dir ', ctat.dir,
                 ' --CPU ', CPU,' --output_dir ', out.dir, 
                 ' 1> ',ifelse(i>1,'>',''),
                 file.path('out/logs', 'star-fusion.out'),
                 ' 2> ',ifelse(i>1,'>',''),
                 file.path('out/logs', 'star-fusion.err'))

    
    system(cmd)
    
    #message(cmd)
    message(paste0('Done STAR-Fusion, Time is: ',Sys.time()))
} 


#======================
# cleanup 
#======================
#review#
cleanup = function (out.dir = 'out/fusion')
{
    
    rm_cmd = paste0('rm -v ! ',out.dir,'star-fusion.fusion_predictions.tsv | ', out.dir, 'star-fusion.fusion_predictions.abridged.tsv')
    
    system(rm_cmd)
    
    message(paste0('Done cleanup of temporary (bam) because they are really big, Time is: ',Sys.time()))
    return('')
}


#=====================
#Wrapper: fusion_rnaseq
#=====================

fusion_rnaseq <- function(ctat.dir = params$ctat.dir,
                          ctat.lib = params$ctat.lib,
                          trim.dir = params$trim.dir,
                          out.dir = params$out.dir,
                          CPU = params$CPU){
    
    #Create CTAT Genome library if not present
    if (!file.exists(ctat.dir)) {
        genome_lib(ctat.dir = ctat.dir,
                   ctat.lib = ctat.lib)
    }
    
    # Name of sequencing files
    sequencing_files <- sequences(trim.dir=trim.dir)
    trim1 <- sequencing_files[[3]]
    trim2 <- sequencing_files[[4]]
    
    for (i in seq_along(sequencing_files[[5]])) {
        #STAR-Fusion command
       fusion(ctat.dir = 'data/ctat/GRCh38_gencode_v37_CTAT_lib_Mar012021.plug-n-play/ctat_genome_lib_build_dir',
              trim1 = trim1[i],
              trim2 = trim2[i],
              out.dir = file.path('out/fusion', sequencing_files[[5]][i]),
              CPU = CPU,
              i=i)
        
        # Message
        message(paste0('--- Done sample, ', sequencing_files[[5]][i], ' Time is: ', Sys.time(), ' ---'))
        
        }
}