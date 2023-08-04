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
    
    cmd1 <- paste0('cd ', ctat.dir)
    cmd2 <- paste0('wget ', ctat.lib)
    cmd3 <- paste0('tar -xvf GRCh38_gencode_v37_CTAT_lib_Mar012021.plug-n-play.tar.gz')
    
    system(cmd1);system(cmd2);system(cmd3)
    
    #message(cmd)
    message(paste0('Done creating CTAT Genome library, Time is: ',Sys.time()))
    
}

#=====================
#STAR-Fusion
#=====================

fusion <- function(ctat.dir = 'data/ctat/GRCh38_gencode_v37_CTAT_lib_Mar012021.plug-n-play/ctat_genome_lib_build_dir',
                   trim1 = sequences()[[3]][1],
                   trim2 = sequences()[[4]][1],
                   fus.dir = 'out/fusion',
                   CPU = 10,
                   i=1){
    
    # Create output directory
    if (!file.exists(fus.dir)) {
        dir.create(fus.dir, recursive = TRUE)
    }
    
    cmd <- paste0(ifelse(Sys.info()['sysname'] == 'Windows','wsl.exe ',''),
                  'STAR-Fusion --left_fq ', trim1, ' --right_fq ', trim2,
                 ' --genome_lib_dir ', ctat.dir,
                 ' --CPU ', CPU,' --output_dir ', fus.dir, 
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
clean_files = function (fus.dir = 'out/fusion')
{
    rm_cmd1 = paste0('rm ',fus.dir,'/*/*Aligned.out.bam*')
    rm_cmd2 = paste0('rm ',fus.dir,'/*/*.cmds')
    rm_cmd3 = paste0('rm -r ',fus.dir,'/*/*_starF_checkpoints')
    rm_cmd4 = paste0('rm -r ',fus.dir,'/*/*tmp_chim_read_mappings_dir')
    rm_cmd5 = paste0('rm -r ',fus.dir,'/*/*star-fusion.preliminary')
    
    system(rm_cmd1);system(rm_cmd2);system(rm_cmd3);system(rm_cmd4);system(rm_cmd5)
    
    message(paste0('Done cleanup of output & directory files because they are really big, Time is: ',Sys.time()))
    return('')
}


#=====================
#Wrapper: fusion_rnaseq
#=====================

fusion_rnaseq <- function(ctat.dir = 'data/ctat',
                          ctat.lib = 'https://data.broadinstitute.org/Trinity/CTAT_RESOURCE_LIB/GRCh38_gencode_v37_CTAT_lib_Mar012021.plug-n-play.tar.gz',
                          trim.dir = 'out/fastq.trim',
                          fus.dir = 'out/fusion',
                          CPU = 10){
    
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
       fusion(ctat.dir = file.path(ctat.dir,'GRCh38_gencode_v37_CTAT_lib_Mar012021.plug-n-play/ctat_genome_lib_build_dir'),
              trim1 = trim1[i],
              trim2 = trim2[i],
              fus.dir = file.path(fus.dir, sequencing_files[[5]][i]),
              CPU = CPU,
              i=i)
        
        # Message
        message(paste0('--- Done sample, ', sequencing_files[[5]][i], ' Time is: ', Sys.time(), ' ---'))
        
    }
    
    #final cleanup
    clean_files(fus.dir = fus.dir)
    
}