# kallisto function
source('R/counts_functions.R')

# Index command
index <- function(ref.transcriptome = 'data/reference_transcriptome/chr1.fasta.gz') { #reference transcriptome: full transcriptome can be found here: https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_43/gencode.v43.transcripts.fa.gz

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

# Kallisto command
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

# Wrapper: kallisto_rnaseq
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
        #Kallisto quant command
        kallisto(trim1 = trim1[i],
                 trim2 = trim2[i],
                 quant.dir =  file.path('out/kallisto', sequencing_files[[5]][i]),
                 ref.transcriptome = ref.transcriptome,
                 i = i)
        
        # Message
        message(paste0('--- Done sample, ', sequencing_files[[5]][i], ' Time is: ', Sys.time(), ' ---'))
    }
}
