suppressWarnings(library(optparse))
source('R/counts_functions.R')
source('R/counts_mini.R')
source('R/kallisto.R')
source('R/fusion.R')
source('R/cnv.R')

#For info, in Windows, you can call Rscript as such: C:'\Program Files\'R\R-4.2.2\bin\x64\Rscript.exe script.R

#options (QC, counts, kallisto, fusion and cnv)
option_list = list(
    make_option( "--fqdir",type="character", default='data/fastq2', 
                 help="fastq directory [default %default]", metavar="character"),
    
    make_option( "--trimdir",type="character", default='out/fastq.trim', 
                 help="trimmed fastq directory [default %default]", metavar="character"),
    
    make_option( "--genomeindexdir",type="character", default='data/reference_genome/chr1_small_index', 
                 help="reference genome index directory [default %default]", metavar="character"),
    
    make_option( "--annotationgtf",type="character", default='data/reference_genome/gencode.v43.primary_assembly.annotation_small.gtf', 
                 help="reference genome gtf file [default %default]", metavar="character"),
    
    make_option( "--genomefasta",type="character", default='data/reference_genome/chr1_small.fa', 
                 help="reference genome fasta file [default %default]", metavar="character"),
    
    make_option( "--outdir",type="character", default='out/', 
                 help="output directory [default %default]", metavar="character"),
    
    make_option( "--cutadapt",type="character", default='/usr/bin/cutadapt', 
                 help="cutadapt directory [default %default]", metavar="character"),
    
    make_option( "--threads",type="integer", default=12, 
                 help="Nb of Threads to use [default %default]", metavar="character"),
    
    make_option( "--nbfiles",type="character", default='all', 
                 help="specify how many files to process [default %default] or an integer or two integers seperated by a comma (e.g. 2,4)", metavar="character"),
    
    ###QC specific inputs
    make_option( "--qcdir",type="character", default='fastqc_results', 
                 help="QC argument: output directory [default %default]", metavar="character"),
    
    make_option( "--metadata",type="character", default='data/librairies_1_a_70_et_RIN.xlsx', 
                 help="QC argument: output directory [default %default]", metavar="character"),
    
    make_option( "--fastqc",type="character", default='fastqc', 
                 help="QC argument: where is fastqc installed [default %default]", metavar="character"),
    
    ###kallisto specific inputs
    make_option( "--reftranscriptome",type="character", default='data/reference_transcriptome/chr1.fasta.gz', 
                 help="reference transcriptome  [default %default]", metavar="character"),
    
    make_option( "--quantdir",type="character", default='out/kallisto', 
                 help="quantifications output directory [default %default]", metavar="character"),

    make_option( "--txdbdir",type="character", default='data/reference_transcriptome/gencode.v43.annotation.gtf', 
                 help="txdb directory [default %default]", metavar="character"),
    
    make_option( "--txdbfile",type="character", default='data/reference_transcriptome/gencode.v43.annotation.sqlite', 
                 help="txdb SQL directory [default %default]", metavar="character"),

    ###fusion specific inputs
    make_option( "--fusdir",type="character", default='out/fusion', 
                 help="fusion directory [default %default]", metavar="character"),
    
    make_option( "--ctatdir",type="character", default='data/ctat', 
                 help="CTAT Genome library directory [default %default]", metavar="character"),  
    
    make_option( "--ctatlib",type="character", default='https://data.broadinstitute.org/Trinity/CTAT_RESOURCE_LIB/GRCh38_gencode_v37_CTAT_lib_Mar012021.plug-n-play.tar.gz', 
                 help="CTAT genome library plug-and-play [default %default]", metavar="character"),
    
    make_option( "--CPU",type="integer", default=10, 
                 help="Nb of CPU to use [default %default]", metavar="character"),
    
    ###cnv specific inputs
    make_option( "--cytobanddir",type="character", default='data/cnv/cytoband/cytoBand.txt',
                 help="Original Cytoband information directory [default %default]", metavar="character"),
    
    make_option( "--centromeredir",type="character", default='data/cnv/centromere/cytoBand.txt', 
                 help="Original Centromere information directory [default %default]", metavar="character"),
    
    make_option( "--cnvdir",type="character", default='out/cnv', 
                 help="cnv output directory [default %default]", metavar="character"),
    
    make_option( "--countsdir",type="character", default='out/gene_counts',
                 help="Count files directory [default %default]", metavar="character"),
    
    make_option( "--cytodir",type="character", default='data/cnv/cytoband/cytoband.txt',
                 help="Final Cytoband information directory [default %default]", metavar="character"),
    
    make_option( "--centrodir",type="character", default='data/cnv/centromere/centromere.txt', 
                 help="Final Centromere information directory [default %default]", metavar="character")
    
) 

#
parser <- OptionParser(usage = "%prog [options] QC/counts", option_list=option_list)
arguments <- parse_args(parser,positional_arguments = 1)

if(arguments$args == 'QC' | arguments$args == 'counts' | arguments$args == 'kallisto' | arguments$args == 'mini' | arguments$args == 'fusion' | arguments$args == 'cnv') {
    sprintf("Running command ( %s )", arguments$args)
} else {
    stop(sprintf("Specified command ( %s ) does not exist", arguments$args))
}



#running counts_rnaseq wrapper
if(arguments$args == 'counts') {
    counts_rnaseq(fq.dir = arguments$options$fqdir,
                  trim.dir= arguments$options$trimdir,
                  genomedir = arguments$options$genomeindexdir,
                  annotation.gtf = arguments$options$annotationgtf,
                  genomefasta = arguments$options$genomefasta,
                  out.dir = arguments$options$outdir,
                  cutadapt=arguments$options$cutadapt,
                  threads= arguments$options$threads,
                  nbfiles= arguments$options$nbfiles
    )
}

#running QC notebook
if(arguments$args == 'QC') {
    dir.create(arguments$options$outdir)
    library(rmarkdown)
    rmarkdown::render('./Rmarkdown/fastqc_reports.Rmd',params = list(fq.dir = arguments$options$fqdir,
                                                                     qc.dir = file.path(arguments$options$outdir,arguments$options$qcdir),
                                                                     threads =  arguments$options$threads,
                                                                     metadata = arguments$options$metadata,
                                                                     fastqc.path = arguments$options$fastqc,
                                                                     adapters.dir = NULL),
                      output_dir = arguments$options$outdir
    )
}

#running kallisto_rnaseq
if (arguments$args == 'kallisto') {
    kallisto_rnaseq(trim.dir = arguments$options$trimdir,
                    quant.dir = arguments$options$quantdir,
                    ref.transcriptome = arguments$options$reftranscriptome,
                    threads = arguments$options$threads,
                    txdb.dir = arguments$options$txdbdir,
                    txdb.file = arguments$options$txdbfile)
}
#running fusion_rnaseq
if (arguments$args == 'fusion') {
    fusion_rnaseq(ctat.dir = arguments$options$ctatdir,
                  ctat.lib = arguments$options$ctatlib,
                  trim.dir = arguments$options$trimdir,
                  fus.dir = arguments$options$fusdir,
                  CPU = arguments$options$CPU)
}

#running cnv_rnaseq
if (arguments$args == 'cnv') {
    cnv_rnaseq(cyto.dir = arguments$options$cytodir,
               cytoband.dir = arguments$options$cytobanddir,
               centro.dir = arguments$options$centrodir,
               centromere.dir = arguments$options$centromeredir,
               counts.dir = arguments$options$countsdir,
               cnv.dir = arguments$options$cnvdir)    
}



#running test alignments
if (arguments$args == 'mini') {
    counts_mini_rnaseq(fq.dir = arguments$options$fqdir,
                    trim.dir= arguments$options$trimdir,
                    genomedir = arguments$options$genomeindexdir,
                    annotation.gtf = arguments$options$annotationgtf,
                    genomefasta = arguments$options$genomefasta,
                    out.dir = arguments$options$outdir,
                    threads= arguments$options$threads,
                    nbfiles= arguments$options$nbfiles,
                    head = 4000000)    
}
