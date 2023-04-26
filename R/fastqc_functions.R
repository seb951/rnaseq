#fastqc function


#--------------
###hack for windows
#--------------
fastqc_windows = function (fq.dir = getwd(), qc.dir = NULL, threads = 8, fastqc.path = "", adapters.dir = NULL) 
{
  fq.dir = sub('C:','/mnt/c',fq.dir)
  qc.dir = sub('C:','/mnt/c',qc.dir)
  adapters.dir = sub('C:','/mnt/c',adapters.dir)
  fastqc.path = sub('C:','/mnt/c',fastqc.path)
  
  if (file.exists(qc.dir)==F) {dir.create(qc.dir)}
  
  cmd <- paste0('wsl.exe ',fastqc.path, " ", fq.dir, "/*R1.fastq.gz ",fq.dir,"/*R2.fastq.gz --threads ", 
                threads, " --outdir ", qc.dir, ifelse(is.null(adapters.dir),"",paste0(' --adapters ',adapters.dir)))
  system(cmd)
}

#--------------
####regular function
#--------------
fastqc2 = function (fq.dir = getwd(), qc.dir = NULL, threads = 8, fastqc.path = "", adapters.dir = NULL) 
{
  if (file.exists(qc.dir)==F) {dir.create(qc.dir)}
  
  cmd <- paste0(fastqc.path, " ", fq.dir, "/*R1.fastq.gz ",fq.dir,"/*R2.fastq.gz --threads ", 
                threads, " --outdir ", qc.dir, ifelse(is.null(adapters.dir),"",paste0(' --adapters ',adapters.dir)))
  system(cmd)
}


#--------------
###wrapper
#--------------
fastqc_wrapper = function(fq.dir = 'path',
                           qc.dir = 'path', 
                           threads = 8,
                           fastqc.path = 'path',
                           metadata.dir = NULL,
                           adapters.dir = NULL){
  
  if(length(list.files(qc.dir))==0) {
    if(Sys.info()['sysname'] != 'Windows') {
      fastqc2(qc.dir = qc.dir ,fq.dir = fq.dir,threads = threads,fastqc.path = fastqc.path, adapters.dir = adapters.dir)}
    
    if(Sys.info()['sysname'] == 'Windows') {
      message('waning: Running a modified version of fastqc() for Windows');
      fastqc_windows(qc.dir = qc.dir,
                     fq.dir = fq.dir,
                     threads = threads,
                     fastqc.path = fastqc.path,
                     adapters.dir = adapters.dir)}
  } else {message(paste0(Sys.time(),' --- Directory ',qc.dir, ' is not empty. Will assume that fastqc has already ran.'))}
  
  #zip the whole thing
  zip_cmd = paste0('zip -j ',file.path(qc.dir,"../fastqc_individual_reports.zip "),file.path(qc.dir,"*fastqc.zip"))
  if(file.exists(file.path(fq.dir,"../fastqc_individual_reports.zip"))==F) system(zip_cmd)
  
  message(paste0(Sys.time(),
                 ' --- Zip archive is stored here: ',
                 file.path(fq.dir,"../fastqc_individual_reports.zip")))
  
  #aggregate results
  qc = suppressMessages(qc_aggregate(qc.dir = qc.dir, progressbar = F))
  
  # Generates a summary of qc_aggregate
  qc_summary = summary(qc)
  
  # General statistics of fastqc reports.
  qc_stats = qc_stats(qc)
  qc_stats$ends = sapply(strsplit(qc_stats$sample,split = '_'),tail, 1)
  qc_stats$bio_sample = gsub('_R[12]','',qc_stats$sample)
  
  #fix (shorten) names
  if(length(grep('i5.',qc_stats$bio_sample[1]))==1) {
    qc_stats$bio_sample = sapply(strsplit(qc_stats$bio_sample,split ='i5.',fixed = T),"[[", 2)
  }
  
  #qc_read_collection
  qc.files <- list.files(qc.dir,pattern = '.zip', full.names = TRUE)
  qc.names = gsub("_fastqc.zip","",list.files(qc.dir,pattern = '.zip',full.names = F))
  
  qc_read_collection <- suppressMessages(qc_read_collection(qc.files,sample_names = qc.names,modules = c("all")))
  
  mean_quality_perfile = qc_read_collection$per_base_sequence_quality %>% select(c('sample','Mean')) %>% group_by(sample) %>% summarise(mean_quality= mean(Mean))
  
  mean_quality_perfile$mean_quality = signif(mean_quality_perfile$mean_quality,4)
  
  qc_stats = merge(qc_stats,mean_quality_perfile)
 
  # metadata
  if(!is.null(metadata.dir) & file.exists(metadata.dir)){
    metadata = suppressMessages(read_xlsx(metadata.dir))
    rin = metadata %>% select(c(Nom...1,RIN))
    colnames(rin)[1] = 'bio_sample'
    qc_stats =  merge(qc_stats,rin)
  }
  
  qc_stats$type = ifelse(regexpr('umeur',qc_stats$sample)>0,'Tumeur','Sain')
  qc_stats$tot.seq =  as.numeric(qc_stats$tot.seq) / 1000000
  
  message(paste0(Sys.time(),' --- Done preparing QC data'))
  
  return(list(qc_stats,qc_summary,qc_read_collection))
}