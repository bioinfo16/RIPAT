#' @title Make data files for RIPAT.
#' 
#' @description
#' Download datafiles for running RIPAT.
#'              
#' @usage 
#' makeData(organism = 'GRCh37', dataType = 'cpg')
#' 
#' @param organism a single character. Two versions of organism such as GRCh37, GRCh38 (Human).\cr
#'                 Default is 'GRCh37'.
#' @param dataType a single character. Data type what user needs (cpg, repeat and variant).\cr
#'                 Default is 'cpg'.
#' 
#' @examples 
#' makeData(organism = 'GRCh37')
#' 
#' @export
makeData = function(organism = 'GRCh37', dataType = 'cpg'){
  message('----- Make data files for RIPAT. (Time : ', date(), ')')
  message('- Validate options')
  if(length(which(c('GRCh37', 'GRCh38') %in% organism)) == 0){
    stop("[ERROR] Please use GRCh37/GRCh38 data.\n----- This process is halted. (Time : ", date(), ")\n")
  }
  message('- OK!')
  outPath = system.file("extdata", package = 'RIPAT')
  cat('+ The data file path : ', outPath, '\n')
  if(organism == 'GRCh37'){otherkey = 'hg19'}else if(organism == 'GRCh38'){otherkey = 'hg38'}
  if(length(which(dataType == c('cpg', 'repeat'))) != 0){
    message('- Load UCSC data')
    UCSCSession = rtracklayer::browserSession("UCSC")
    rtracklayer::genome(UCSCSession) <- otherkey
    chr_info = readRDS(file = system.file("extdata", paste0(organism, '_chrom.rds'), package = 'RIPAT'))
    gr_chr = GenomicRanges::GRanges(seqnames = paste0('chr', chr_info$chrom),
                                    ranges = IRanges::IRanges(start = 1, end = chr_info$length))
    if(dataType == 'cpg'){
      ctab_list = lapply(c(1:length(gr_chr)), function(a){
        rtracklayer::getTable(rtracklayer::ucscTableQuery(UCSCSession, track = "cpgIslandExt", range = gr_chr[a], table = "cpgIslandExt"))
      })
      ctab = data.frame(do.call(rbind, ctab_list), stringsAsFactors = FALSE)[,-1]
      names(ctab) = c("chrom", "start", "end", "name", "length", "cpgNum", "gcNum", "perCpg", "perGc", "obsExp")
      message('- Save CpG island data')
      saveRDS(ctab, file = paste0(outPath, '/', organism, '_cpg.rds'))
    } else if(dataType == 'repeat'){
      rtab_list = lapply(c(1:length(gr_chr)), function(a){
        rtracklayer::getTable(rtracklayer::ucscTableQuery(UCSCSession, track = "rmsk", range = gr_chr[a], table = "rmsk"))
      })
      rtab = data.frame(do.call(rbind, rtab_list), stringsAsFactors = FALSE)[,c(6:8,10:13)]
      names(rtab) = c("genoName", "genoStart", "genoEnd", "strand", "repName", "repClass", "repFamily")
      rtab = subset(rtab, rtab$repClass != 'Simple_repeat')
      rtab = subset(rtab, rtab$repClass != 'Unknown')
      rtab = subset(rtab, !stringr::str_detect(rtab$repClass, "[?]"))
      rtab = subset(rtab, !stringr::str_detect(rtab$repFamily, "[?]"))
      rtab = subset(rtab, !stringr::str_detect(rtab$repName, '[?]'))
      mtab_list = lapply(c(1:length(gr_chr)), function(a){
        rtracklayer::getTable(rtracklayer::ucscTableQuery(UCSCSession, track = "microsat", range = gr_chr[a], table = "microsat"))
      })
      mtab = data.frame(do.call(rbind, mtab_list), stringsAsFactors = FALSE)[,-1]
      names(mtab) = c("chrom", "chromStart", "chromEnd", "name")
      message('- Save repeat and microsatellite data')
      saveRDS(rtab, file = paste0(outPath, '/', organism, '_repeat.rds'))
      saveRDS(mtab, file = paste0(outPath, '/', organism, '_microsat.rds'))
    }
    message('- OK!')
  } else if(dataType == 'variant'){
    message('- Load NCBI Clinvar data')
    utils::download.file(url = "http://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz",
                         destfile = paste0(outPath, '/variant_summary.txt.gz'))
    vtab_raw = utils::read.delim(gzfile(paste0(outPath, '/variant_summary.txt.gz')), header = TRUE, stringsAsFactors = FALSE)
    vtab = subset(vtab_raw, vtab_raw$Assembly == organism)
    vtab = subset(vtab, vtab$Chromosome != 'MT')
    vtab1 = subset(vtab, vtab$ReviewStatus %in% c('reviewed by expert panel', 'criteria provided, multiple submitters, no conflicts'))
    vtab1 = subset(vtab, vtab$NumberSubmitters >= 2)
    vtab1$Chromosome = paste0('chr', vtab1$Chromosome)
    message('- Save clinical variant data')
    saveRDS(vtab1, file = paste0(outPath, '/', organism, '_clinvar.rds'))
    message('- OK!')
  } else {
    stop("[ERROR] Please enter cpg, repeat and variant.\n----- This process is halted. (Time : ", date(), ")\n")
  }
  message('----- Finish. (Time : ', date(), ')\n')
}