#' @title Make data files for RIPAT
#' 
#' @description \preformatted{
#' RIPAT use specific format of data files for retroviral vector 
#' integration pattern analysis. This function made these data files from
#' well-known genomic databases easily. The data file location is extdata folder of RIPAT.
#' }
#'              
#' @usage 
#' makeData(organism = 'GRCh37', dataType = 'gene')
#' 
#' @param organism a single character. 2 versions of organisms such as GRCh37, GRCh38 (Human). Default is 'GRCh37'.
#' @param dataType a single character. Data type what user needs (gene, cpg, repeat and variant). Default is 'gene'.
#' 
#' @examples 
#'
#' makeData(organism = 'GRCh37')
#' 
#' @export
makeData = function(organism = 'GRCh37', dataType = 'gene'){
  message('----- Make data files for RIPAT. (Time : ', date(), ')')
  message('- Validate options')
  if(length(which(c('GRCh37', 'GRCh38') %in% organism)) == 0){
    stop("[ERROR] Please use GRCh37/GRCh38 data.\n----- This process is halted. (Time : ", date(), ")\n")
  }
  message('- OK!')
  outPath = system.file("extdata", package = 'RIPAT')
  cat('+ The data file path : ', outPath, '\n')
  if(organism == 'GRCh37'){otherkey = 'hg19'}else if(organism == 'GRCh38'){otherkey = 'hg38'}
  if(dataType == 'gene'){
    message('- Load ensembl data (Gene, TSS)')
    if(organism == 'GRCh37'){
      usedMart = biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL", host = "grch37.ensembl.org", path = "/biomart/martservice")
      del_i1 = c(2737,5285,5633,5642,6579,6845,6292,6843,6293,8668,8930,9349,6262,10267,5359,10376,6235,11311,13654,16282,17113,17404,6845,17514,17588,17536,17614,19303,19310,17872,19305,17871,19311,4783,20265,22151,6591,25877,27762,6590,22267,27761,6589,22266,25875,22265,25874,27759,31248,24484)
      del_i2 = c(1382,1995,2126,2268,2413,3930,4596,4742,4978,5314,5740,6595,6773,7155,8937,9632,9866,12069,12659,13219,14343,15367,15720,16541,17756,18890,19718,19785,19908,21941,21942,21943,25336,30823)
    } else if(organism == 'GRCh38'){
      usedMart = biomaRt::useMart("ensembl")
      del_i1 = c(1180,1372,2284,5620,10266,11649,12247,12797,15042,17070,22217,24767,25275,25436,3135,35456,36936,37072,37759,37932,38006)
      del_i2 = c(7521,18038,18467,23294)
    }
    usedDataset = biomaRt::useDataset(dataset = "hsapiens_gene_ensembl", mart = usedMart)
    gtab_raw = biomaRt::getBM(attributes = c("hgnc_symbol" ,"ensembl_gene_id", "gene_biotype", "chromosome_name", "start_position", "end_position", "strand", "percentage_gene_gc_content", "description"), mart = usedDataset)
    gtab = subset(gtab_raw, gtab_raw$hgnc_symbol != '')
    gtab = subset(gtab, gtab$chromosome_name != 'MT')
    gtab = subset(gtab, !stringr::str_detect(gtab$chromosome_name, '_'))
    gtab = subset(gtab, !stringr::str_detect(gtab$chromosome_name, 'GL'))
    gtab = subset(gtab, !stringr::str_detect(gtab$chromosome_name, 'KI'))
    row.names(gtab) = c(1:nrow(gtab))
    gtab1 = gtab[-unique(del_i1),]; row.names(gtab1) = c(1:nrow(gtab1))
    gtab2 = gtab1[-unique(del_i2),]; gtab2$chromosome_name = paste0('chr', gtab2$chromosome_name)
    strand1 = gtab2$strand; strand1[which(strand1 == -1)] = '-'; strand1[which(strand1 == 1)] = "+"
    gtab2$strand = strand1
    ttab_raw = biomaRt::getBM(attributes = c("hgnc_symbol" ,"ensembl_gene_id", "ensembl_transcript_id", "transcript_biotype", "chromosome_name", "transcription_start_site", "transcription_start_site", "strand", "percentage_gene_gc_content", "description"), mart = usedDataset)
    ttab = subset(ttab_raw, ttab_raw$hgnc_symbol != '')
    ttab = subset(ttab, ttab$chromosome_name != 'MT')
    ttab = subset(ttab, !stringr::str_detect(ttab$chromosome_name, '_'))
    ttab = subset(ttab, !stringr::str_detect(ttab$chromosome_name, 'GL'))
    ttab = subset(ttab, !stringr::str_detect(ttab$chromosome_name, 'KI'))
    ttab1 = ttab[which(paste0(ttab$hgnc_symbol, '--', ttab$ensembl_gene_id) %in% paste0(gtab2$hgnc_symbol, '--', gtab2$ensembl_gene_id)),]
    ttab2 = ttab1[,-2]; names(ttab2)[5:6] = c('start_position', 'end_position')
    ttab2$chromosome_name = paste0('chr', ttab2$chromosome_name)
    strand2 = ttab2$strand; strand2[which(strand2 == -1)] = '-'; strand2[which(strand2 == 1)] = '+'
    ttab2$strand = strand2
    message('- Save gene and TSS data')
    saveRDS(gtab2, file = paste0(outPath, '/', organism, '_gene.rds'))
    saveRDS(ttab2, file = paste0(outPath, '/', organism, '_TSS.rds'))
    message('- OK!')
  } else if(length(which(dataType == c('cpg', 'repeat'))) != 0){
    message('- Load UCSC data (Repeat, Microsatellite, CpG island)')
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
    vtab_raw = read.delim(gzfile(paste0(outPath, '/variant_summary.txt.gz')), header = TRUE, stringsAsFactors = FALSE)
    vtab = subset(vtab_raw, vtab_raw$Assembly == organism)
    vtab = subset(vtab, vtab$Chromosome != 'MT')
    vtab1 = subset(vtab, vtab$ReviewStatus %in% c('reviewed by expert panel', 'criteria provided, multiple submitters, no conflicts'))
    vtab1$Chromosome = paste0('chr', vtab1$Chromosome)
    message('- Save clinical variant data')
    saveRDS(vtab1, file = paste0(outPath, '/', organism, '_clinvar.rds'))
    message('- OK!')
  } else {
    stop("[ERROR] Please enter gene, cpg, repeat and variant.\n----- This process is halted. (Time : ", date(), ")\n")
  }
   message('----- Finish. (Time : ', date(), ')\n')
}


