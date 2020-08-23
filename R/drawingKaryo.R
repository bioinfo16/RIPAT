#' @title Draw the karyogram plot.
#' 
#' @description
#' Draw a karyogram plot and show integration site.
#' 
#' @usage 
#' drawingKaryo(hits, feature, organism = 'GRCh37',
#'              includeUndecided = FALSE, outPath = getwd(), 
#'              outFileName = paste0('RIPAT', round(unclass(Sys.time()))))
#' 
#' @param hits a GR object. This object made from \code{makeInputObj} function.
#' @param feature a GR object. This object made from annotation function. 
#' @param organism a character vector. This function serves 2 versions of organisms
#'                 such as GRCh37, GRCh38 (Human). Default is 'GRCh37'.
#' @param includeUndecided TRUE or FALSE. If user want to use undecided hits in analysis,
#'                         enter TRUE. Default is FALSE.
#' @param outPath a string vector. Type path to save a plot.
#' @param outFileName a character vector. This value used when saving the idegoram image file.
#' 
#' @return Return the ideogram plot and object.
#' 
#' @examples
#' data(blast_obj)
#' data(blast_gene)
#' drawingKaryo(hits = blast_obj, feature = blast_gene$Gene_data, outFileName = 'blast_res')
#' 
#' @export
drawingKaryo = function(hits, feature, organism = 'GRCh37', includeUndecided = FALSE, outPath = getwd(), outFileName = paste0('RIPAT', round(unclass(Sys.time())))){
  message('----- Draw an ideogram. (Time : ', date(), ')')
  message('- Validate options')
  if(length(which(c('GRCh37', 'GRCh38') %in% organism)) == 0){
    stop("[ERROR] Please use GRCh37/GRCh38 data.\n----- This process is halted. (Time : ", date(), ")\n")
  }
  if(stringr::str_ends(outPath, pattern = '/')){outPath = stringr::word(outPath, start = 1, end = nchar(outPath), sep = '')}
  message('- OK!')
  message('- Load the dataset.')
  chr_length = readRDS(file = system.file("extdata", paste0(organism, '_chrom.rds'), package = "RIPAT"))
  feats = as.data.frame(feature, stringsAsFactors = FALSE)
  feats = feats[,c(1,2,3)]
  feats = data.frame(cbind(feats, 'gieStain' = rep('feats', nrow(feats))), stringsAsFactors = FALSE)
  tmp_len = data.frame(data.frame('chrom' = paste0('chr', chr_length[,1]), 'start' = 1, 'end' = chr_length[,2], stringsAsFactors = FALSE), 'gieStain' = rep('chr', nrow(chr_length)), stringsAsFactors = FALSE)
  names(feats) = names(tmp_len)
  feats = data.frame(rbind(feats, tmp_len), stringsAsFactors = FALSE)
  message("- OK!")
  message("- Draw an ideogram plot.")
  grDevices::png(filename = paste0(outPath, '/Ideogram_', outFileName, '_', organism, '_decided.png'), width = 1200, height = 750)
  ideo = karyoploteR::plotKaryotype(genome = regioneR::toGRanges(tmp_len), cytobands = feats)
  exp_data = karyoploteR::kpPlotMarkers(ideo, hits[[1]], r1 = 0.5, labels = rep('O', length(hits[[1]])), line.color = 'blue', label.color = 'blue', cex = 0.5)
  grDevices::dev.off()
  if(includeUndecided){
    grDevices::png(filename = paste0(outPath, '/Ideogram_', outFileName, '_', organism, '_undecided.png'), width = 1200, height = 750)
    ideo = karyoploteR::plotKaryotype(genome = regioneR::toGRanges(tmp_len), cytobands = feats)
    exp_data2 = karyoploteR::kpPlotMarkers(ideo, hits[[2]], r1 = 0.5, labels = rep('O', length(hits[[2]])), line.color = 'forestgreen', label.color = 'forestgreen', cex = 0.5)
    grDevices::dev.off()
  }
  message('----- Finish. (Time : ', date(), ')')
  return(ideo)
}
