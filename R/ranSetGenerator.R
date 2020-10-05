#' @title Random data of viral vector integration site
#' 
#' @description
#' Generate viral vector integration site random data
#' 
#' @usage 
#' ranSetGenerator(organism = 'GRCh37', randomSize = if(doRandom){10000}else{NULL}, outPath = getwd(),
#'               outFileName = paste0('RIPAT', round(unclass(Sys.time()))))
#' 
#' @param organism a single character. This function can run by two versions of organisms
#'                 such as GRCh37, GRCh38 (Human). Default is 'GRCh37'.
#' @param randomSize an integer vector. A random set size. Default is 10000.
#' @param outPath an string vector. Output file is saved in this path. Default value is R home directory.
#' @param outFileName a character vector. Attached ID to the result file name.
#' 
#' @return Return a list of integration site positions on chromosome as an R object and text file.
#'         
#' @examples 
#' ranSetGenerator(organism = 'GRCh37')
#' 
#' @export
ranSetGenerator = function(organism = 'GRCh37', randomSize = if(doRandom){10000}else{NULL}, outPath = getwd(), outFileName = paste0('RIPAT', round(unclass(Sys.time())))){
  message('----- Generate integration sites. (Time : ', date(), ')')
  message('- Validate options')
  if(length(which(c('GRCh37', 'GRCh38') %in% organism)) == 0){
    stop("[ERROR] Please use GRCh37/GRCh38 data.\n----- This process is halted. (Time : ", date(), ")\n")
  }
  if(stringr::str_ends(outPath, pattern = '/')){outPath = stringr::word(outPath, start = 1, end = nchar(outPath), sep = '')}
  message('- OK!')
  message('- Generate random data')
  ch_size = readRDS(file = system.file("extdata", paste0(organism, '_chrom.rds'), package = "RIPAT"))
  ch_size_num = as.numeric(ch_size$length)
  random_num = sort(sample.int(sum(ch_size_num), size = randomSize, replace = FALSE, prob = NULL), decreasing = FALSE)
  flag_num = c(0, cumsum(ch_size_num)[1:23])
  random_set = data.frame(c(1:randomSize), do.call('rbind', lapply(random_num, function(a){
    tmp_chr = paste0('chr', as.character(ch_size$chrom)[max(which(a - flag_num > 0))])
    tmp_pos = a-flag_num[max(which(a - flag_num > 0))]
    return(c(tmp_chr, tmp_pos))
  })), stringsAsFactors = FALSE); names(random_set) = c('Random', 'Random_chr', 'Random_pos')
  ran_tab = random_set[,c(2,3,3)]
  names(ran_tab) = c('chr', 'start', 'end')
  utils::write.table(ran_tab, file = paste0(outPath, '/', outFileName, '_Random_set_', organism, '_gene.txt'), quote = FALSE, append = FALSE, sep = '\t', na = '', row.names = FALSE, col.names = TRUE)
  gr_random = GenomicRanges::makeGRangesFromDataFrame(ran_tab)
  message('- OK!')
  message(' ----- Finish. (Time : ', date(), ')')
  return(gr_random)
}