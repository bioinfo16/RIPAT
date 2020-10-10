#' @title Annotate integration sites by clinical variants.
#' 
#' @description
#' Annotate vector integration sites by clinical variant data.
#' 
#' @usage 
#' annoByVar(hits, ran_hits = NULL,
#'           mapTool = 'blast',
#'           organism = 'GRCh37',
#'           interval = 5000,
#'           range = c(-20000, 20000),
#'           outPath = getwd(),
#'           outFileName = paste0('RIPAT', round(unclass(Sys.time()))))
#' 
#' @param hits GR object. This object made by \code{\link{makeExpSet}} function.
#' @param ran_hits GR object or list. This object is output of \code{\link{makeRanSet}} function.
#' @param mapTool Character. Function uses two types of object\cr from BLAST and BLAT.
#'                Default is 'blast'. If you want to use BLAT result, use 'blat'.
#' @param organism Character. This function can run by two versions of organisms\cr
#'                 such as GRCh37, GRCh38 (Human). Default is 'GRCh37'.
#' @param interval Integer. This number means interval number\cr for
#'                 distribution analysis. Default is 5000.
#' @param range Integer array. The range of highlight region for analysis.\cr
#'              Default range is \code{c(-20000, 20000)}.
#' @param outPath String. Plots are saved in this path. \cr Default value is R home directory.
#' @param outFileName Character. Attached ID to the result file name.
#' 
#' @return Return a result list that is made up of insertion and distribution result tables
#'         and GenomicRange object of clinical variant data.
#'         
#' @examples
#' data(blast_obj); data(var_exam_db)
#' saveRDS(var_exam_db,
#'         paste0(system.file("extdata", package = 'RIPAT'),
#'         '/GRCh37_clinvar.rds'))
#'
#' blast_clivar = annoByVar(hits = blast_obj, ran_hits = NULL,
#'                          outFileName = 'blast_res')
#' 
#' @export
annoByVar = function(hits, ran_hits = NULL, mapTool = 'blast', organism = 'GRCh37', interval = 5000, range = c(-20000, 20000), outPath = getwd(), outFileName = paste0('RIPAT', round(unclass(Sys.time())))){
  message('----- Annotate integration sites. (Time : ', date(), ')')
  message('- Validate options')
  if(length(which(c('blast', 'blat') %in% mapTool)) == 0){stop("[ERROR] Please confirm the alignment tool name.\n----- This process is halted. (Time : ", date(), ")\n")}
  if(length(which(c('GRCh37', 'GRCh38') %in% organism)) == 0){stop("[ERROR] Please use GRCh37/GRCh38 data.\n----- This process is halted. (Time : ", date(), ")\n")}
  if(stringr::str_ends(outPath, pattern = '/')){outPath = stringr::word(outPath, start = 1, end = nchar(outPath), sep = '')}
  if(range[1] - range[2] >= 0){
    stop("[ERROR] Please check distribution range.\n----- This process is halted. (Time : ", date(), ")\n")
  } else {ranges = seq(from = range[1], to = range[2], by = interval)}
  message('- OK!')
  message('- Load the dataset.')
  dbPath = system.file("extdata", package = "RIPAT")
  dbFile = paste0('/', organism, '_clinvar.rds')
  dataTable = data.frame(readRDS(paste0(dbPath, dbFile)), stringsAsFactors = FALSE)
  dataTable = dataTable[,c(19:21,2,3,5,14,15,7,25,26)]
  colnames(dataTable) = c('chrom', 'start', 'end', 'type', 'name', 'symbol', 'phenotypeList', 'origin', 'clinicalSignificance', 'reviewstatus', 'numberSubmitters')
  gr_clin = GenomicRanges::makeGRangesFromDataFrame(dataTable, keep.extra.columns = TRUE, ignore.strand = TRUE)
  message('- OK!')
  message('- Find integration sites located in variants.')
  inside_clin_only = as.data.frame(GenomicRanges::findOverlaps(hits[[1]], gr_clin, type = 'any', ignore.strand = TRUE), stringsAsFactors = FALSE)
  only_res_que = data.frame(hits[[1]][inside_clin_only$queryHits,], stringsAsFactors = FALSE)[,c(1:3,6:8)]
  only_res_sub = dataTable[inside_clin_only$subjectHits,]
  inside_tab_only = unique(cbind(only_res_que, only_res_sub)[,c(4, 1:3, 7:17)])
  colnames(only_res_que) = c('IR_chr', 'IR_start', 'IR_end', 'IR_ID', 'IR_identity', 'IR_align_length')
  names(inside_tab_only) = c('q_name', 'q_chr', 'q_start', 'q_end', 'chr', 'start', 'end', 'type', 'name', 'symbol', 'phenotypeList', 'origin', 'clinicalSignificance', 'reviewstatus', 'numberSubmitters')
  if(length(hits[[2]]) != 0){
    inside_clin_dup = as.data.frame(GenomicRanges::findOverlaps(hits[[2]], gr_clin, type = 'any', ignore.strand = TRUE), stringsAsFactors = FALSE)
    dup_res_que = data.frame(hits[[2]][inside_clin_dup$queryHits,], stringsAsFactors = FALSE)[,c(1:3,6:8)]
    dup_res_sub = dataTable[inside_clin_dup$subjectHits,]
    inside_tab_dup = unique(cbind(dup_res_que, dup_res_sub)[,c(4, 1:3, 7:17)])
    colnames(dup_res_que) = c('IR_chr', 'IR_start', 'IR_end', 'IR_ID', 'IR_identity', 'IR_align_length')
    names(inside_tab_dup) = c('q_name', 'q_chr', 'q_start', 'q_end', 'chr', 'start', 'end', 'type', 'name', 'symbol', 'phenotypeList', 'origin', 'clinicalSignificance', 'reviewstatus', 'numberSubmitters')
  } else {inside_tab_dup = NULL}
  message('- OK!')
  message('- Calculate distance.')
  only_hits_tab = data.frame(hits[[1]], stringsAsFactors = FALSE)
  dist_only = lapply(c(1:nrow(only_hits_tab)), function(a){
    x = only_hits_tab$start[a]; y = only_hits_tab$query[a]; z = only_hits_tab$seqnames[a]
    cal1 = dataTable$start-x; cal2 = dataTable$end-x
    cal1_i = intersect(which(cal1 <= abs(range[1])), which(cal1 > 0)); cal2_i = intersect(which(abs(cal2) <= range[2]), which(cal2 < 0))
    dat = data.frame(dataTable[c(cal1_i, cal2_i),], dist = -c(cal1[cal1_i], cal2[cal2_i]))
    dat = unique(data.frame('query' = rep(y, nrow(dat)), dat))
    dat = dat[which(dat$chr == z),]
    dat = dat[order(abs(dat$dist), decreasing = FALSE),][1,]
    return(dat)
  })
  hist_only_v = hist(unlist(lapply(dist_only, function(x){x$dist/1000})), breaks = ranges/1000, plot = FALSE)
  if(length(hits[[2]]) != 0){
    dup_hits_tab = data.frame(hits[[2]], stringsAsFactors = FALSE)
    dist_dup = lapply(c(1:nrow(dup_hits_tab)), function(a){
      x = dup_hits_tab$start[a]; y = dup_hits_tab$query[a]; z = dup_hits_tab$seqnames[a]
      cal1 = dataTable$start-x; cal2 = dataTable$end-x
      cal1_i = intersect(which(cal1 <= abs(range[1])), which(cal1 > 0)); cal2_i = intersect(which(abs(cal2) <= range[2]), which(cal2 < 0))
      dat = data.frame(dataTable[c(cal1_i, cal2_i),], dist = -c(cal1[cal1_i], cal2[cal2_i]))
      dat = unique(data.frame('query' = rep(y, nrow(dat)), dat))
      dat = dat[which(dat$chr == z),]
      dat = dat[order(abs(dat$dist), decreasing = FALSE),][1,]
      return(dat)
    })
    hist_dup_v = hist(unlist(lapply(dist_dup, function(x){x$dist/1000})), breaks = ranges/1000, plot = FALSE)
  }
  message('- OK!')
  if(!is.null(ran_hits)){
    message('- Do random set analysis.')
    randomSize = length(ran_hits); gr_random = ran_hits
    random_set = data.frame(c(1:randomSize), data.frame(ran_hits), stringsAsFactors = FALSE)[,c(1:3)]; names(random_set) = c('Random', 'Random_chr', 'Random_pos')
    inside_cl_ran = as.data.frame(GenomicRanges::findOverlaps(gr_random, gr_clin, type = 'any', ignore.strand = TRUE), stringsAsFactors = FALSE)
    a = as.data.frame(gr_random[inside_cl_ran$queryHits,], stringsAsFactors = FALSE)
    b = dataTable[inside_cl_ran$subjectHits,]
    inside_ran_tab = unique(cbind(a,b)[,c(1:3, 4:14)])
    names(inside_ran_tab) = c('q_chr', 'q_start', 'q_end', 'chr', 'start', 'end', 'type', 'name', 'symbol', 'phenotypeList', 'origin', 'clinicalSignificance', 'reviewstatus', 'numberSubmitters')
    cl_dist_ran = lapply(c(1:nrow(ran_tab)), function(a){
      x = as.numeric(random_set$Random_pos[a]); y = random_set$Random[a]; z = random_set$Random_chr[a]
      cal1 = dataTable$start-x; cal2 = dataTable$end-x
      cal1_i = intersect(which(cal1 <= abs(range[1])), which(cal1 > 0)); cal2_i = intersect(which(abs(cal2) <= range[2]), which(cal2 < 0))
      dat = data.frame(dataTable[c(cal1_i, cal2_i),], dist = -c(cal1[cal1_i], cal2[cal2_i]))
      dat = unique(data.frame('query' = rep(y, nrow(dat)), dat))
      dat = dat[which(dat$chr == z),]
      dat = dat[order(abs(dat$dist), decreasing = FALSE),][1,]
      return(dat)
    })
    all_dist_ran = unlist(lapply(cl_dist_ran, function(x){x$dist}))
    hist_obj_ran = hist(all_dist_ran, plot = FALSE, breaks = ranges)
    message('- OK!')
  } else {message('[WARN] Skip random set analysis.')}
  message('- Draw histograms.')
  all_dist_only = unlist(lapply(dist_only, function(x){x$dist}))
  if(length(hits[[2]]) != 0){
    all_dist_dup = unlist(lapply(dist_dup, function(x){x$dist}))
  } else {all_dist_dup = NULL}
  hist_obj = hist(all_dist_only, plot = FALSE, breaks = ranges)
  inside_tab = inside_tab_only
  cl_dist = list('Decided' = dist_only)
  if(!is.null(ran_hits)){
    count_site = hist_obj$counts; count_site_ran = hist_obj_ran$counts
    count_all = nrow(only_hits_tab)
    count_data = data.frame('Range' = factor(rep(ranges[ranges != 0]/1000, 2), levels = ranges[ranges != 0]/1000), 'Group' = c(rep('Observed', length(count_site)), rep('Random', length(count_site_ran))), 'Count' = c(count_site, count_site_ran), 'Freq' = c(count_site/count_all, count_site_ran/randomSize))
  } else {
    count_site = hist_obj$counts;
    count_all = nrow(only_hits_tab)
    count_data = data.frame('Range' = factor(ranges[ranges != 0]/1000, levels = ranges[ranges != 0]/1000), 'Group' = rep('Observed', length(count_site)), 'Count' = count_site, 'Freq' = count_site/count_all)
  }
  grDevices::png(paste0(outPath, '/', outFileName, '_distribution_var_', organism, '.png'), width = 1200, height = 750)
  cl_plot = ggplot2::ggplot(count_data) + ggplot2::geom_bar(ggplot2::aes(x = Range, y = Freq, fill = Group), stat = "identity", position = "dodge", width = 0.5) +
    ggplot2::lims(y = c(0, max(count_data$Freq)*1.5)) + ggplot2::ggtitle(label = "Random distribution (Clinical variant)") +
    ggplot2::xlab('Intervals (Kbs)') + ggplot2::ylab("Ratio of Integration Events") + ggplot2::scale_fill_manual(values = c('mediumspringgreen', 'mediumpurple')) +
    ggplot2::theme(panel.background = ggplot2::element_rect(fill="white", colour = "white"),
                   panel.grid.major = ggplot2::element_line(size = 0.5, linetype = 'dotted', colour = 'black'),
                   axis.line = ggplot2::element_line(colour = "darkgrey"), legend.title = ggplot2::element_blank(),
                   legend.key.size = ggplot2::unit(0.7, "cm"), plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 20),
                   legend.text = ggplot2::element_text(size = 18), axis.text = ggplot2::element_text(size = 17),
                   axis.title = ggplot2::element_text(size = 18))
  print(cl_plot)
  grDevices::dev.off()
  message('- OK!')
  result_list = list(inside_tab, cl_dist, count_data, gr_clin, organism)
  names(result_list) = c('Variant_inside', 'Variant_distribution', 'Variant_plot_data', 'Variant_data', 'Target_ver')
  if(!is.null(ran_hits)){
    result_list = c(result_list, list(inside_ran_tab, cl_dist_ran))
    names(result_list)[6:7] = c('Variant_inside_ran', 'Random_distribution')
  }
  result_list = c(result_list, list('Variant_hist_plot' = cl_plot))
  message('----- Finish. (Time : ', date(), ')')
  return(result_list)
}
