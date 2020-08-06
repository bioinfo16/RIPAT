#' @title Annotate the vector integration site by repeat sequence data.
#' 
#' @description \preformatted{
#' This function uses repeat sequence data to search the feature of integration regions.
#' User can get query sequence inserted in repeat sites and distribution from repeat coordinations.
#' Plus, users can do random distribution analysis by this function.
#' }
#' 
#' @usage 
#' annoByRepeat(hits, mapTool = 'blast', organism = 'GRCh37', interval = 5000, 
#'              range = c(-20000, 20000), doRandom = TRUE, randomSize = if(doRandom){10000}else{NULL}, 
#'              includeUndecided = FALSE, outPath = getwd(),
#'              outFileName = paste0('RIPAT', round(unclass(Sys.time()))))
#' 
#' @param hits a GR object. This object made by \code{makeInputObj} function.
#' @param mapTool a single character. Function serves two types of object such as outputs from BLAST and BLAT.
#'                Default is 'blast'. If you want to use BLAT result, use 'blat'.
#' @param organism a single character. This function can run by 2 versions of organisms such as GRCh37, GRCh38 (Human). Default is 'GRCh37'.
#' @param interval an integer vector. This number means interval number for distribution analysis. Default is 5000.
#' @param range an integer array. It means the range for highlight region of this analysis. Default range is c(-20000, 20000).
#' @param doRandom TRUE or FALSE. If user types TRUE, random set is generated and do random distribution analysis.
#'                 If this value is FALSE, random distribution analysis is not executed. Default is TRUE.
#' @param randomSize an integer vector. A random set size. Default is 10000.
#' @param includeUndecided TRUE or FALSE. If user want to use undecided hits in analysis, enter TRUE.
#'                         Default is FALSE.
#' @param outPath an string vector. Plots are saved in this path. Default value is R home directory.
#' @param outFileName a character vector. Attached ID to the result file name.
#' 
#' @return Return a result list that is made up of insertion and distribution result tables and GenomicRange object of Repeat data.
#' @examples 
#' data(blast_obj)
#' makeData(organism = 'GRCh37', dataType = 'repeat')
#' blast_repeat = annoByRepeat(hits = blast_obj, doRandom = FALSE, outFileName = 'blast_res')
#' 
#'           
#' @export
annoByRepeat = function(hits, mapTool = 'blast', organism = 'GRCh37', interval = 5000, range = c(-20000, 20000), doRandom = TRUE, randomSize = if(doRandom){10000}else{NULL}, includeUndecided = FALSE, outPath = getwd(), outFileName = paste0('RIPAT', round(unclass(Sys.time())))){
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
  dbFile = paste0('/', organism, '_repeat.rds')
  dataTable = data.frame(readRDS(paste0(dbPath, dbFile)), stringsAsFactors = FALSE)
  colnames(dataTable) = c('chr', 'start', 'end', 'strand','repName', 'repClass', 'repFamily')
  dbFile2 = paste0('/', organism, '_microsat.rds')
  dataTable_micro = data.frame(readRDS(paste0(dbPath, dbFile2)), stringsAsFactors = FALSE)  
  colnames(dataTable_micro) = c('chr', 'start', 'end', 'name')
  dataTable$chr = as.character(dataTable$chr); dataTable_micro$chr = as.character(dataTable_micro$chr)
  gr_repeat = GenomicRanges::makeGRangesFromDataFrame(dataTable, keep.extra.columns = TRUE, ignore.strand = FALSE)
  gr_micro = GenomicRanges::makeGRangesFromDataFrame(dataTable_micro, keep.extra.columns = TRUE, ignore.strand = TRUE)
  message('- OK!')
  message('- Find integration sites located in repeats.')
  inside_repeat_only = as.data.frame(GenomicRanges::findOverlaps(hits[[1]], gr_repeat, type = 'any', ignore.strand = TRUE), stringsAsFactors = FALSE)
  inside_micro_only = as.data.frame(GenomicRanges::findOverlaps(hits[[1]], gr_micro, type = 'any', ignore.strand = TRUE), stringsAsFactors = FALSE)
  only_res_que = data.frame(hits[[1]][inside_repeat_only$queryHits,], stringsAsFactors = FALSE)[,c(1:3,6:8)]
  only_res_sub = dataTable[inside_repeat_only$subjectHits,]
  only_res_que_m = data.frame(hits[[1]][inside_micro_only$queryHits,], stringsAsFactors = FALSE)[,c(1:3,6:8)]
  only_res_sub_m = dataTable[inside_micro_only$subjectHits,]
  inside_tab_only = unique(cbind(only_res_que, only_res_sub)[,c(4,1:3,5,6,11,7:10,12,13)])
  inside_tab_only_m = unique(cbind(only_res_que_m, only_res_sub_m)[,c(4,1:3,5,6,11,7:10,12,13)])
  names(inside_tab_only) = c('q_name', 'q_chr', 'q_start', 'q_end', 'identity', 'align_length', 'rep_name', 'chr', 'start', 'end', 'strand', 'repClass', 'repFamily')
  names(inside_tab_only_m) = c('q_name', 'q_chr', 'q_start', 'q_end', 'identity', 'align_length', 'rep_name', 'chr', 'start', 'end', 'strand', 'repClass', 'repFamily')
  if(length(hits[[2]]) != 0){
    inside_repeat_dup = as.data.frame(GenomicRanges::findOverlaps(hits[[2]], gr_repeat, type = 'any', ignore.strand = TRUE), stringsAsFactors = FALSE)
    inside_micro_dup = as.data.frame(GenomicRanges::findOverlaps(hits[[2]], gr_micro, type = 'any', ignore.strand = TRUE), stringsAsFactors = FALSE)
    dup_res_que = data.frame(hits[[2]][inside_repeat_dup$queryHits,], stringsAsFactors = FALSE)[,c(1:3,6:8)]
    dup_res_sub = dataTable[inside_repeat_dup$subjectHits,]
    dup_res_que_m = data.frame(hits[[2]][inside_micro_dup$queryHits,], stringsAsFactors = FALSE)[,c(1:3,6:8)]
    dup_res_sub_m = dataTable[inside_micro_dup$subjectHits,]
    inside_tab_dup = unique(cbind(dup_res_que, dup_res_sub)[,c(4,1:3,5,6,11,7:10,12:13)])
    inside_tab_dup_m = unique(cbind(dup_res_que_m, dup_res_sub_m)[,c(4,1:3,5,6,11,7:10,12,13)])
    names(inside_tab_dup) = c('q_name', 'q_chr', 'q_start', 'q_end', 'identity', 'align_length', 'rep_name', 'chr', 'start', 'end', 'strand', 'repClass', 'repFamily')
    names(inside_tab_dup_m) = c('q_name', 'q_chr', 'q_start', 'q_end', 'identity', 'align_length', 'rep_name', 'chr', 'start', 'end', 'strand', 'repClass', 'repFamily')
  } else {inside_tab_dup = NULL}
  message('- OK!')
  message('- Calculate distance.')
  strp = dataTable[which(dataTable$strand == '+'),]; strn = dataTable[which(dataTable$strand == '-'),]
  only_hits_tab = data.frame(hits[[1]], stringsAsFactors = FALSE)
  dist_only = lapply(c(1:nrow(only_hits_tab)), function(a){
    x = only_hits_tab$start[a]; y = only_hits_tab$query[a]; z = as.character(only_hits_tab$seqnames)[a]
    cal1p = strp$start-x; cal2p = strp$end-x; cal1n = strn$end-x; cal2n = strn$start-x
    cal1_ip = intersect(which(cal1p <= abs(range[1])), which(cal1p > 0)); cal2_ip = intersect(which(abs(cal2p) <= range[2]), which(cal2p < 0))
    cal1_in = intersect(which(cal1n >= range[1]), which(cal1n < 0)); cal2_in = intersect(which(cal2n <= range[2]), which(cal2n > 0))
    dat = data.frame(rbind(strp[c(cal1_ip, cal2_ip),], strn[c(cal1_in, cal2_in),]), dist = c(-c(cal1p[cal1_ip], cal2p[cal2_ip]), cal1n[cal1_in], cal2n[cal2_in]))
    dat = unique(data.frame('query' = rep(y, nrow(dat)), dat))
    dat = dat[which(dat$chr == z),]
    dat = dat[order(abs(dat$dist), decreasing = FALSE),][1,]
    cal1m = dataTable_micro$start-x; cal2m = dataTable_micro$end-x
    cal1_im = intersect(which(cal1m <= abs(range[1])), which(cal1m > 0))
    cal2_im = intersect(which(abs(cal2m) <= range[2]), which(cal2m < 0))
    dat_m = data.frame(dataTable[c(cal1_im, cal2_im),], dist = -c(cal1m[cal1_im], cal2m[cal2_im]))
    dat_m = unique(data.frame('query' = rep(y, nrow(dat_m)), dat_m))
    dat_m = dat_m[which(dat_m$chr == z),]
    dat_m = dat_m[order(abs(dat_m$dist), decreasing = FALSE),][1,]
    return(list(dat, dat_m))
  })
  hist_only_r = hist(unlist(lapply(dist_only, function(x){x[[1]]$dist/1000})), breaks = ranges/1000, plot = FALSE)
  hist_only_m = hist(unlist(lapply(dist_only, function(x){x[[2]]$dist/1000})), breaks = ranges/1000, plot = FALSE)
  if(length(hits[[2]]) != 0){
    dup_hits_tab = data.frame(hits[[2]], stringsAsFactors = FALSE)
    dist_dup = lapply(c(1:nrow(dup_hits_tab)), function(a){
      x = dup_hits_tab$start[a]; y = dup_hits_tab$query[a]; z = as.character(dup_hits_tab$seqnames)[a]
      cal1p = strp$start-x; cal2p = strp$end-x; cal1n = strn$end-x; cal2n = strn$start-x
      cal1_ip = intersect(which(cal1p <= abs(range[1])), which(cal1p > 0)); cal2_ip = intersect(which(abs(cal2p) <= range[2]), which(cal2p < 0))
      cal1_in = intersect(which(cal1n >= range[1]), which(cal1n < 0)); cal2_in = intersect(which(cal2n <= range[2]), which(cal2n > 0))
      dat = data.frame(rbind(strp[c(cal1_ip, cal2_ip),], strn[c(cal1_in, cal2_in),]), dist = c(-c(cal1p[cal1_ip], cal2p[cal2_ip]), cal1n[cal1_in], cal2n[cal2_in]))
      dat = unique(data.frame('query' = rep(y, nrow(dat)), dat))
      dat = dat[which(dat$chr == z),]
      dat = dat[order(abs(dat$dist), decreasing = FALSE),][1,]
      cal1m = dataTable_micro$start-x; cal2m = dataTable_micro$end-x
      cal1_im = intersect(which(cal1m <= abs(range[1])), which(cal1m > 0))
      cal2_im = intersect(which(abs(cal2m) <= range[2]), which(cal2m < 0))
      dat_m = data.frame(dataTable[c(cal1_im, cal2_im),], dist = -c(cal1m[cal1_im], cal2m[cal2_im]))
      dat_m = unique(data.frame('query' = rep(y, nrow(dat_m)), dat_m))
      dat_m = dat_m[which(dat_m$chr == z),]
      dat_m = dat_m[order(abs(dat_m$dist), decreasing = FALSE),][1,]
      return(list(dat, dat_m))
    })
    hist_dup_r = hist(unlist(lapply(dist_dup, function(x){x[[1]]$dist/1000})), breaks = ranges/1000, plot = FALSE)
    hist_dup_m = hist(unlist(lapply(dist_dup, function(x){x[[2]]$dist/1000})), breaks = ranges/1000, plot = FALSE)
  }
  message('- OK!')
  if(doRandom){
    message('- Do random set analysis.')
    seed = round(unclass(Sys.time()))
    ch_size = readRDS(file = system.file("extdata", paste0(organism, '_chrom.rds'), package = "RIPAT"))
    ch_size_num = as.numeric(ch_size$length)
    ch_start = cumsum(ch_size_num) - ch_size_num + 1
    ch_end = cumsum(ch_size_num)
    ch_ratio = (ch_size_num / sum(ch_size_num))
    random_set = sample(ch_size$chrom, size = randomSize, replace = TRUE, prob = ch_ratio)
    count_ch = plyr::count(random_set); row.names(count_ch) = count_ch$x
    count_ch = count_ch[ch_size$chrom,]
    count_ch = cbind(count_ch, ch_size_num)
    ran_set = apply(count_ch, 1, function(x){sample(c(1:x[3]), size = x[2], replace = FALSE, prob = NULL)})
    chr_ran = rep(paste0('chr', count_ch$x), count_ch$freq)
    pos_ran = unlist(ran_set)
    ran_tab = data.frame('Random' = c(1:randomSize), 'Random_chr' = chr_ran, 'Random_pos' = pos_ran, stringsAsFactors = FALSE)
    tmp = ran_tab[,c(2,3,3)]
    names(tmp) = c('chr', 'start', 'end')
    gr_random = GenomicRanges::makeGRangesFromDataFrame(tmp)
    utils::write.table(ran_tab, file = paste0(outPath, '/', outFileName, '_Random_set_', organism, '_repeat.txt'), quote = FALSE, append = FALSE, sep = '\t', na = '', row.names = FALSE, col.names = TRUE)
    inside_repeat_ran = as.data.frame(GenomicRanges::findOverlaps(gr_random, gr_repeat, type = 'any', ignore.strand = TRUE), stringsAsFactors = FALSE)
    a = as.data.frame(gr_random[inside_repeat_ran$queryHits,], stringsAsFactors = FALSE)
    b = dataTable[inside_repeat_ran$subjectHits,]
    inside_repeat_ran_tab = unique(cbind(a,b)[,c(1:3, 10, 6:9, 11, 12)])
    names(inside_repeat_ran_tab) = c('q_chr', 'q_start', 'q_end', 'repeat', 'chr', 'start', 'end', 'strand', 'repClass', 'repFamily')
    inside_micro_ran = as.data.frame(GenomicRanges::findOverlaps(gr_random, gr_micro, type = 'any', ignore.strand = TRUE), stringsAsFactors = FALSE)
    a2 = as.data.frame(gr_random[inside_micro_ran$queryHits,], stringsAsFactors = FALSE)
    b2 = dataTable_micro[inside_micro_ran$subjectHits,]
    inside_micro_ran_tab = unique(cbind(a2,b2)[,c(1:3, 9, 6:8)])
    names(inside_micro_ran_tab) = c('q_chr', 'q_start', 'q_end', 'microsatellite', 'chr', 'start', 'end')
    dist_ran = lapply(c(1:nrow(ran_tab)), function(a){
      x = ran_tab$Random_pos[a]; y = ran_tab$Random[a]; z = ran_tab$Random_chr[a]
      cal1p = strp$start-x; cal2p = strp$end-x; cal1n = strn$end-x; cal2n = strn$start-x
      cal1_ip = intersect(which(cal1p <= abs(range[1])), which(cal1p > 0)); cal2_ip = intersect(which(abs(cal2p) <= range[2]), which(cal2p < 0))
      cal1_in = intersect(which(cal1n >= range[1]), which(cal1n < 0)); cal2_in = intersect(which(cal2n <= range[2]), which(cal2n > 0))
      dat = data.frame(rbind(strp[c(cal1_ip, cal2_ip),], strn[c(cal1_in, cal2_in),]), dist = c(-c(cal1p[cal1_ip], cal2p[cal2_ip]), cal1n[cal1_in], cal2n[cal2_in]))
      dat = unique(data.frame('query' = rep(y, nrow(dat)), dat))
      dat = dat[which(dat$chr == z),]
      dat = dat[order(abs(dat$dist), decreasing = FALSE),][1,]
      cal1m = dataTable_micro$start-x; cal2m = dataTable_micro$end-x
      cal1_im = intersect(which(cal1m <= abs(range[1])), which(cal1m > 0))
      cal2_im = intersect(which(abs(cal2m) <= range[2]), which(cal2m < 0))
      dat_m = data.frame(dataTable[c(cal1_im, cal2_im),], dist = -c(cal1m[cal1_im], cal2m[cal2_im]))
      dat_m = unique(data.frame('query' = rep(y, nrow(dat_m)), dat_m))
      dat_m = dat_m[which(dat_m$chr == z),]
      dat_m = dat_m[order(abs(dat_m$dist), decreasing = FALSE),][1,]
      return(list(dat, dat_m))
    })
    all_dist_ran = unlist(lapply(dist_ran, function(x){x[[1]]$dist}))
    all_dist_m_ran = unlist(lapply(dist_ran, function(x){x[[2]]$dist}))
    hist_obj_ran = hist(all_dist_ran, plot = FALSE, breaks = ranges)
    hist_obj_m_ran = hist(all_dist_m_ran, plot = FALSE, breaks = ranges)
    message('- OK!')
  } else {message('[WARN] Skip random set analysis.')}
  message('- Draw histograms.')
  all_dist_only = unlist(lapply(dist_only, function(x){x[[1]]$dist}))
  all_dist_m_only = unlist(lapply(dist_only, function(x){x[[2]]$dist}))
  if(length(hits[[2]]) != 0){
    all_dist_dup = unlist(lapply(dist_dup, function(x){x[[1]]$dist}))
    all_dist_m_dup = unlist(lapply(dist_dup, function(x){x[[2]]$dist}))
  } else {all_dist_dup = NULL; all_dist_t_dup = NULL}
  if(includeUndecided){
    hist_obj = hist(c(all_dist_only, all_dist_dup), plot = FALSE, breaks = ranges)
    hist_obj_m = hist(c(all_dist_m_only, all_dist_dup_m), plot = FALSE, breaks = ranges)
    inside_tab = data.frame(rbind(inside_tab_only, inside_tab_dup), stringsAsFactors = FALSE)
    inside_tab_m = data.frame(rbind(inside_tab_only_m, inside_tab_dup_m), stringsAsFactors = FALSE)
    r_dist = list('Decided' = lapply(dist_only, function(x){x[[1]]}), 'Undecided' = lapply(dist_dup, function(x){x[[1]]}))
    m_dist = list('Decided' = lapply(dist_only, function(x){x[[2]]}), 'Undecided' = lapply(dist_dup, function(x){x[[2]]}))
  } else {
    hist_obj = hist(all_dist_only, plot = FALSE, breaks = ranges)
    hist_obj_m = hist(all_dist_m_only, plot = FALSE, breaks = ranges)
    inside_tab = inside_tab_only; inside_tab_m = inside_tab_only_m;
    r_dist = list('Decided' = lapply(dist_only, function(x){x[[1]]}))
    m_dist = list('Decided' = lapply(dist_only, function(x){x[[2]]}))
  }
  if(doRandom){
    count_site = hist_obj$counts; count_site_ran = hist_obj_ran$counts
    count_m_site = hist_obj_m$counts; count_m_site_ran = hist_obj_m_ran$counts
    if(includeUndecided){count_all = sum(length(all_dist_only), length(all_dist_dup))} else {count_all = length(all_dist_only)}
    count_data = data.frame('Range' = factor(rep(ranges[ranges != 0]/1000, 2), levels = ranges[ranges != 0]/1000), 'Group' = c(rep('Observed', length(count_site)), rep('Random', length(count_site_ran))), 'Count' = c(count_site, count_site_ran), 'Freq' = c(count_site/count_all, count_site_ran/randomSize))
    count_m_data = data.frame('Range' = factor(rep(ranges[ranges != 0]/1000, 2), levels = ranges[ranges != 0]/1000), 'Group' = c(rep('Observed', length(count_m_site)), rep('Random', length(count_m_site_ran))), 'Count' = c(count_m_site, count_m_site_ran), 'Freq' = c(count_m_site/count_all, count_m_site_ran/randomSize))
  } else {
    count_site = hist_obj$counts; count_m_site = hist_obj_m$counts
    if(includeUndecided){count_all = sum(length(all_dist_only), length(all_dist_dup))} else {count_all = length(all_dist_only)}
    count_data = data.frame('Range' = factor(ranges[ranges != 0]/1000, levels = ranges[ranges != 0]/1000), 'Group' = rep('Observed', length(count_site)), 'Count' = count_site, 'Freq' = count_site/count_all)
    count_m_data = data.frame('Range' = factor(ranges[ranges != 0]/1000, levels = ranges[ranges != 0]/1000), 'Group' = rep('Observed', length(count_m_site)), 'Count' = count_m_site, 'Freq' = count_m_site/count_all)
  }
  grDevices::png(paste0(outPath, '/', outFileName, '_distribution_repeat_', organism, '.png'), width = 1200, height = 750)
  r_plot = ggplot2::ggplot(count_data) + ggplot2::geom_bar(ggplot2::aes(x = Range, y = Freq, fill = Group), stat = "identity", position = "dodge", width = 0.5) +
    ggplot2::lims(y = c(0, max(count_data$Freq)*1.5)) + ggplot2::ggtitle(label = "Random distribution (Repeat)") +
    ggplot2::xlab('Intervals (Kbs)') + ggplot2::ylab("Ratio of Integration Events") + ggplot2::scale_fill_manual(values = c('mediumspringgreen', 'mediumpurple')) +
    ggplot2::theme(panel.background = ggplot2::element_rect(fill="white", colour = "white"), panel.grid.major = ggplot2::element_line(size = 0.5, linetype = 'dotted', colour = 'black'),
                   axis.line = ggplot2::element_line(colour = "darkgrey"), legend.title = ggplot2::element_blank(),
                   legend.key.size = ggplot2::unit(0.7, "cm"), plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 15),
                   legend.text = ggplot2::element_text(size = 13), axis.text.x = ggplot2::element_text(size = 12),
                   axis.text.y = ggplot2::element_text(size = 12), axis.title.x = ggplot2::element_text(size = 13), axis.title.y = ggplot2::element_text(size = 13))
  print(r_plot)
  grDevices::dev.off()
  grDevices::png(paste0(outPath, '/', outFileName, '_distribution_microsatellite_', organism, '.png'), width = 1200, height = 750)
  m_plot = ggplot2::ggplot(count_m_data) + ggplot2::geom_bar(ggplot2::aes(x = Range, y = Freq, fill = Group), stat = "identity", position = "dodge", width = 0.5) +
    ggplot2::lims(y = c(0, max(count_m_data$Freq)*1.5)) + ggplot2::ggtitle(label = "Random distribution (Microsatellite)") +
    ggplot2::xlab('Intervals (Kbs)') + ggplot2::ylab("Ratio of Integration Events") + ggplot2::scale_fill_manual(values = c('mediumspringgreen', 'mediumpurple')) +
    ggplot2::theme(panel.background = ggplot2::element_rect(fill="white", colour = "white"), panel.grid.major = ggplot2::element_line(size = 0.5, linetype = 'dotted', colour = 'black'),
                   axis.line = ggplot2::element_line(colour = "darkgrey"), legend.title = ggplot2::element_blank(),
                   legend.key.size = ggplot2::unit(0.7, "cm"), plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 15),
                   legend.text = ggplot2::element_text(size = 13), axis.text.x = ggplot2::element_text(size = 12),
                   axis.text.y = ggplot2::element_text(size = 12), axis.title.x = ggplot2::element_text(size = 13), axis.title.y = ggplot2::element_text(size = 13))
  print(m_plot)
  grDevices::dev.off()
  message('- OK!')
  result_list = list(inside_tab, r_dist, count_data, gr_repeat, inside_tab_m, m_dist, count_m_data, gr_micro, organism)
  names(result_list) = c('Repeat_inside', 'Repeat_distribution', 'Repeat_plot_data', 'Repeat_data', 'Micro_inside', 'Micro_distribution', 'Micro_plot_data', 'Micro_data', 'Target_ver')
  if(doRandom){
    result_list = c(result_list, list(inside_repeat_ran_tab, inside_micro_ran_tab, dist_ran))
    names(result_list)[10:12] = c('Repeat_inside_ran', 'Micro_inside_ran', 'Random_distribution')
  }
  result_list = c(result_list, list('Repeat_hist_plot' = r_plot, 'Microsatellite_hist_plot' = m_plot))
  message('----- Finish. (Time : ', date(), ')')
  return(result_list)
}
