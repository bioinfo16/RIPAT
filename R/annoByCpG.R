#' @title Annotate integration sites by CpG sites.
#' 
#' @description
#' Annotate vector integration sites by CpG site data.
#'  
#' @usage 
#' annoByCpG(hits, mapTool = 'blast', organism = 'GRCh37', interval = 5000, 
#'           range = c(-20000, 20000), doRandom = TRUE,
#'           randomSize = if(doRandom){10000}else{NULL}, 
#'           includeUndecided = FALSE, outPath = getwd(),
#'           outFileName = paste0('RIPAT', round(unclass(Sys.time()))))
#' 
#' @param hits a GR object. This object made by \code{makeInputObj} function.
#' @param mapTool a single character. Function serves two types of object
#'                such as outputs from BLAST and BLAT.
#'                Default is 'blast'. If you want to use BLAT result, use 'blat'.
#' @param organism a single character. This function can run by two versions of organisms
#'                 such as GRCh37, GRCh38 (Human). Default is 'GRCh37'.
#' @param interval an integer vector. This number means interval number for
#'                 distribution analysis. Default is 5000.
#' @param range an integer array. The range of highlight region for analysis.
#'              Default range is c(-20000, 20000).
#' @param doRandom TRUE or FALSE. If user types TRUE, random set is generated
#'                 and user can do random distribution analysis. Default is TRUE.
#'                 If this value is FALSE, random distribution analysis is not executed.
#' @param randomSize an integer vector. A random set size. Default is 10000.
#' @param includeUndecided TRUE or FALSE. If user want to use undecided hits in analysis,
#'                         enter TRUE. Default is FALSE.
#' @param outPath an string vector. Plots are saved in this path. Default value is R home directory.
#' @param outFileName a character vector. Attached ID to the result file name.
#' 
#' @return Return a result list that is made up of insertion and distribution result tables
#'         and GenomicRange object of CpG data.
#'         
#' @examples 
#' data(blast_obj); data(cpg_exam_db)
#' saveRDS(cpg_exam_db, paste0(system.file("extdata", package = 'RIPAT'), '/GRCh37_cpg.rds'))
#' 
#' blast_cpg = annoByCpG(hits = blast_obj, doRandom = FALSE, outFileName = 'blast_res')
#'            
#' @export
annoByCpG = function(hits, mapTool = 'blast', organism = 'GRCh37', interval = 5000, range = c(-20000, 20000), doRandom = TRUE, randomSize = if(doRandom){10000}else{NULL}, includeUndecided = FALSE, outPath = getwd(), outFileName = paste0('RIPAT', round(unclass(Sys.time())))){
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
  dbFile = paste0('/', organism, '_cpg.rds')
  dataTable = data.frame(readRDS(paste0(dbPath, dbFile)), stringsAsFactors = FALSE)
  dataTable$chrom = as.character(dataTable$chrom)
  gr_cpgs = GenomicRanges::makeGRangesFromDataFrame(dataTable, keep.extra.columns = TRUE, ignore.strand = TRUE)
  message('- OK!')
  message('- Find integration sites located in CpG sites.')
  inside_cpg_only = as.data.frame(GenomicRanges::findOverlaps(hits[[1]], gr_cpgs, type = 'any', ignore.strand = TRUE), stringsAsFactors = FALSE)
  only_res_que = data.frame(hits[[1]][inside_cpg_only$queryHits,], stringsAsFactors = FALSE)[,c(1:3,6:8)]
  only_res_sub = dataTable[inside_cpg_only$subjectHits,]
  inside_tab_only = unique(cbind(only_res_que, only_res_sub)[,c(4,1:3,5,6,11,8:10,12,13,15,16)])
  names(inside_tab_only) = c('q_name', 'q_chr', 'q_start', 'q_end', 'identity', 'align_length', 'CpG_name', 'chr', 'start', 'end', 'length', 'cpgNum', 'perCpg', 'perGc')
  if(length(hits[[2]]) != 0){
    inside_cpg_dup = as.data.frame(GenomicRanges::findOverlaps(hits[[2]], gr_cpgs, type = 'any', ignore.strand = TRUE), stringsAsFactors = FALSE)
    dup_res_que = data.frame(hits[[2]][inside_cpg_dup$queryHits,], stringsAsFactors = FALSE)[,c(1:3,6:8)]
    dup_res_sub = dataTable[inside_cpg_dup$subjectHits,]
    inside_tab_dup = unique(cbind(dup_res_que, dup_res_sub)[,c(4,1:3,5,6,11,8:10,12,13,15,16)])
    names(inside_tab_dup) = c('q_name', 'q_chr', 'q_start', 'q_end', 'identity', 'align_length', 'CpG_name', 'chr', 'start', 'end', 'length', 'cpgNum', 'perCpg', 'perGc')
  } else {inside_tab_dup = NULL}
  message('- OK!')
  message('- Calculate distance.')
  only_hits_tab = data.frame(hits[[1]], stringsAsFactors = FALSE)
  dist_only = lapply(c(1:nrow(only_hits_tab)), function(a){
    x = only_hits_tab$start[a]; y = only_hits_tab$query[a]; z = as.character(only_hits_tab$seqnames)[a]
    cal1 = dataTable$start-x; cal2 = dataTable$end-x
    cal1_i = intersect(which(cal1 <= abs(range[1])), which(cal1 > 0)); cal2_i = intersect(which(abs(cal2) <= range[2]), which(cal2 < 0))
    dat = data.frame(dataTable[c(cal1_i, cal2_i),], dist = -c(cal1[cal1_i], cal2[cal2_i]))
    dat = unique(data.frame('query' = rep(y, nrow(dat)), dat))
    dat = dat[which(dat$chr == z),]
    dat = dat[order(abs(dat$dist), decreasing = FALSE),][1,]
    return(dat)
  })
  hist_only_c = hist(unlist(lapply(dist_only, function(x){x$dist/1000})), breaks = ranges/1000, plot = FALSE)
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
    hist_dup_c = hist(unlist(lapply(dist_dup, function(x){x$dist/1000})), breaks = ranges/1000, plot = FALSE)
  }
  message('- OK!')
  if(doRandom){
    message('- Do random set analysis.')
    seed = round(unclass(Sys.time()))
    ch_size = readRDS(file = system.file("extdata", paste0(organism, '_chrom.rds'), package = "RIPAT"))
    ch_size_num = as.numeric(ch_size$length)
    ch_start = cumsum(ch_size_num) - ch_size_num + 1
    ch_end = cumsum(ch_size_num)
    ch_ratio = ch_size_num / sum(ch_size_num)
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
    utils::write.table(ran_tab, file = paste0(outPath, '/', outFileName, '_Random_set_', organism, '_cpg.txt'), quote = FALSE, append = FALSE, sep = '\t', na = '', row.names = FALSE, col.names = TRUE)
    inside_cpg_ran = as.data.frame(GenomicRanges::findOverlaps(gr_random, gr_cpgs, type = 'any', ignore.strand = TRUE), stringsAsFactors = FALSE)
    a = as.data.frame(gr_random[inside_cpg_ran$queryHits,], stringsAsFactors = FALSE)
    b = dataTable[inside_cpg_ran$subjectHits,]
    inside_ran_tab = unique(cbind(a,b)[,c(1:3,10,7:9,11,12,14,15)])
    names(inside_ran_tab) = c('q_chr', 'q_start', 'q_end', 'CpG_name', 'chr', 'start', 'end', 'length', 'cpgNum', 'perCpg', 'perGc')
    c_dist_ran = lapply(c(1:nrow(ran_tab)), function(a){
      x = ran_tab$Random_pos[a]; y = ran_tab$Random[a]; z = ran_tab$Random_chr[a]
      cal1 = dataTable$start-x; cal2 = dataTable$end-x
      cal1_i = intersect(which(cal1 <= abs(range[1])), which(cal1 > 0)); cal2_i = intersect(which(abs(cal2) <= range[2]), which(cal2 < 0))
      dat = data.frame(dataTable[c(cal1_i, cal2_i),], dist = -c(cal1[cal1_i], cal2[cal2_i]))
      dat = unique(data.frame('query' = rep(y, nrow(dat)), dat))
      dat = dat[which(dat$chr == z),]
      dat = dat[order(abs(dat$dist), decreasing = FALSE),][1,]
      return(dat)
    })
    all_dist_ran = unlist(lapply(c_dist_ran, function(x){x$dist}))
    hist_obj_ran = hist(all_dist_ran, plot = FALSE, breaks = seq(from = range[1], to = range[2], by = interval))
    message('- OK!')
  } else {message('[WARN] Skip random set analysis.')}
  message('- Draw histograms.')
  all_dist_only = unlist(lapply(dist_only, function(x){x$dist}))
  if(length(hits[[2]]) != 0){
    all_dist_dup = unlist(lapply(dist_dup, function(x){x$dist}))
  } else {all_dist_dup = NULL}
  if(includeUndecided){
    hist_obj = hist(c(all_dist_only, all_dist_dup), plot = FALSE, breaks = ranges)
    inside_tab = data.frame(rbind(inside_tab_only, inside_tab_dup), stringsAsFactors = FALSE)
    c_dist = list('Decided' = dist_only, 'Undecided' = dist_dup)
  } else {
    hist_obj = hist(all_dist_only, plot = FALSE, breaks = ranges)
    inside_tab = inside_tab_only
    c_dist = list('Decided' = dist_only)
  }
  if(doRandom){
    count_site = hist_obj$counts; count_site_ran = hist_obj_ran$counts
    if(includeUndecided){count_all = sum(length(all_dist_only), length(all_dist_dup))} else {count_all = length(all_dist_only)}
    count_data = data.frame('Range' = factor(rep(ranges[ranges != 0]/1000, 2), levels = ranges[ranges != 0]/1000), 'Group' = c(rep('Observed', length(count_site)), rep('Random', length(count_site_ran))), 'Count' = c(count_site, count_site_ran), 'Freq' = c(count_site/count_all, count_site_ran/randomSize))
  } else {
    count_site = hist_obj$counts
    if(includeUndecided){count_all = sum(length(all_dist_only), length(all_dist_dup))} else {count_all = length(all_dist_only)}
    count_data = data.frame('Range' = factor(ranges[ranges != 0]/1000, levels = ranges[ranges != 0]/1000), 
                            'Group' = rep('Observed', length(count_site)), 'Count' = count_site, 'Freq' = count_site/count_all)
  }
  grDevices::png(paste0(outPath, '/', outFileName, '_distribution_CpG_', organism, '.png'), width = 1200, height = 750)
  c_plot = ggplot2::ggplot(count_data) + ggplot2::geom_bar(ggplot2::aes(x = Range, y = Freq, fill = Group), stat = "identity", position = "dodge", width = 0.5) +
    ggplot2::lims(y = c(0, max(count_data$Freq)*1.5)) + ggplot2::ggtitle(label = "Random distribution (CpG)") +
    ggplot2::xlab('Intervals (Kbs)') + ggplot2::ylab("Ratio of Integration Events") + ggplot2::scale_fill_manual(values = c('mediumspringgreen', 'mediumpurple')) +
    ggplot2::theme(panel.background = ggplot2::element_rect(fill="white", colour = "white"), panel.grid.major = ggplot2::element_line(size = 0.5, linetype = 'dotted', colour = 'black'),
                   axis.line = ggplot2::element_line(colour = "darkgrey"), legend.title = ggplot2::element_blank(),
                   legend.key.size = ggplot2::unit(0.7, "cm"), plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 15),
                   legend.text = ggplot2::element_text(size = 13), axis.text.x = ggplot2::element_text(size = 12),
                   axis.text.y = ggplot2::element_text(size = 12), axis.title.x = ggplot2::element_text(size = 13), axis.title.y = ggplot2::element_text(size = 13))
  print(c_plot)
  grDevices::dev.off()
  message('- OK!')
  result_list = list(inside_tab, c_dist, count_data, gr_cpgs, organism)
  names(result_list) = c('CpG_inside', 'CpG_distribution', 'CpG_plot_data', 'CpG_data', 'Target_ver')
  if(doRandom){
    result_list = c(result_list, list(inside_ran_tab, c_dist_ran))
    names(result_list)[6:7] = c('CpG_inside_ran', 'Random_distribution')
  }
  result_list = c(result_list, list('CpG_hist_plot' = c_plot))
  message('----- Finish. (Time : ', date(), ')')
  return(result_list)
}
