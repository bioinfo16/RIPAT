#' @title Annotate integration sites by genes and TSSs.
#' 
#' @description
#' Annotate vector integration sites by gene data.
#' 
#' @usage 
#' annoByGene(hits, ran_hits = NULL,
#'            mapTool = 'blast',
#'            organism = 'GRCh37',
#'            interval = 5000,
#'            range = c(-20000, 20000),
#'            outPath = getwd(),
#'            outFileName = paste0('RIPAT', round(unclass(Sys.time()))))
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
#' @return Return a result list that is made up of insertion and distribution result tables\cr and GenomicRange object of Gene and TSS data.
#'         
#' @examples 
#' data(blast_obj); data(gene_exam_db); data(tss_exam_db)
#' saveRDS(gene_exam_db, 
#'         paste0(system.file("extdata", package = 'RIPAT'), 
#'                '/GRCh37_gene.rds'))
#' saveRDS(tss_exam_db, 
#'         paste0(system.file("extdata", package = 'RIPAT'),
#'                '/GRCh37_TSS.rds'))
#'
#' blast_gene = annoByGene(hits = blast_obj, ran_hits = NULL, 
#'                         outFileName = 'blast_res')
#' 
#' @export
annoByGene = function(hits, ran_hits = NULL, mapTool = 'blast', organism = 'GRCh37', interval = 5000, range = c(-20000, 20000), outPath = getwd(), outFileName = paste0('RIPAT', round(unclass(Sys.time())))){
  message('----- Annotate integration sites. (Time : ', date(), ')')
  message('- Validate options')
  if(length(which(c('blast', 'blat') %in% mapTool)) == 0){
    stop("[ERROR] Please confirm the alignment tool name.\n----- This process is halted. (Time : ", date(), ")\n")
  }
  if(length(which(c('GRCh37', 'GRCh38') %in% organism)) == 0){
    stop("[ERROR] Please use GRCh37/GRCh38 data.\n----- This process is halted. (Time : ", date(), ")\n")
  }
  if(stringr::str_ends(outPath, pattern = '/')){outPath = stringr::word(outPath, start = 1, end = nchar(outPath), sep = '')}
  if(range[1] - range[2] >= 0){
    stop("[ERROR] Please check distribution range.\n----- This process is halted. (Time : ", date(), ")\n")
  } else {ranges = seq(from = range[1], to = range[2], by = interval)}
  message('- OK!')
  message('- Load the dataset.')
  dbPath = system.file("extdata", package = "RIPAT")
  dbFile = paste0('/', organism, '_gene.rds')
  dataTable = data.frame(readRDS(paste0(dbPath, dbFile)), stringsAsFactors = FALSE)
  dataTable = cbind(dataTable, 'g_len' = abs(dataTable$end - dataTable$start) + 1)
  dbFile2 = paste0('/', organism, '_TSS.rds')
  dataTable_tss = data.frame(readRDS(paste0(dbPath, dbFile2)), stringsAsFactors = FALSE)
  gr_genes = GenomicRanges::makeGRangesFromDataFrame(dataTable, keep.extra.columns = TRUE, ignore.strand = FALSE)
  gr_tss = GenomicRanges::makeGRangesFromDataFrame(dataTable_tss, keep.extra.columns = TRUE, ignore.strand = FALSE)
  message('- OK!')
  message('- Find integration sites located in genes.')
  inside_gene_only = as.data.frame(GenomicRanges::findOverlaps(hits[[1]], gr_genes, type = 'any', ignore.strand = TRUE), stringsAsFactors = FALSE)
  only_res_que = data.frame(hits[[1]][inside_gene_only$queryHits,], stringsAsFactors = FALSE)[,c(1:3,6:8)]
  only_res_sub = dataTable[inside_gene_only$subjectHits,]
  inside_tab = unique(cbind(only_res_que, only_res_sub)[,c(5,1:4,6:14)])
  names(inside_tab)[c(7:14)] = c('Target_gene', 'Target_ensembl_id', 'Target_gene_type', 'Target_chr', 'Target_start', 'Target_end', 'Target_strand', 'Target_gene_length')
  message('- OK!')
  message('- Calculate distance.')
  strp = dataTable[which(dataTable$strand == '+'),]; strn = dataTable[which(dataTable$strand == '-'),]
  only_hits_tab = data.frame(hits[[1]], stringsAsFactors = FALSE)
  dist_only = lapply(c(1:nrow(only_hits_tab)), function(a){
    x = only_hits_tab$start[a]; y = only_hits_tab$query[a]; z = only_hits_tab$seqnames[a]
    cal1p = strp$start-x; cal2p = strp$end-x; cal1n = strn$end-x; cal2n = strn$start-x
    cal1_ip = intersect(which(cal1p <= abs(range[1])), which(cal1p > 0)); cal2_ip = intersect(which(abs(cal2p) <= range[2]), which(cal2p < 0))
    cal1_in = intersect(which(cal1n >= range[1]), which(cal1n < 0)); cal2_in = intersect(which(cal2n <= range[2]), which(cal2n > 0))
    dat = data.frame(rbind(strp[c(cal1_ip, cal2_ip),], strn[c(cal1_in, cal2_in),]), dist = c(-c(cal1p[cal1_ip], cal2p[cal2_ip]), cal1n[cal1_in], cal2n[cal2_in]))
    dat = unique(data.frame('query' = rep(y, nrow(dat)), dat))
    dat = dat[which(dat$chr == z),]
    dat = dat[order(abs(dat$dist), decreasing = FALSE),][1,]
    cal_tss = dataTable_tss$start-x
    cal_i1_tss = intersect(which(cal_tss <= abs(range[1])), which(cal_tss > 0))
    cal_i2_tss = intersect(which(abs(cal_tss) <= range[2]), which(cal_tss < 0))
    dat_tss = cbind(dataTable_tss[c(cal_i1_tss, cal_i2_tss),], dist = -c(cal_tss[cal_i1_tss], cal_tss[cal_i2_tss]))
    dat_tss = unique(data.frame('query' = rep(y, nrow(dat_tss)), dat_tss))
    dat_tss = dat_tss[which(dat_tss$chr == z),]
    dat_tss = dat_tss[order(abs(dat_tss$dist), decreasing = FALSE),][1,]
    return(list(dat, dat_tss))
  })
  hist_only_g = hist(unlist(lapply(dist_only, function(x){x[[1]]$dist/1000})), breaks = ranges/1000, plot = FALSE)
  hist_only_t = hist(unlist(lapply(dist_only, function(x){x[[2]]$dist/1000})), breaks = ranges/1000, plot = FALSE)
  message('- OK!')
  if(!is.null(ran_hits)){
    message('- Do random set analysis.')
    randomSize = length(ran_hits); gr_random = ran_hits
    random_set = data.frame(c(1:randomSize), data.frame(ran_hits), stringsAsFactors = FALSE)[,c(1:3)]; names(random_set) = c('Random', 'Random_chr', 'Random_pos')
    inside_gene_ran = as.data.frame(GenomicRanges::findOverlaps(gr_random, gr_genes, type = 'any', ignore.strand = TRUE), stringsAsFactors = FALSE)
    a = as.data.frame(gr_random[inside_gene_ran$queryHits,], stringsAsFactors = FALSE)
    b = dataTable[inside_gene_ran$subjectHits,]
    inside_ran_tab = unique(cbind(a,b)[,c(1,2,3,6:13)])
    names(inside_ran_tab) = c('q_chr', 'q_start', 'q_end', 'symbol', 'ensembl_id', 'gene_type', 'chr', 'start', 'end', 'strand', 'percentage_gc_content')
    dist_ran = lapply(c(1:nrow(random_set)), function(a){
      x = as.numeric(random_set$Random_pos[a]); y = random_set$Random[a]; z = random_set$Random_chr[a]
      cal1p = strp$start-x; cal2p = strp$end-x; cal1n = strn$end-x; cal2n = strn$start-x
      cal1_ip = intersect(which(cal1p <= abs(range[1])), which(cal1p > 0)); cal2_ip = intersect(which(abs(cal2p) <= range[2]), which(cal2p < 0))
      cal1_in = intersect(which(cal1n >= range[1]), which(cal1n < 0)); cal2_in = intersect(which(cal2n <= range[2]), which(cal2n > 0))
      dat = data.frame(rbind(strp[c(cal1_ip, cal2_ip),], strn[c(cal1_in, cal2_in),]), dist = c(-c(cal1p[cal1_ip], cal2p[cal2_ip]), cal1n[cal1_in], cal2n[cal2_in]))
      dat = unique(data.frame('query' = rep(y, nrow(dat)), dat))
      dat = dat[which(dat$chr == z),]
      dat = dat[order(abs(dat$dist), decreasing = FALSE),][1,]
      cal_tss = dataTable_tss$start-x
      cal_i1_tss = intersect(which(cal_tss <= abs(range[1])), which(cal_tss > 0))
      cal_i2_tss = intersect(which(abs(cal_tss) <= range[2]), which(cal_tss < 0))
      dat_tss = cbind(dataTable_tss[c(cal_i1_tss, cal_i2_tss),], dist = c(cal_tss[cal_i1_tss], cal_tss[cal_i2_tss]))
      dat_tss = unique(data.frame('query' = rep(y, nrow(dat_tss)), dat_tss))
      dat_tss = dat_tss[which(dat_tss$chr == z),]
      dat_tss = dat_tss[order(abs(dat_tss$dist), decreasing = FALSE),][1,]
      return(list(dat, dat_tss))
    })
    all_dist_ran = unlist(lapply(dist_ran, function(x){x[[1]]$dist})); all_dist_t_ran = unlist(lapply(dist_ran, function(x){x[[2]]$dist}))
    hist_obj_ran = hist(all_dist_ran, plot = TRUE, breaks = ranges); hist_t_obj_ran = hist(all_dist_t_ran, plot = T, breaks = ranges)
    message('- OK!')
  } else {message('[WARN] Skip random set analysis.')}
  message('- Draw histograms.')
  all_dist_only = unlist(lapply(dist_only, function(x){x[[1]]$dist}))
  all_dist_t_only = unlist(lapply(dist_only, function(x){x[[2]]$dist}))
  hist_obj = hist(all_dist_only, plot = FALSE, breaks = ranges)
  hist_t_obj = hist(all_dist_t_only, plot = FALSE, breaks = ranges)
  g_dist = list('Decided' = lapply(dist_only, function(x){x[[1]]})); t_dist = list('Decided' = lapply(dist_only, function(x){x[[2]]}))
  if(!is.null(ran_hits)){
    count_site = hist_obj$counts; count_site_ran = hist_obj_ran$counts
    count_t_site = hist_t_obj$counts; count_t_site_ran = hist_t_obj_ran$counts
    count_all = nrow(only_hits_tab)
    count_data = data.frame('Range' = factor(rep(ranges[ranges != 0]/1000, 2), levels = ranges[ranges != 0]/1000), 'Group' = c(rep('Observed', length(count_site)), rep('Random', length(count_site_ran))), 'Count' = c(count_site, count_site_ran), 'Freq' = c(count_site/count_all, count_site_ran/randomSize))
    count_t_data = data.frame('Range' = factor(rep(ranges[ranges != 0]/1000, 2), levels = ranges[ranges != 0]/1000), 'Group' = c(rep('Observed', length(count_t_site)), rep('Random', length(count_t_site_ran))), 'Count' = c(count_t_site, count_t_site_ran), 'Freq' = c(count_t_site/count_all, count_t_site_ran/randomSize))
  } else {
    count_site = hist_obj$counts; count_t_site = hist_t_obj$counts
    count_all = nrow(only_hits_tab)
    count_data = data.frame('Range' = factor(ranges[ranges != 0]/1000, levels = ranges[ranges != 0]/1000), 'Group' = rep('Observed', length(count_site)), 'Count' = count_site, 'Freq' = count_site/count_all)
    count_t_data = data.frame('Range' = factor(ranges[ranges != 0]/1000, levels = ranges[ranges != 0]/1000), 'Group' = rep('Observed', length(count_t_site)), 'Count' = count_t_site, 'Freq' = count_t_site/count_all)
  }
  grDevices::png(paste0(outPath, '/', outFileName, '_distribution_gene_', organism, '.png'), width = 1200, height = 750)
  g_plot = ggplot2::ggplot(count_data) + ggplot2::geom_bar(ggplot2::aes(x = Range, y = Freq, fill = Group), stat = "identity", position = "dodge", width = 0.5) +
    ggplot2::lims(y = c(0, max(count_data$Freq)*1.5)) + ggplot2::ggtitle(label = "Random distribution (Gene)") +
    ggplot2::xlab('Intervals (Kbs)') + ggplot2::ylab("Ratio of Integration Events") + ggplot2::scale_fill_manual(values = c('mediumspringgreen', 'mediumpurple')) +
    ggplot2::theme(panel.background = ggplot2::element_rect(fill="white", colour = "white"), panel.grid.major = ggplot2::element_line(size = 0.5, linetype = 'dotted', colour = 'black'),
                   axis.line = ggplot2::element_line(colour = "darkgrey"), legend.title = ggplot2::element_blank(), 
                   legend.key.size = ggplot2::unit(0.7, "cm"), plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 20),
                   legend.text = ggplot2::element_text(size = 18), axis.text = ggplot2::element_text(size = 17),
                   axis.title = ggplot2::element_text(size = 18))
  print(g_plot)
  grDevices::dev.off()
  grDevices::png(paste0(outPath, '/', outFileName, '_distribution_tss_', organism, '.png'), width = 1200, height = 750)
  t_plot = ggplot2::ggplot(count_t_data) + ggplot2::geom_bar(ggplot2::aes(x = Range, y = Freq, fill = Group), stat = "identity", position = "dodge", width = 0.5) +
    ggplot2::lims(y = c(0, max(count_t_data$Freq)*1.5)) + ggplot2::ggtitle(label = "Random distribution (TSS)") +
    ggplot2::xlab('Intervals (Kbs)') + ggplot2::ylab("Ratio of Integration Events") + ggplot2::scale_fill_manual(values = c('mediumspringgreen', 'mediumpurple')) +
    ggplot2::theme(panel.background = ggplot2::element_rect(fill="white", colour = "white"), panel.grid.major = ggplot2::element_line(size = 0.5, linetype = 'dotted', colour = 'black'),
                   axis.line = ggplot2::element_line(colour = "darkgrey"), legend.title = ggplot2::element_blank(), 
                   legend.key.size = ggplot2::unit(0.7, "cm"), plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 20),
                   legend.text = ggplot2::element_text(size = 18), axis.text = ggplot2::element_text(size = 17),
                   axis.title = ggplot2::element_text(size = 18))
  print(t_plot)
  grDevices::dev.off()
  message('- OK!')
  result_list = list(inside_tab, g_dist, count_data, gr_genes, t_dist, count_t_data, gr_tss, organism)
  names(result_list) = c('Gene_inside', 'Gene_distribution', 'Gene_plot_data', 'Gene_data', 'TSS_distribution', 'TSS_plot_data', 'TSS_data', 'Target_ver')
  if(!is.null(ran_hits)){
    result_list = c(result_list, list(inside_ran_tab, dist_ran))
    names(result_list)[9:10] = c('Gene_inside_ran', 'Random_distribution')
  }
  result_list = c(result_list, list('Gene_hist_plot' = g_plot, 'TSS_hist_plot' = t_plot))
  message('----- Finish. (Time : ', date(), ')')
  return(result_list)
}