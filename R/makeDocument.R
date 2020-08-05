#' @title Make result documents as excel files.
#' 
#' @description \preformatted{
#' User can make the excel file for output data preservation by this function.
#' }
#' 
#' @usage 
#' makeDocument(res, resType, includeUndecided = FALSE, interval = 5000, 
#'              range = c(-20000, 20000), outPath = getwd(),
#'              outFileName = paste0('RIPAT', round(unclass(Sys.time()))))
#' 
#' @param res a GR object. This object is output of \code{annoByGene, annoByCpG, annoByRepeat, annoByVar} function.
#' @param resType a character vector. User enter the annotation type of input
#'                such as gene, cpg, repeat and variant.
#' @param includeUndecided TRUE or FALSE. If user want to use undecided hits in analysis, enter TRUE.
#'                         Default is FALSE.
#' @param interval an integer vector. This number means interval number for distribution analysis. Default is 5000.
#' @param range an integer array. It means the range for highlight region of this analysis. Default range is c(-20000, 20000).
#' @param outPath an string vector. Plots are saved in this path. Default value is R home directory.
#' @param outFileName a character vector. Attached ID to the result file name.
#' @return make excel file about vector integration sites and proportion test result.
#' 
#' @examples 
#' data(blast_gene)
#' makeDocument(res = blast_gene, resType = 'gene', interval = 5000, 
#'              range = c(-20000, 20000), outFileName = 'blast_gene_res')
#' 
#'
#' @export
makeDocument = function(res, resType, includeUndecided = FALSE, interval = 5000, range = c(-20000, 20000), outPath = getwd(), outFileName = paste0('RIPAT', round(unclass(Sys.time())))){
  message('----- Make result documents. (Time : ', date(), ')')
  if(stringr::str_ends(outPath, pattern = '/')){
    outPath = stringr::word(outPath, start = 1, end = nchar(outPath), sep = '')
  } else {}
  if(range[1] - range[2] >= 0){
    stop("[ERROR] Please check distribution range.\n----- This process is halted. (Time : ", date(), ")\n")
  } else {ranges = seq(from = range[1], to = range[2], by = interval)}
  message('- OK!')
  message('- Edit experimental data.')
  if(resType == 'gene'){
    dist_dat = data.frame(do.call("rbind", res$Gene_distribution$Decided), stringsAsFactors = FALSE)
    dist_dat2 = data.frame(do.call("rbind", res$TSS_distribution$Decided), stringsAsFactors = FALSE)
    dist_dat = dist_dat[-which(is.na(dist_dat$query)),]; dist_dat2 = dist_dat2[-which(is.na(dist_dat2$query)),]
    if(includeUndecided){
      dist_dat_u = data.frame(do.call("rbind", res$Gene_distribution$Undecided), stringsAsFactors = FALSE)
      dist_dat_u2 = data.frame(do.call("rbind", res$TSS_distribution$Undecided), stringsAsFactors = FALSE)
      dist_dat_u = dist_dat_u[-which(is.na(dist_dat_u$query)),]; dist_dat_u2 = dist_dat_u2[-which(is.na(dist_dat_u2$query)),]
      dist_dat = data.frame(rbind(dist_dat, dist_dat_u), stringsAsFactors = FALSE)
      dist_dat2 = data.frame(rbind(dist_dat2, dist_dat_u2), stringsAsFactors = FALSE)
    } else {}
    types = unique(dist_dat$gene_type)
    type_tab = lapply(types, function(a){subset(dist_dat, dist_dat$gene_type == a)})
    names(type_tab) = types
    type_len = sort(unlist(lapply(type_tab, nrow)), decreasing = TRUE)
    type_tab = type_tab[names(type_len)]
  } else if(resType == 'repeat'){
    dist_dat = data.frame(do.call("rbind", res$Repeat_distribution$Decided), stringsAsFactors = FALSE)
    dist_dat2 = data.frame(do.call("rbind", res$Micro_distribution$Decided), stringsAsFactors = FALSE)
    dist_dat = dist_dat[-which(is.na(dist_dat$query)),]; dist_dat2 = dist_dat2[-which(is.na(dist_dat2$query)),]
    if(includeUndecided){
      dist_dat_u = data.frame(do.call("rbind", res$Repeat_distribution$Undecided), stringsAsFactors = FALSE)
      dist_dat_u2 = data.frame(do.call("rbind", res$Micro_distribution$Undecided), stringsAsFactors = FALSE)
      dist_dat_u = dist_dat_u[-which(is.na(dist_dat_u$query)),]; dist_dat_u2 = dist_dat_u2[-which(is.na(dist_dat_u2$query)),]
      dist_dat = data.frame(rbind(dist_dat, dist_dat_u), stringsAsFactors = FALSE)
      dist_dat2 = data.frame(rbind(dist_dat2, dist_dat_u2), stringsAsFactors = FALSE)
    } else {}
    types = unique(dist_dat$repClass)
    type_tab = lapply(types, function(a){subset(dist_dat, dist_dat$repClass == a)})
    names(type_tab) = types
    type_len = sort(unlist(lapply(type_tab, nrow)), decreasing = TRUE)
    type_tab = type_tab[names(type_len)]
  } else {
    dist_dat = data.frame(do.call("rbind", res[[2]]$Decided), stringsAsFactors = FALSE)
    dist_dat = dist_dat[-which(is.na(dist_dat$query)),]
    if(includeUndecided){
      dist_dat_u = data.frame(do.call("rbind", res[[2]]$Undecided), stringsAsFactors = FALSE)
      dist_dat_u = dist_dat_u[-which(is.na(dist_dat_u$query)),]
      dist_dat = data.frame(rbind(dist_dat, dist_dat_u), stringsAsFactors = FALSE)
    }
  }
  message('- OK!')
  message('- Edit random set data.')
  if(!is.null(which(names(res) == 'Random_distribution'))){
    if(length(res$Random_distribution[[1]]) == 2){
      dist_dat_ran1 = lapply(res$Random_distribution, function(a){return(a[[1]])})
      dist_dat_ran2 = lapply(res$Random_distribution, function(a){return(a[[2]])})
      ran1 = data.frame(do.call("rbind", dist_dat_ran1), stringsAsFactors = FALSE)
      ran2 = data.frame(do.call("rbind", dist_dat_ran2), stringsAsFactors = FALSE)
      ran1 = ran1[-which(is.na(ran1$query)),]; ran2 = ran2[-which(is.na(ran2$query)),]
      if(resType == 'gene'){
        ran1_type = unique(ran1$gene_type)
        ran1_type_tab = lapply(ran1_type, function(a){subset(ran1, ran1$gene_type == a)})
      } else if(resType == 'repeat'){
        ran1_type = unique(ran1$repClass)
        ran1_type_tab = lapply(ran1_type, function(a){subset(ran1, ran1$repClass == a)})
      } else {}
      names(ran1_type_tab) = ran1_type
      ran1_type_len = sort(unlist(lapply(ran1_type_tab, nrow)), decreasing = TRUE)
      ran1_type_tab = ran1_type_tab[names(ran1_type_len)]
    } else if(length(res$Random_distribution[[1]]) == 1){
      ran1 = data.frame(do.call('rbind', res$Random_distribution), stringsAsFactors = FALSE)
    } else {dist_dat_ran = NULL}
  } else {dist_dat_ran = NULL}
  message('- OK!')
  message('- Write result file(s).')
  # create excel file object
  docs = openxlsx::createWorkbook()
  docs_random = openxlsx::createWorkbook()
  # cell style
  hs = openxlsx::createStyle(fgFill = "#228B22", halign = "CENTER", textDecoration = "Bold", fontColour = "white")
  hs_random = openxlsx::createStyle(fgFill = "#4B0082", halign = "CENTER", textDecoration = "Bold", fontColour = "white")
  # add sheets
  # write data tables
  if(resType == 'gene'){
    # add sheets
    openxlsx::addWorksheet(wb = docs, sheetName = "Inside_gene")
    openxlsx::addWorksheet(wb = docs_random, sheetName = "Inside_gene")
    lapply(names(type_len), function(a){openxlsx::addWorksheet(wb = docs, sheetName = a);
      return(cat('+ Add sheet to file : \"', a, '\"\n'))})
    lapply(names(ran1_type_len), function(a){openxlsx::addWorksheet(wb = docs_random, sheetName = a);
      return(cat('+ Add sheet to file : \"', a, '\"\n'))})
    openxlsx::addWorksheet(wb = docs, sheetName = "Dist_TSS")
    openxlsx::addWorksheet(wb = docs_random, sheetName = "Dist_TSS")
    # add data
    openxlsx::writeDataTable(wb = docs, sheet = "Inside_gene", x = res$Gene_inside, headerStyle = hs, rowNames = FALSE, colNames = TRUE, sep = '\t')
    openxlsx::writeDataTable(wb = docs_random, sheet = "Inside_gene", x = res$Gene_inside_ran, headerStyle = hs_random, rowNames = FALSE, colNames = TRUE, sep = '\t')
    lapply(names(type_len), function(a){openxlsx::writeDataTable(wb = docs, sheet = a, x = type_tab[[a]], headerStyle = hs, rowNames = FALSE, colName = TRUE, sep = '\t')})
    lapply(names(ran1_type_len), function(a){openxlsx::writeDataTable(wb = docs_random, sheet = a, x = ran1_type_tab[[a]], headerStyle = hs_random, rowNames = FALSE, colName = TRUE, sep = '\t')})
    openxlsx::writeDataTable(wb = docs, sheet = "Dist_TSS", x = dist_dat2, headerStyle = hs, rowNames = FALSE, colNames = TRUE, sep = '\t')
    openxlsx::writeDataTable(wb = docs_random, sheet = "Dist_TSS", x = ran2, headerStyle = hs_random, rowNames = FALSE, colNames = TRUE, sep = '\t')
  } else if(resType == 'repeat'){
    # add sheets
    openxlsx::addWorksheet(wb = docs, sheetName = "Inside_repeat")
    openxlsx::addWorksheet(wb = docs_random, sheetName = "Inside_repeat")
    lapply(names(type_len), function(a){openxlsx::addWorksheet(wb = docs, sheetName = a);
      return(cat('+ Add sheet to file : \"', a, '\"\n'))})
    lapply(names(ran1_type_len), function(a){openxlsx::addWorksheet(wb = docs_random, sheetName = a);
      return(cat('+ Add sheet to file : \"', a, '\"\n'))})
    openxlsx::addWorksheet(wb = docs, sheetName = "Inside_microsatellite")
    openxlsx::addWorksheet(wb = docs_random, sheetName = "Inside_microsatellite")
    openxlsx::addWorksheet(wb = docs, sheetName = "Dist_microsatellite")
    openxlsx::addWorksheet(wb = docs_random, sheetName = "Dist_microsatellite")
    # add data
    openxlsx::writeDataTable(wb = docs, sheet = "Inside_repeat", x = res$Repeat_inside, headerStyle = hs, rowNames = FALSE, colNames = TRUE, sep = '\t')
    openxlsx::writeDataTable(wb = docs, sheet = "Inside_repeat", x = res$Repeat_inside_ran, headerStyle = hs_random, rowNames = FALSE, colNames = TRUE, sep = '\t')
    lapply(names(type_len), function(a){openxlsx::writeDataTable(wb = docs, sheet = a, x = type_tab[[a]], headerStyle = hs, rowNames = FALSE, colName = TRUE, sep = '\t')})
    lapply(names(ran1_type_len), function(a){openxlsx::writeDataTable(wb = docs_random, sheet = a, x = ran1_type_tab[[a]], headerStyle = hs_random, rowNames = FALSE, colName = TRUE, sep = '\t')})
    openxlsx::writeDataTable(wb = docs, sheet = "Inside_microsatellite", x = res$Micro_inside, headerStyle = hs, rowNames = FALSE, colNames = TRUE, sep = '\t')
    openxlsx::writeDataTable(wb = docs_random, sheet = "Inside_microsatellite", x = res$Micro_inside_ran, headerStyle = hs_random, rowNames = FALSE, colNames = TRUE, sep = '\t')
    openxlsx::writeDataTable(wb = docs, sheet = "Dist_microsatellite", x = dist_dat2, headerStyle = hs, rowNames = FALSE, colNames = TRUE, sep = '\t')
    openxlsx::writeDataTable(wb = docs_random, sheet = "Dist_microsatellite", x = ran2, headerStyle = hs_random, rowNames = FALSE, colNames = TRUE, sep = '\t')
  } else {
    # add sheet
    openxlsx::addWorksheet(wb = docs, sheetName = paste0("Inside_", resType))
    openxlsx::addWorksheet(wb = docs_random, sheetName = paste0("Inside_", resType))
    openxlsx::addWorksheet(wb = docs, sheetName = paste0("Dist_", resType))
    openxlsx::addWorksheet(wb = docs_random, sheetName = paste0("Dist_", resType))
    # add data
    openxlsx::writeDataTable(wb = docs, sheet = paste0("Inside_", resType), x = res[[1]], headerStyle = hs, rowNames = FALSE, colNames = TRUE, sep = '\t')
    openxlsx::writeDataTable(wb = docs_random, sheet = paste0("Inside_", resType), x = res[[6]], headerStyle = hs_random, rowNames = FALSE, colNames = TRUE, sep = '\t')
    openxlsx::writeDataTable(wb = docs, sheet = paste0("Dist_", resType), x = dist_dat, headerStyle = hs, rowNames = FALSE, colNames = TRUE, sep = '\t')
    openxlsx::writeDataTable(wb = docs_random, sheet = paste0("Dist_", resType), x = ran1, headerStyle = hs_random, rowNames = FALSE, colNames = TRUE, sep = '\t')
  }
  # write file
  openxlsx::saveWorkbook(wb = docs, file = paste0(outPath, '/', outFileName, '_exprimental.xlsx'), overwrite = FALSE)
  openxlsx::saveWorkbook(wb = docs_random, file = paste0(outPath, '/', outFileName, '_random.xlsx'), overwrite = FALSE)
  message('- OK!')
  message('- Test data set by the feature type.')
  if(resType %in% c('gene', 'repeat')){
    obs = type_tab[names(type_len)]
    ran = ran1_type_tab[names(type_len)]
    hist_pal = unique(c(RColorBrewer::brewer.pal(11, 'Spectral')[1:5], 'yellow', RColorBrewer::brewer.pal(11, 'Spectral')[7:11],
                        RColorBrewer::brewer.pal(8, 'Dark2')))
    obs_hist = lapply(obs, function(a){hist(a$dist, breaks = ranges, plot = FALSE)})
    grDevices::png(paste0(outPath, '/', outFileName, '_histogram_observed.png'), width = 1200, height = 1200)
    par(mfrow = c(4, ceiling(length(obs)/4)))
    lapply(c(1:length(obs_hist)), function(a){plot(obs_hist[[a]], col = hist_pal[a], main = names(type_len)[a], cex.main = 3)})
    grDevices::dev.off()
    ran_hist = lapply(ran, function(a){hist(a$dist, breaks = ranges, plot = FALSE)})
    grDevices::png(paste0(outPath, '/', outFileName, '_histogram_random.png'), width = 1200, height = 1200)
    par(mfrow = c(4, ceiling(length(ran)/4)))
    lapply(c(1:length(ran_hist)), function(a){plot(ran_hist[[a]], col = hist_pal[a], main = names(type_len)[a], cex.main = 3)})
    grDevices::dev.off()
    
    test_res = lapply(c(1:length(obs_hist)), function(a){
      lapply(c(1:length(obs_hist[[a]]$counts)), function(b){
        suppressWarnings(stats::prop.test(x = c(obs_hist[[a]]$counts[b], ran_hist[[a]]$counts[b]),
                                          n = c(sum(obs_hist[[a]]$counts), sum(ran_hist[[a]]$counts)),
                                          alternative = "two.sided", correct = FALSE))
      })
    })
    test_p_value = lapply(test_res, function(a){unlist(lapply(a, function(b){b$p.value}))})
    gg_tab_test = lapply(c(1:length(obs_hist)), function(a){
      obs_count = obs_hist[[a]]$counts; ran_count = ran_hist[[a]]$counts
      obs_num = obs_hist[[a]]$counts/sum(obs_hist[[a]]$counts); ran_num = ran_hist[[a]]$counts/sum(ran_hist[[a]]$counts)
      group = paste0('I', c(1:length(obs_num)))
      p = test_p_value[[a]]
      convert_p = -log10(test_p_value[[a]]+0.00001)*3
      type = names(type_len)[a]
      return(data.frame(obs_count, ran_count, obs_num, ran_num, group, p, convert_p, type))
    })
    gg_tab_test_all = do.call("rbind", gg_tab_test)
    set_pal = hist_pal[1:length(type_len)]; names(set_pal) = names(type_len)
    valid_p = rep('No', nrow(gg_tab_test_all))
    valid_p[which(gg_tab_test_all$p <= 0.05)] = 'Significant difference'
    gg_tab_test_all = cbind(gg_tab_test_all, valid_p)
    grDevices::png(paste0(outPath, '/', outFileName, '_pvalue_plot.png'), width = 800, height = 600)
    p_plot = ggplot2::ggplot(gg_tab_test_all, 
                    ggplot2::aes(x = group, y = factor(type, levels = names(type_len)[c(length(type_len):1)]), 
                                 size = convert_p, color = type, shape = factor(valid_p, levels = unique(valid_p)[c(2,1)]))) +
             ggplot2::geom_point() + 
             ggplot2::scale_size(range = c(min(gg_tab_test_all$convert_p[-which(is.na(gg_tab_test_all$convert_p))]), 
                                           max(gg_tab_test_all$convert_p[-which(is.na(gg_tab_test_all$convert_p))])), guide = FALSE) +
             ggplot2::scale_color_manual(values = set_pal, guide = FALSE) +
             ggplot2::theme(panel.background = ggplot2::element_blank(),
                            panel.grid = ggplot2::element_line(linetype = "dotted", color = "black"),
                            axis.text.x = ggplot2::element_text(hjust = 1, angle = 45, size = 15),
                            axis.text.y = ggplot2::element_text(size = 15),
                            axis.title.y = ggplot2::element_blank(),
                            legend.title = ggplot2::element_blank(),
                            legend.key = ggplot2::element_blank(),
                            legend.text = ggplot2::element_text(size = 15)) +
             ggplot2::scale_x_discrete(labels = ranges[-which(ranges == 0)]/1000) +
             ggplot2::scale_shape_manual(values = c(16, 1)) + ggplot2::labs(x = "Distance") +
             ggplot2::guides(shape = ggplot2::guide_legend(override.aes = list(size = 5)))
    suppressWarnings(print(p_plot))
    grDevices::dev.off()} else {cat('[NOTICE] Do not run this part by CpG/Variant data.\n')}
  message('- OK!')
  message('----- Finish. (Time : ', date(), ')')
  return(gg_tab_test_all)
}
