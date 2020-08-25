#' @title Make the result object and document.
#' 
#' @description
#' Rearrange the result from annotation functions.
#' 
#' @usage 
#' makeDocument(res, dataType, excelOut = TRUE, 
#'              includeUndecided = FALSE, outPath = getwd(),
#'              outFileName = paste0('RIPAT', round(unclass(Sys.time()))))
#' 
#' @param res a GR object. This object is output of \code{annoByGene, annoByCpG, annoByRepeat, annoByVar} function.
#' @param dataType a character vector. User enter the annotation type of input
#'                 such as gene, cpg, repeat and variant.
#' @param excelOut TRUE or FALSE. If user want to make excel file, enter TRUE.
#'                 Default is TRUE.
#' @param includeUndecided TRUE or FALSE. If user want to use undecided hits in analysis, enter TRUE.
#'                         Default is FALSE.
#' @param outPath an string vector. Plots are saved in this path. Default value is R home directory.
#' @param outFileName a character vector. Attached ID to the result file name.
#' 
#' @return Make output table and excel files about vector integration sites and proportion test result.
#' 
#' @examples 
#' data(blast_gene)
#' makeDocument(res = blast_gene, dataType = 'gene', outFileName = 'blast_gene_res')
#'
#' @export
makeDocument = function(res, dataType ='gene', excelOut = TRUE, includeUndecided = FALSE, outPath = getwd(), outFileName = paste0('RIPAT', round(unclass(Sys.time())))){
  if(!is.null(res)){
    message('----- Make result documents. (Time : ', date(), ')')
    message('- Validate options')
    if(!is.logical(excelOut)){stop("[ERROR] 'excelOut' should be logical constants.\n----- This process is halted. (Time : ", date(), ")\n")} else{}
    if(stringr::str_ends(outPath, pattern = '/')){
      outPath = stringr::word(outPath, start = 1, end = nchar(outPath), sep = '')
    } else {}
    message('- Edit experimental data')
    # Edit result
    inside_tab = data.frame(res[[1]], stringsAsFactors = FALSE)
    dist_tab = data.frame(do.call('rbind', res[[2]]$Decided), stringsAsFactors = FALSE)
    if(length(which(is.na(dist_tab$query))) != 0){dist_tab = dist_tab[-which(is.na(dist_tab$query)),]}
    if(dataType == 'repeat'){
      inside_tab2 = data.frame(res[[5]], stringsAsFactors = FALSE)
      dist_tab2 = data.frame(do.call('rbind', res[[6]]$Decided), stringsAsFactors = FALSE)
      if(length(which(is.na(dist_tab2$query))) != 0){dist_tab2 = dist_tab2[-which(is.na(dist_tab2$query)),]}
      exp_types = unique(dist_tab$repClass)
      type_tab1 = lapply(exp_types, function(a){subset(dist_tab, dist_tab$repClass == a)})
      names(type_tab1) = exp_types
    } else if(dataType == 'gene'){
      dist_tab2 = data.frame(do.call('rbind', res[[5]]$Decided), stringsAsFactors = FALSE)
      if(length(which(is.na(dist_tab2$query))) != 0){dist_tab2 = dist_tab2[-which(is.na(dist_tab2$query)),]}
      exp_types = unique(dist_tab$gene_type); exp_types2 = unique(dist_tab2$gene_type)
      type_tab1 = lapply(exp_types, function(a){subset(dist_tab, dist_tab$gene_type == a)})
      type_tab2 = lapply(exp_types2, function(a){subset(dist_tab2, dist_tab2$gene_type == a)})
      names(type_tab1) = exp_types; names(type_tab2) = exp_types2
    }
    if(includeUndecided){
      inside_u_tab = data.frame(do.call('rbind', res[[1]]$Undecided), stringsAsFactors = FALSE)
      dist_u_tab = data.frame(do.call('rbind', res[[2]]$Undecided), stringsAsFactors = FALSE)
      if(dataType == 'repeat'){
        inside_u_tab2 = data.frame(do.call('rbind', res[[5]]$Undecided), stringsAsFactors = FALSE)
        dist_u_tab2 = data.frame(do.call('rbind', res[[6]]$Undecided), stringsAsFactors = FALSE)
        dist_u_tab2 = dist_u_tab2[-which(is.na(dist_u_tab2$query)),]
      } else if(dataType == 'gene'){
        dist_u_tab2 = data.frame(do.call('rbind', res[[5]]$Undecided), stringsAsFactors = FALSE)
        dist_u_tab2 = dist_u_tab2[-which(is.na(dist_u_tab2$query)),]
      }}
    message('- Edit random data')
    if(length(which(names(res) == 'Random_distribution')) != 0){
      if(dataType == 'repeat'){
        rinside_tab = data.frame(do.call('rbind', res[[10]]), stringsAsFactors = FALSE)
        rdist_list = lapply(res[[12]], function(a){return(a[[1]])})
        rdist_list2 = lapply(res[[12]], function(a){return(a[[2]])})
        rdist_tab = data.frame(do.call('rbind', rdist_list), stringsAsFactors = FALSE)
        rdist_tab2 = data.frame(do.call('rbind', rdist_list2), stringsAsFactors = FALSE)
        if(length(which(is.na(rdist_tab$query))) != 0){rdist_tab = rdist_tab[-which(is.na(rdist_tab$query)),]}
        if(length(which(is.na(rdist_tab2$query))) != 0){rdist_tab2 = rdist_tab2[-which(is.na(rdist_tab2$query)),]}
        ran_types = unique(rdist_tab$repClass)
        ran_type_tab = lapply(ran_types, function(a){subset(rdist_tab, rdist_tab$repClass == a)})
        names(ran_type_tab) = ran_types
      } else if(dataType == 'gene'){
        rinside_tab = data.frame(res[[9]], stringsAsFactors = FALSE)
        rdist_list = lapply(res[[10]], function(a){return(a[[1]])})
        rdist_list2 = lapply(res[[10]], function(a){return(a[[2]])})
        rdist_tab = data.frame(do.call('rbind', rdist_list), stringsAsFactors = FALSE)
        rdist_tab2 = data.frame(do.call('rbind', rdist_list2), stringsAsFactors = FALSE)
        if(length(which(is.na(rdist_tab$query))) != 0){rdist_tab = rdist_tab[-which(is.na(rdist_tab$query)),]}
        if(length(which(is.na(rdist_tab2$query))) != 0){rdist_tab2 = rdist_tab2[-which(is.na(rdist_tab2$query)),]}
        ran_types = unique(rdist_tab$gene_type); ran_types2 = unique(rdist_tab2$gene_type)
        ran_type_tab = lapply(ran_types, function(a){subset(rdist_tab, rdist_tab$gene_type == a)})
        ran_type_tab2 = lapply(ran_types2, function(a){subset(rdist_tab2, rdist_tab2$gene_type == a)})
        names(ran_type_tab) = ran_types; names(ran_type_tab2) = ran_types2
      } else if(dataType %in% c('cpg', 'variant')){
        rinside_tab = data.frame(res[[6]], stringsAsFactors = FALSE)
        rdist_tab = data.frame(do.call('rbind', res[[7]]), stringsAsFactors = FALSE)
        if(length(which(is.na(rdist_tab$query))) != 0){rdist_tab = rdist_tab[-which(is.na(rdist_tab$query)),]}
      } else {stop("[ERROR] 'dataType' needs the data type used in annotation.\n----- This process is halted. (Time : ", date(), ")\n")}
    } else {}
    message('- Compare experimental data and random data')
    if(dataType %in% c('gene', 'repeat')){
      range = as.numeric(as.character(res[[3]][,1])) * 1000
      ranges = seq(min(range), max(range), by = range[2]-range[1])
      obs_hist1 = lapply(type_tab1, function(a){hist(a$dist, breaks = ranges, plot = FALSE)})
      if(dataType == 'gene'){
        obs_hist2 = lapply(type_tab2, function(a){hist(a$dist, breaks = ranges, plot = FALSE)})
        hist_file1 = '_histogram_observed_gene.png'; hist_file2 = '_histogram_observed_tss.png'
        hist_file3 = '_pvalue_plot_gene.png'; hist_file4 = '_pvalue_plot_tss.png'
        obs1 = res$Gene_plot_data[which(res$Gene_plot_data$Group == 'Observed'),]; ran1 = res$Gene_plot_data[which(res$Gene_plot_data$Group == 'Random'),]
        obs2 = res$TSS_plot_data[which(res$TSS_plot_data$Group == 'Observed'),]; ran2 = res$TSS_plot_data[which(res$TSS_plot_data$Group == 'Random'),]
        if(length(which(names(res) == 'Random_distribution')) != 0){
          ran_hist2 = lapply(ran_type_tab2, function(a){hist(a$dist, breaks = ranges, plot = FALSE)})
          obs_all1 = obs1$Count[1]/obs1$Freq[1]; ran_all1 = ran1$Count[1]/ran1$Freq[1]
          test_all_res1 = lapply(c(1:nrow(obs1)), function(a){suppressWarnings(stats::prop.test(x = c(obs1$Count[a], ran1$Count[a]), n = c(obs_all1, ran_all1), alternative = 'two.sided', correct = FALSE))})
          test_all_res2 = lapply(c(1:nrow(obs2)), function(a){suppressWarnings(stats::prop.test(x = c(obs2$Count[a], ran2$Count[a]), n = c(obs_all1, ran_all1), alternative = 'two.sided', correct = FALSE))})
          test_all_p_value1 = unlist(lapply(test_all_res1, function(a){a$p.value})); test_all_p_value2 = unlist(lapply(test_all_res2, function(a){a$p.value}))
          htab = res$Gene_plot_data; htab2 = res$TSS_plot_data
          htab = data.frame(htab[c(1:(nrow(htab)/2)),], htab[c((nrow(htab)/2+1):nrow(htab)),c(2:4)], stringsAsFactors = FALSE)
          htab2 = data.frame(htab2[c(1:(nrow(htab2)/2)),], htab2[c((nrow(htab2)/2+1):nrow(htab2)),c(2:4)], stringsAsFactors = FALSE)
          colnames(htab) = c('Range', paste0(colnames(htab)[2:4], rep(c(1, 2), c(3,3))))
          colnames(htab2) = c('Range', paste0(colnames(htab2)[2:4], rep(c(1, 2), c(3,3))))
          htab_new1 = data.frame(htab, 'test_pvalue' = test_all_p_value1); htab_new2 = data.frame(htab2, 'test_pvalue' = test_all_p_value2)
          names(htab_new1) = names(htab_new2) = c('Range', 'Group1', 'Count1', 'Freq1', 'Group2', 'Count2', 'Freq2',  'test_pvalue')
          test_res2 = lapply(c(1:length(obs_hist2)), function(a){
            lapply(c(1:length(obs_hist2[[a]]$counts)), function(b){
              suppressWarnings(stats::prop.test(x = c(obs_hist2[[a]]$counts[b], ran_hist2[[a]]$counts[b]),
                                                n = c(sum(obs_hist2[[a]]$counts), sum(ran_hist2[[a]]$counts)),
                                                alternative = "two.sided", correct = FALSE))})})
          test_p_value2 = lapply(test_res2, function(a){unlist(lapply(a, function(b){b$p.value}))})
          gg_tab_test2 = lapply(c(1:length(obs_hist2)), function(a){
            obs_count = obs_hist2[[a]]$counts; ran_count = ran_hist2[[a]]$counts
            obs_num = obs_hist2[[a]]$counts/sum(obs_hist2[[a]]$counts); ran_num = ran_hist2[[a]]$counts/sum(ran_hist2[[a]]$counts)
            group = paste0('I', c(1:(length(ranges)-1)))
            p = test_p_value2[[a]]
            convert_p = -log10(test_p_value2[[a]]+0.00001)*3
            type = names(type_tab2)[a]
            return(data.frame(obs_count, ran_count, obs_num, ran_num, group, p, convert_p, type))
          })
          gg_tab_test_all2 = do.call("rbind", gg_tab_test2); valid_p2 = rep('No', nrow(gg_tab_test_all2))
          valid_p2[which(gg_tab_test_all2$p <= 0.05)] = 'Significant difference'
          gg_tab_test_all2 = cbind(gg_tab_test_all2, valid_p2)
          gg_tab_test_all2$group = as.character(rep(ranges[-which(ranges == 0)]/1000, nrow(gg_tab_test_all2)/(length(ranges)-1)))
          gg_tab_test_all2$group = factor(gg_tab_test_all2$group, levels = as.character(ranges[-which(ranges == 0)]/1000))
          grDevices::png(paste0(outPath, '/', outFileName, hist_file4), width = 1200, height = 800)
          p_plot2 = ggplot2::ggplot(gg_tab_test_all2, 
                                    ggplot2::aes(x = group, y = factor(type, levels = names(obs_hist2)[c(length(obs_hist2):1)]), 
                                                 size = convert_p, shape = factor(valid_p2, levels = unique(valid_p2)[c(2,1)]))) +
            ggplot2::geom_point(color = 'navy') + 
            ggplot2::scale_size(range = c(min(gg_tab_test_all2$convert_p[-which(is.na(gg_tab_test_all2$convert_p))]), 
                                          max(gg_tab_test_all2$convert_p[-which(is.na(gg_tab_test_all2$convert_p))])), guide = FALSE) +
            ggplot2::theme(panel.background = ggplot2::element_blank(), 
                           panel.grid = ggplot2::element_line(linetype = "dotted", color = "black"),
                           axis.text.x = ggplot2::element_text(hjust = 1, angle = 45, size = 15), axis.text.y = ggplot2::element_text(size = 15),
                           axis.title.y = ggplot2::element_blank(), 
                           legend.title = ggplot2::element_blank(), legend.key = ggplot2::element_blank(), legend.text = ggplot2::element_text(size = 15),
                           axis.line = ggplot2::element_line(color = "black")) +
            ggplot2::scale_x_discrete(labels = ranges[-which(ranges == 0)]/1000) +
            ggplot2::scale_shape_manual(values = c(16, 1)) + ggplot2::labs(x = "Distance") +
            ggplot2::guides(shape = ggplot2::guide_legend(override.aes = list(size = 5)))
          suppressWarnings(print(p_plot2))
          grDevices::dev.off()
        }
        grDevices::png(paste0(outPath, '/', outFileName, hist_file2), width = 3000, height = 2000)
        par(mfrow = c(4, ceiling(length(exp_types2)/4)))
        lapply(c(1:length(obs_hist2)), function(a){plot(obs_hist2[[a]], col = '#008080', main = names(type_tab2)[a], cex.main = 3)})
        grDevices::dev.off()
        grDevices::png(paste0(outPath, '/', outFileName, '_random', hist_file2), width = 3000, height = 2000)
        par(mfrow = c(6, ceiling(length(ran_hist2)/6)))
        lapply(c(1:length(ran_hist2)), function(a){plot(ran_hist2[[a]], col = '#4682B4', main = names(ran_type_tab2)[a], cex.main = 3)})
        grDevices::dev.off()
      } else {
        hist_file1 = '_histogram_observed_repeat.png'
        hist_file3 = '_pvalue_plot_repeat.png'
        obs1 = res$Repeat_plot_data[which(res$Repeat_plot_data$Group == 'Observed'),]; ran1 = res$Repeat_plot_data[which(res$Repeat_plot_data$Group == 'Random'),]
        obs2 = res$Micro_plot_data[which(res$Micro_plot_data$Group == 'Observed'),]; ran2 = res$Micro_plot_data[which(res$Micro_plot_data$Group == 'Random'),]
        if(length(which(names(res) == 'Random_distribution')) != 0){
          check_key = which(obs1$Count != 0)[1]; check_key_ran = which(ran1$Count != 0)[1]
          obs_all1 = obs1$Count[check_key]/obs1$Freq[check_key]; ran_all1 = ran1$Count[check_key_ran]/ran1$Freq[check_key_ran]
          test_all_res1 = lapply(c(1:nrow(obs1)), function(a){suppressWarnings(stats::prop.test(x = c(obs1$Count[a], ran1$Count[a]), n = c(obs_all1, ran_all1), alternative = 'two.sided', correct = FALSE))})
          test_all_res2 = lapply(c(1:nrow(obs2)), function(a){suppressWarnings(stats::prop.test(x = c(obs2$Count[a], ran2$Count[a]), n = c(obs_all1, ran_all1), alternative = 'two.sided', correct = FALSE))})
          test_all_p_value1 = unlist(lapply(test_all_res1, function(a){a$p.value}))
          test_all_p_value2 = unlist(lapply(test_all_res2, function(a){a$p.value}))
          htab = res$Repeat_plot_data; htab2 = res$Micro_plot_data
          htab = data.frame(htab[c(1:(nrow(htab)/2)),], htab[c((nrow(htab)/2+1):nrow(htab)),c(2:4)], stringsAsFactors = FALSE)
          htab2 = data.frame(htab2[c(1:(nrow(htab2)/2)),], htab2[c((nrow(htab2)/2+1):nrow(htab2)),c(2:4)], stringsAsFactors = FALSE)
          colnames(htab) = c('Range', paste0(colnames(htab)[2:4], rep(c(1,2), c(3,3))))
          colnames(htab2) = c('Range', paste0(colnames(htab2)[2:4], rep(c(1,2), c(3,3))))
          htab_new1 = data.frame(htab, 'test_pvalue' = test_all_p_value1); htab_new2 = data.frame(htab2, 'test_pvalue' = test_all_p_value2)
          names(htab_new1) = names(htab_new2) = c('Range', 'Group1', 'Count1', 'Freq1', 'Group2', 'Count2', 'Freq2',  'test_pvalue')
        }
      }
      grDevices::png(paste0(outPath, '/', outFileName, hist_file1), width = 2000, height = 1000)
      par(mfrow = c(3, ceiling(length(exp_types)/3)))
      lapply(c(1:length(obs_hist1)), function(a){plot(obs_hist1[[a]], col = '#008080', main = names(type_tab1)[a], cex.main = 3)})
      grDevices::dev.off()
      ran_hist1 = lapply(ran_type_tab, function(a){hist(a$dist, breaks = ranges, plot = FALSE)})
      grDevices::png(paste0(outPath, '/', outFileName, '_random', hist_file1), width = 2000, height = 1000)
      par(mfrow = c(4, ceiling(length(ran_hist1)/4)))
      lapply(c(1:length(ran_hist1)), function(a){plot(ran_hist1[[a]], col = '#4682B4', main = names(ran_type_tab)[a], cex.main = 3)})
      grDevices::dev.off()
      if(length(which(names(res) == 'Random_distribution')) != 0){
        test_res1 = lapply(c(1:length(obs_hist1)), function(a){
          lapply(c(1:length(obs_hist1[[a]]$counts)), function(b){
            suppressWarnings(stats::prop.test(x = c(obs_hist1[[a]]$counts[b], ran_hist1[[a]]$counts[b]),
                                              n = c(sum(obs_hist1[[a]]$counts), sum(ran_hist1[[a]]$counts)),
                                              alternative = "two.sided", correct = FALSE))})})
        test_p_value1 = lapply(test_res1, function(a){unlist(lapply(a, function(b){b$p.value}))})
        gg_tab_test1 = lapply(c(1:length(obs_hist1)), function(a){
          obs_count = obs_hist1[[a]]$counts; ran_count = ran_hist1[[a]]$counts
          obs_num = obs_hist1[[a]]$counts/sum(obs_hist1[[a]]$counts); ran_num = ran_hist1[[a]]$counts/sum(ran_hist1[[a]]$counts)
          group = paste0('I', c(1:(length(ranges)-1)))
          p = test_p_value1[[a]]
          convert_p = -log10(test_p_value1[[a]]+0.00001)*3
          type = names(type_tab1)[a]
          return(data.frame(obs_count, ran_count, obs_num, ran_num, group, p, convert_p, type))
        })
        gg_tab_test_all1 = do.call("rbind", gg_tab_test1)
        valid_p1 = rep('No', nrow(gg_tab_test_all1))
        valid_p1[which(gg_tab_test_all1$p <= 0.05)] = 'Significant difference'
        gg_tab_test_all1 = cbind(gg_tab_test_all1, valid_p1)
        gg_tab_test_all1$group = as.character(rep(ranges[-which(ranges == 0)]/1000, nrow(gg_tab_test_all1)/(length(ranges)-1)))
        gg_tab_test_all1$group = factor(gg_tab_test_all1$group, levels = as.character(ranges[-which(ranges == 0)]/1000))
        grDevices::png(paste0(outPath, '/', outFileName, hist_file3), width = 1000, height = 800)
        p_plot1 = ggplot2::ggplot(gg_tab_test_all1, 
                                  ggplot2::aes(x = group, y = factor(type, levels = names(obs_hist1)[c(length(obs_hist1):1)]), 
                                               size = convert_p, shape = factor(valid_p1, levels = unique(valid_p1)[c(2,1)]))) +
          ggplot2::geom_point(color = '#663399') + 
          ggplot2::scale_size(range = c(min(gg_tab_test_all1$convert_p[-which(is.na(gg_tab_test_all1$convert_p))]), 
                                        max(gg_tab_test_all1$convert_p[-which(is.na(gg_tab_test_all1$convert_p))])), guide = FALSE) +
          ggplot2::theme(panel.background = ggplot2::element_blank(), 
                         panel.grid = ggplot2::element_line(linetype = "dotted", color = "black"),
                         axis.line = ggplot2::element_line(color = "black"),
                         axis.text.x = ggplot2::element_text(hjust = 1, angle = 45, size = 15), 
                         axis.text.y = ggplot2::element_text(size = 15),
                         axis.title.y = ggplot2::element_blank(), legend.title = ggplot2::element_blank(), legend.key = ggplot2::element_blank(), legend.text = ggplot2::element_text(size = 15)) +
          ggplot2::scale_x_discrete(labels = ranges[-which(ranges == 0)]/1000) +
          ggplot2::scale_shape_manual(values = c(16, 1)) + ggplot2::labs(x = "Distance") +
          ggplot2::guides(shape = ggplot2::guide_legend(override.aes = list(size = 5)))
        suppressWarnings(print(p_plot1))
        grDevices::dev.off()
      }
    } else {
      if(dataType == 'cpg'){htab = res$CpG_plot_data} else if(dataType == 'variant'){htab = res$Variant_plot_data}
      obs1 = htab[which(htab$Group == 'Observed'),]; ran1 = htab[which(htab$Group == 'Random'),]
      obs_all1 = obs1$Count[1]/obs1$Freq[1]; ran_all1 = ran1$Count[1]/ran1$Freq[1]
      test_all_res1 = lapply(c(1:nrow(obs1)), function(a){suppressWarnings(stats::prop.test(x = c(obs1$Count[a], ran1$Count[a]), n = c(obs_all1, ran_all1), alternative = 'two.sided', correct = FALSE))})
      test_all_p_value1 = unlist(lapply(test_all_res1, function(a){a$p.value}))
      htab = data.frame(htab[c(1:(nrow(htab)/2)),], htab[c((nrow(htab)/2+1):nrow(htab)),c(2:4)], stringsAsFactors = FALSE)
      htab_new = data.frame(htab, 'test_pvalue' = test_all_p_value1)
      names(htab_new) = c('Range', 'Group1', 'Count1', 'Freq1', 'Group2', 'Count2', 'Freq2',  'test_pvalue')
    }
    message('- OK!')
    if(excelOut){
      message('- Write result file(s).')
      docs = openxlsx::createWorkbook(); docs2 = openxlsx::createWorkbook()
      test_docs = openxlsx::createWorkbook()
      hs = openxlsx::createStyle(fgFill = "forestgreen", halign = "CENTER", textDecoration = "Bold", fontColour = "white")
      hs_random = openxlsx::createStyle(fgFill = "mediumpurple", halign = "CENTER", textDecoration = "Bold", fontColour = "white")
      hs_test = openxlsx::createStyle(fgFill = "dodgerblue", halign = "CENTER", textDecoration = "Bold", fontColour = "white")
      if(dataType == 'gene'){
        cat('+ Add sheet to file : \"Exp_inside_gene\"\n'); openxlsx::addWorksheet(wb = docs, sheetName = "Exp_inside_gene")
        lapply(names(type_tab1), function(a){openxlsx::addWorksheet(wb = docs, sheetName = paste0('Exp_', stringr::str_sub(a, start = 1, end = 26)));
          return(cat(paste0('+ Add sheet to file : \"Exp_', a, '\"\n')))})
        cat('+ Add sheet to file : \"Ran_inside_gene\"\n'); openxlsx::addWorksheet(wb = docs, sheetName = "Ran_inside_gene")
        lapply(names(ran_type_tab), function(a){openxlsx::addWorksheet(wb = docs, sheetName = paste0('Ran_', stringr::str_sub(a, start = 1, end = 26)));
          return(cat(paste0('+ Add sheet to file : \"Ran_', a, '\"\n')))})
        lapply(names(type_tab2), function(a){openxlsx::addWorksheet(wb = docs2, sheetName = paste0('Exp_', stringr::str_sub(a, start = 1, end = 26)));
          return(cat(paste0('+ Add sheet to file : \"Exp_', a, '\"\n')))})
        lapply(names(ran_type_tab2), function(a){openxlsx::addWorksheet(wb = docs2, sheetName = paste0('Ran_', stringr::str_sub(a, start = 1, end = 26)));
          return(cat(paste0('+ Add sheet to file : \"Ran_', a, '\"\n')))})
        openxlsx::writeDataTable(wb = docs, sheet = "Exp_inside_gene", x = inside_tab, headerStyle = hs, rowNames = FALSE, colNames = TRUE, sep = '\t')
        openxlsx::writeDataTable(wb = docs, sheet = "Ran_inside_gene", x = rinside_tab, headerStyle = hs_random, rowNames = FALSE, colNames = TRUE, sep = '\t')
        lapply(names(type_tab1), function(a){openxlsx::writeDataTable(wb = docs, sheet = paste0('Exp_', stringr::str_sub(a, start = 1, end = 26)), x = type_tab1[[a]], headerStyle = hs, rowNames = FALSE, colName = TRUE, sep = '\t')})
        lapply(names(ran_type_tab), function(a){openxlsx::writeDataTable(wb = docs, sheet = paste0('Ran_', stringr::str_sub(a, start = 1, end = 26)), x = ran_type_tab[[a]], headerStyle = hs_random, rowNames = FALSE, colName = TRUE, sep = '\t')})
        lapply(names(type_tab2), function(a){openxlsx::writeDataTable(wb = docs2, sheet = paste0('Exp_', stringr::str_sub(a, start = 1, end = 26)), x = type_tab2[[a]], headerStyle = hs, rowNames = FALSE, colName = TRUE, sep = '\t')})
        lapply(names(ran_type_tab2), function(a){openxlsx::writeDataTable(wb = docs2, sheet = paste0('Ran_', stringr::str_sub(a, start = 1, end = 26)), x = ran_type_tab2[[a]], headerStyle = hs_random, rowNames = FALSE, colName = TRUE, sep = '\t')})
        cat(paste0('+ Add sheet to file : \"Test_all_type_Gene\"\n')); openxlsx::addWorksheet(wb = test_docs, sheetName = "Test_all_type_Gene")
        cat(paste0('+ Add sheet to file : \"Test_all_type_TSS\"\n')); openxlsx::addWorksheet(wb = test_docs, sheetName = "Test_all_type_TSS")
        cat(paste0('+ Add sheet to file : \"Test_by_types_Gene\"\n')); openxlsx::addWorksheet(wb = test_docs, sheetName = "Test_by_types_Gene")
        cat(paste0('+ Add sheet to file : \"Test_by_types_TSS\"\n')); openxlsx::addWorksheet(wb = test_docs, sheetName = "Test_by_types_TSS")
        openxlsx::writeDataTable(wb = test_docs, sheet = "Test_all_type_Gene", x = htab_new1, headerStyle = hs_test, rowNames = FALSE, colNames = TRUE, sep = '\t')
        openxlsx::writeDataTable(wb = test_docs, sheet = "Test_all_type_TSS", x = htab_new2, headerStyle = hs_test, rowNames = FALSE, colNames = TRUE, sep = '\t')
        openxlsx::writeDataTable(wb = test_docs, sheet = "Test_by_types_Gene", x = gg_tab_test_all1, headerStyle = hs_test, rowNames = FALSE, colNames = TRUE, sep = '\t')
        openxlsx::writeDataTable(wb = test_docs, sheet = "Test_by_types_TSS", x = gg_tab_test_all2, headerStyle = hs_test, rowNames = FALSE, colNames = TRUE, sep = '\t')
        openxlsx::saveWorkbook(wb = docs, file = paste0(outPath, '/', outFileName, '_Gene.xlsx'), overwrite = TRUE)
        openxlsx::saveWorkbook(wb = docs2, file = paste0(outPath, '/', outFileName, '_TSS.xlsx'), overwrite = TRUE)
        openxlsx::saveWorkbook(wb = test_docs, file = paste0(outPath, '/', outFileName, '_test_result.xlsx'), overwrite = TRUE)
        res_list = list('Test_result_gene' = htab_new1, 'Test_result_TSS' = htab_new2, 'Test_result_by_gene_types' = gg_tab_test_all1, 'Test_result_by_TSS_types' = gg_tab_test_all2)
      } else if(dataType == 'repeat'){
        cat(paste0('+ Add sheet to file : \"Exp_inside_repeat\"\n')); openxlsx::addWorksheet(wb = docs, sheetName = "Exp_inside_repeat")
        lapply(names(type_tab1), function(a){openxlsx::addWorksheet(wb = docs, sheetName = paste0('Exp_', stringr::str_sub(a, start = 1, end = 26)));
          return(cat(paste0('+ Add sheet to file : \"Exp_', a, '\"\n')))})
        cat(paste0('+ Add sheet to file : \"Ran_inside_repeat\"\n')); openxlsx::addWorksheet(wb = docs, sheetName = "Ran_inside_repeat")
        lapply(names(ran_type_tab), function(a){openxlsx::addWorksheet(wb = docs, sheetName = paste0('Ran_', stringr::str_sub(a, start = 1, end = 26)));
          return(cat(paste0('+ Add sheet to file : \"Ran_', a, '\"\n')))})
        cat(paste0('+ Add sheet to file : \"Exp_inside_microsatellite\"\n')); openxlsx::addWorksheet(wb = docs2, sheetName = "Exp_inside_microsatellite")
        cat(paste0('+ Add sheet to file : \"Exp_dist_microsatellite\"\n')); openxlsx::addWorksheet(wb = docs2, sheetName = "Exp_dist_microsatellite")
        cat(paste0('+ Add sheet to file : \"Ran_inside_microsatellite\"\n')); openxlsx::addWorksheet(wb = docs2, sheetName = "Ran_inside_microsatellite")
        cat(paste0('+ Add sheet to file : \"Ran_dist_microsatellite\"\n')); openxlsx::addWorksheet(wb = docs2, sheetName = "Ran_dist_microsatellite")
        openxlsx::writeDataTable(wb = docs, sheet = "Exp_inside_repeat", x = res$Repeat_inside, headerStyle = hs, rowNames = FALSE, colNames = TRUE, sep = '\t')
        openxlsx::writeDataTable(wb = docs, sheet = "Ran_inside_repeat", x = res$Repeat_inside_ran, headerStyle = hs_random, rowNames = FALSE, colNames = TRUE, sep = '\t')
        lapply(names(type_tab1), function(a){openxlsx::writeDataTable(wb = docs, sheet = paste0('Exp_', a), x = type_tab1[[a]], headerStyle = hs, rowNames = FALSE, colName = TRUE, sep = '\t')})
        lapply(names(ran_type_tab), function(a){openxlsx::writeDataTable(wb = docs, sheet = paste0('Ran_', a), x = ran_type_tab[[a]], headerStyle = hs_random, rowNames = FALSE, colName = TRUE, sep = '\t')})
        openxlsx::writeDataTable(wb = docs2, sheet = "Exp_inside_microsatellite", x = res$Micro_inside, headerStyle = hs, rowNames = FALSE, colNames = TRUE, sep = '\t')
        openxlsx::writeDataTable(wb = docs2, sheet = "Ran_inside_microsatellite", x = res$Micro_inside_ran, headerStyle = hs_random, rowNames = FALSE, colNames = TRUE, sep = '\t')
        openxlsx::writeDataTable(wb = docs2, sheet = "Exp_dist_microsatellite", x = dist_tab2, headerStyle = hs, rowNames = FALSE, colNames = TRUE, sep = '\t')
        openxlsx::writeDataTable(wb = docs2, sheet = "Ran_dist_microsatellite", x = ran2, headerStyle = hs_random, rowNames = FALSE, colNames = TRUE, sep = '\t')
        cat(paste0('+ Add sheet to file : \"Test_all_type_Repeat\"\n')); openxlsx::addWorksheet(wb = test_docs, sheetName = "Test_all_type_Repeat")
        cat(paste0('+ Add sheet to file : \"Test_all_type_Micro\"\n')); openxlsx::addWorksheet(wb = test_docs, sheetName = "Test_all_type_Micro")
        cat(paste0('+ Add sheet to file : \"Test_by_types_Repeat\"\n')); openxlsx::addWorksheet(wb = test_docs, sheetName = "Test_by_types_Repeat")
        openxlsx::writeDataTable(wb = test_docs, sheet = "Test_all_type_Repeat", x = htab_new1, headerStyle = hs_test, rowNames = FALSE, colNames = TRUE, sep = '\t')
        openxlsx::writeDataTable(wb = test_docs, sheet = "Test_all_type_Micro", x = htab_new2, headerStyle = hs_test, rowNames = FALSE, colNames = TRUE, sep = '\t')
        openxlsx::writeDataTable(wb = test_docs, sheet = "Test_by_types_Repeat", x = gg_tab_test_all1, headerStyle = hs_test, rowNames = FALSE, colNames = TRUE, sep = '\t')
        openxlsx::saveWorkbook(wb = docs, file = paste0(outPath, '/', outFileName, '_Repeat.xlsx'), overwrite = TRUE)
        openxlsx::saveWorkbook(wb = docs2, file = paste0(outPath, '/', outFileName, '_Microsatellite.xlsx'), overwrite = TRUE)
        openxlsx::saveWorkbook(wb = test_docs, file = paste0(outPath, '/', outFileName, '_test_result.xlsx'), overwrite = TRUE)
        res_list = list('Test_result_repeat' = htab_new1, 'Test_result_microsatellite' = htab_new2, 'Test_result_by_repeat_types' = gg_tab_test_all1)
      } else {
        if(dataType == 'cpg'){key = 'CpG'} else if(dataType == 'variant'){key = 'Var'}
        cat(paste0('+ Add sheet to file : \"Exp_inside_', key, '\"\n')); openxlsx::addWorksheet(wb = docs, sheetName = paste0("Exp_inside_", key))
        cat(paste0('+ Add sheet to file : \"Exp_dist_', key, '\"\n')); openxlsx::addWorksheet(wb = docs, sheetName = paste0("Exp_dist_", key))
        cat(paste0('+ Add sheet to file : \"Ran_inside_', key, '\"\n')); openxlsx::addWorksheet(wb = docs, sheetName = paste0("Ran_inside_", key))
        cat(paste0('+ Add sheet to file : \"Ran_dist_', key, '\"\n')); openxlsx::addWorksheet(wb = docs, sheetName = paste0("Ran_dist_", key))
        openxlsx::writeDataTable(wb = docs, sheet = paste0("Exp_inside_", key), x = res[[1]], headerStyle = hs, rowNames = FALSE, colNames = TRUE, sep = '\t')
        openxlsx::writeDataTable(wb = docs, sheet = paste0("Ran_inside_", key), x = res[[6]], headerStyle = hs_random, rowNames = FALSE, colNames = TRUE, sep = '\t')
        openxlsx::writeDataTable(wb = docs, sheet = paste0("Exp_dist_", key), x = dist_tab, headerStyle = hs, rowNames = FALSE, colNames = TRUE, sep = '\t')
        openxlsx::writeDataTable(wb = docs, sheet = paste0("Ran_dist_", key), x = rdist_tab, headerStyle = hs_random, rowNames = FALSE, colNames = TRUE, sep = '\t')
        cat(paste0('+ Add sheet to file : \"Test_all_type_', key, '\"\n')); openxlsx::addWorksheet(wb = test_docs, sheetName = "Test_all_type_Gene")
        openxlsx::writeDataTable(wb = test_docs, sheet = "Test_all_type_Gene", x = htab_new, headerStyle = hs_test, rowNames = FALSE, colNames = TRUE, sep = '\t')
        openxlsx::saveWorkbook(wb = docs, file = paste0(outPath, '/', outFileName, '_', key, '.xlsx'), overwrite = TRUE)
        openxlsx::saveWorkbook(wb = test_docs, file = paste0(outPath, '/', outFileName, '_test_result.xlsx'), overwrite = TRUE)
        res_list = list('Test_result' = htab_new)
      }
      message('- OK!')
    }
  } else {
    stop("[ERROR] Please check res value.\n----- This process is halted. (Time : ", date(), ")\n")
  }
  message('----- Finish. (Time : ', date(), ')')
  return(res_list)
}