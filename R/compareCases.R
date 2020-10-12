#' @title Annotate integration sites by genes and TSSs.
#' 
#' @description
#' Annotate vector integration sites by gene data.
#' 
#' @usage 
#' compareCases(res_list, excelOut = TRUE,
#'              outPath = getwd(),
#'              outFileName = paste0('RIPAT', round(unclass(Sys.time()))))
#' 
#' @param res_list List. This is the list of result made by \code{\link{annoByGene}}, \code{\link{annoByCpG}},
#'                 \code{\link{annoByRepeat}}, \code{\link{annoByVar}} and \code{\link{annoByCGC}} function.
#' @param excelOut TRUE or FALSE. If user want to make excel file, enter TRUE.\cr
#'                 Default is TRUE.
#' @param outPath String. Outputs are saved in this path. \cr Default value is R home directory.
#' @param outFileName Character. Attached ID to the result file name.
#' 
#' @return Return a result list that is made up of distribution plot and table.
#'         
#' @examples 
#' data(blast_gene)
#' in_list = list(blast_gene, blast_gene)
#'
#' comp_res = compareCases(res_list = in_list,
#'                         excelOut = FALSE,
#'                         outFileName = 'comp_res')
#' 
#' @export
compareCases = function(res_list, excelOut = TRUE, outPath = getwd(), outFileName = paste0('RIPAT', round(unclass(Sys.time())))){
  message('----- Compare two cases. (Time : ', date(), ')')
  message('- Validate options')
  if(!is.logical(excelOut)){stop("[ERROR] 'excelOut' should be logical constants.\n----- This process is halted. (Time : ", date(), ")\n")}
  get_range = unique(lapply(res_list, function(a){as.character(a[[3]]$Range)}))
  if(length(get_range) != 1){stop("[ERROR] All results should be generated with same range and interval number.\n----- This process is halted. (Time : ", date(), ")\n")} else {}
  if(stringr::str_ends(outPath, pattern = '/')){outPath = stringr::word(outPath, start = 1, end = nchar(outPath), sep = '')}
  if(length(res_list) > 5){stop("[ERROR] Maximum length of input list is 5\n----- This process is halted. (Time : ", date(), ")\n")}
  color = c('#8A2BE2', '#5F9EA0', '#FF7F50', '#008080', '#EE82EE', '#1E90FF', '#FF1493', '#FFA500', '#4682B4', '#9ACD32')
  message('- OK!')
  message('- Edit an input list.')
  iname = paste0(c('Gene', 'TSS', 'CpG', 'Repeat', 'Micro', 'Variant'), '_plot_data')
  used_list = lapply(res_list, function(a){a[iname]})
  used_list = lapply(used_list, function(a){a[names(unlist(lapply(a, nrow)))]})
  dname = paste0(stringr::str_remove_all(unlist(lapply(used_list, names)), pattern = '_plot_data'), '_',
                 paste0('C', rep(c(1:length(used_list)), unlist(lapply(used_list, length)))))
  used_list = do.call(c, lapply(used_list, function(a){lapply(a, function(b){b[which(b$Group == 'Observed'),]})}))
  used_list = lapply(c(1:length(dname)), function(a){tmp = used_list[[a]]; tmp$Group = dname[a]; return(tmp)})
  message('- OK!')
  message('- Merge results.')
  plot_tab = data.frame(do.call(rbind, used_list), stringsAsFactors = FALSE)
  plot_tab$Range = factor(plot_tab$Range, levels = unique(unlist(get_range)))
  plot_tab$Group = factor(plot_tab$Group, levels = unique(plot_tab$Group))
  message('- OK!')
  message('- Draw a histogram.')
  grDevices::png(paste0(outPath, '/', outFileName, '_distribution_', organism, '.png'), width = 1200, height = 750)
  res_plot = ggplot2::ggplot(plot_tab) + ggplot2::geom_bar(ggplot2::aes(x = Range, y = Freq, fill = Group), stat = "identity", position = "dodge", width = 0.7) +
    ggplot2::lims(y = c(0, max(plot_tab$Freq)*1.5)) +
    ggplot2::ggtitle(label = "Distribution") +
    ggplot2::xlab('Intervals (Kbs)') + ggplot2::ylab("Ratio of Integration Events") +
    ggplot2::scale_fill_manual(values = color[c(1:length(used_list))]) +
    ggplot2::theme(panel.background = ggplot2::element_rect(fill="white", colour = "white"), panel.grid.major = ggplot2::element_line(size = 0.5, linetype = 'dotted', colour = 'black'),
                   axis.line = ggplot2::element_line(colour = "darkgrey"), legend.title = ggplot2::element_blank(), 
                   legend.key.size = ggplot2::unit(0.7, "cm"), plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 20),
                   legend.text = ggplot2::element_text(size = 18), axis.text = ggplot2::element_text(size = 17),
                   axis.title = ggplot2::element_text(size = 18))
  print(res_plot)
  grDevices::dev.off()
  message('- OK!')
  message('- Do sample test.')
  vc = do.call(cbind, used_list)[,seq(3,length(used_list)*4,by = 4)]
  vn = unlist(lapply(used_list, function(a){tmp = unlist(as.numeric(a$Count)/as.numeric(a$Freq)); return(unlist(tmp[which(!is.nan(tmp))])[1])}))
  test_res = suppressWarnings(lapply(c(1:ncol(vc)), function(a){stats::prop.test(x = as.numeric(vc[a,]), n = vn, alternative = "two.sided", correct = FALSE)}))
  freq_tab = do.call(cbind, used_list)[,seq(4,length(used_list)*4,by = 4)]
  colnames(freq_tab) = dname
  test_tab = data.frame('Range' = unique(unlist(get_range)), freq_tab, 'p_value' = unlist(lapply(test_res, function(a){a$p.value})), stringsAsFactors = FALSE)
  if(excelOut){
    message('- Write a result file.')
    docs = openxlsx::createWorkbook()
    hs = openxlsx::createStyle(fgFill = "#008b8b", halign = "CENTER", textDecoration = "Bold", fontColour = "white")
    cat('+ Add sheet to file : \"TEST_result\"\n'); openxlsx::addWorksheet(wb = docs, sheetName = "TEST_result")
    openxlsx::writeDataTable(wb = docs, sheet = "TEST_result", x = test_tab, headerStyle = hs, rowNames = FALSE, colNames = TRUE, sep = '\t')
    openxlsx::saveWorkbook(wb = docs, file = paste0(outPath, '/', outFileName, '_comp_res.xlsx'), overwrite = TRUE)
    message('- OK!')
  }
  message('----- Finish. (Time : ', date(), ')')
  return(list('plot' = res_plot, 'test' = test_tab))
}
