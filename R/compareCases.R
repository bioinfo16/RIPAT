#' @title Annotate integration sites by genes and TSSs.
#' 
#' @description
#' Annotate vector integration sites by gene data.
#' 
#' @usage 
#' annoByGene(res1, res2,
#'            outPath = getwd(),
#'            outFileName = paste0('RIPAT', round(unclass(Sys.time()))))
#' 
#' @param res1 GR object. This object made by \code{\link{makeExpSet}} function.
#' @param res2 GR object or list. This object is output of \code{\link{makeRanSet}} function.
#' @param outPath String. Plots are saved in this path. \cr Default value is R home directory.
#' @param outFileName Character. Attached ID to the result file name.
#' 
#' @return Return a result list that is made up of insertion and distribution result tables\cr and GenomicRange object of Gene and TSS data.
#'         
#' @examples 
#' data(blast_res); data(blast_res2)
#'
#' blast_gene = compareCases(res1 = blast_res, res2 = blast_res2,
#'                           outFileName = 'blast_res')
#' 
#' @export
compareCases = function(res1, res2, outPath, outFileName){}