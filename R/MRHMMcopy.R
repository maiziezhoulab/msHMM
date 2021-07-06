#' @title MRHMMcopy
#' @description Find segmentations across sample.
#'
#' @param ematrix The input expression matrix.
#' @param sample.ids The ids of samples.
#' @param tumor.sample.ids The ids of tumors
#' @param annotation A data.frame object which should contain the chromosome
#' name, the start position, and the end position of each gene.
#' @param param The parameter used in HMM. See details in function
#' \code{\link[MRHMM]{generateParam}}
#' @param autosomes Array of LOGICAL values corresponding to the 'chr' argument
#' where an element is TRUE if the chromosome is an autosome, otherwise FALSE.
#' If not provided, will automatically set the following chromosomes to false:
#' "X", "Y", "23", "24", "chrX", chrY", "M", "MT", "chrM". (Copied from HMMcopy)
#' @param maxiter The maximum number of iteration for EM algorithm.
#' @param getparam If \code{TRUE}, return with the parameters.
#' @param verbose If \code{TRUE}, print the messages.
#'
#' @return The segments
#' @export
#'
#' @examples
MRHMMcopy <- function(ematrix, sample.ids, tumor.sample.ids, annotation, param = NULL,
                      autosomes = NULL, maxiter = 50,
                      getparam = FALSE, verbose = TRUE) {
  chr = annotation$chr

  if (!is.factor(chr)) {
    warning("chr is not a factor, converting to factor")
    chr = as.factor(chr)
  }

  if (is.null(autosomes)) {
    autosomes = (chr != "X" & chr != "Y" & chr != "23" & chr != "24" &
                   chr != "chrX" & chr != "chrY" & chr != "M" & chr != "MT" & chr != "chrM")
  }

  if (is.null(param)) {
    param = generateParam(ematrix)
  }

  if (getparam) {
    return(param)
  }

  output.list = EMSegment(ematrix, sample.ids, tumor.sample.ids, chr, autosomes, param, maxiter,
                              verbose)
  for(m in 1:dim(ematrix)[2]){
    output.list[[m]]$segs = processSegments(output.list[[m]]$segs, chr,
                                            annotation$start, annotation$end,
                                            ematrix[, m])
  }
  return(output.list)
}

# ematrix = readRDS(paste0(dir,"ematrix.rds"))
# annotation = readRDS(paste0(dir,"annotation.rds"))
# load("/Users/zhwen/Documents/R_workspace/CaSpER/CaSpER-master/data/yale_meningioma.rda")
#
# sample.ids = colnames(yale_meningioma[["data"]])
# control.sample.ids = yale_meningioma$control.sample.ids
# tumor.sample.ids = setdiff(sample.ids, control.sample.ids)
# removeCentromere = T
#
# if (removeCentromere) {
#   isCentromer = annotation$isCentromer == "no"
#   ematrix = ematrix[isCentromer, ]
#   annotation = annotation[isCentromer, ]
# }
#
# annotation = data.frame(chr = as.factor(annotation$cytoband),
#                         start = annotation$start,
#                         end = annotation$end)
# #####
# param = generateParam(ematrix)
# MRHMMcopy(ematrix, annotation, param, verbose = F)
#####
# each scale --> change scale outsides

# maxiter = 50
# chrom = annotation$chr
# verbose = T
# autosomes <- (chrom != "X" & chrom != "Y" & chrom != "23" & chrom != "24"
#               & chrom != "chrX" & chrom != "chrY" & chrom != "M" & chrom != "MT"
#               & chrom != "chrM")
