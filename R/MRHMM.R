## usethis namespace: start
#' @importFrom Rcpp sourceCpp
#' @useDynLib MRHMM forward_backward
#' @useDynLib MRHMM viterbi
## usethis namespace: end
NULL

#' @title MRHMMsegment
#' @description Find segmentations across sample.
#'
#' @param ematrix The input expression matrix with rows as genes and columns as
#' cells.
#' @param annotation A data.frame object which should contain the chromosome
#' name, the start position, and the end position of each gene.
#' @param tumor.sample.ids The ids of tumors.
#' @param sample.ids The ids of samples.
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
MRHMMsegment <- function(ematrix, annotation, tumor.sample.ids,
                         sample.ids = NULL, param = NULL, autosomes = NULL,
                         maxiter = 50, getparam = FALSE, verbose = TRUE) {
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

  if (is.null(sample.ids)) {
    sample.ids = colnames(ematrix)
    if (is.null(sample.ids))
      stop("No sample ids!")
  }
  if (!all(tumor.sample.ids %in% sample.ids))
    stop("The tumor ids and the sample ids do not match")

  output.list = EMSegment(ematrix, tumor.sample.ids, sample.ids, chr, autosomes,
                          param, maxiter, verbose)
  for(m in 1:dim(ematrix)[2]){
    output.list[[m]]$segs = processSegments(output.list[[m]]$segs, chr,
                                            annotation$start, annotation$end,
                                            ematrix[, m])
  }
  return(output.list)
}
