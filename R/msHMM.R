## usethis namespace: start
#' @importFrom Rcpp sourceCpp
#' @useDynLib msHMM forward_backward
#' @useDynLib msHMM viterbi
## usethis namespace: end
NULL

#' @title msHMMsegment
#' @description Find segmentations across sample.
#'
#' @param ematrix The input expression matrix with rows as genes and columns as
#' cells.
#' @param annotation A data.frame object which should contain the chromosome
#' name, the start position, and the end position of each gene.
#' @param tran.namelist A list containing the cell indexes. The cell indexes in
#' each index of \code{tran.namelist} can only be ether numeric (integer)
#' variables or the names of cells. The cells in the same index of
#' \code{tran.namelist} will share a same transition probability matrix.
#' @param param The parameter used in HMM. See details in function
#' \code{\link[msHMM]{generateParam}}
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
msHMMsegment <- function(ematrix, annotation, tran.namelist = NULL,
                         param = NULL, autosomes = NULL, maxiter = 50,
                         tolerance = 0.1*ncol(ematrix), getparam = FALSE,
                         verbose = TRUE) {
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

  if(is.null(tran.namelist)){
    tran.namelist = list()
    tran.namelist[[1]] = 1:ncol(ematrix)
  }
  if(length(tran.namelist) > ncol(ematrix))
    stop("The type of transition probabilities is larger than number of cells")
  for(mPT in 1:length(tran.namelist)){
    if(class(tran.namelist[[mPT]]) == "numeric" || class(tran.namelist[[mPT]]) == "integer"){
      if(max(tran.namelist[[mPT]]) > ncol(ematrix)){
        stop("The index in tran.namelist is out of bound")
      }
    }else if(!(all(tran.namelist[[mPT]] %in% colnames(ematrix)))){
      stop("The tran.namelist can only be ether numeric (integer) variables or the names of cells.")
    }
  }

  output.list = EMSegment(ematrix, tran.namelist, chr, autosomes,
                          param, maxiter, tolerance, verbose)
  for(m in 1:dim(ematrix)[2]){
    output.list[[m]]$segs = processSegments(output.list[[m]]$segs, chr,
                                            annotation$start, annotation$end,
                                            ematrix[, m])
  }
  return(output.list)
}
