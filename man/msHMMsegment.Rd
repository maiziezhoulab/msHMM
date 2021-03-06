% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/msHMM.R
\name{msHMMsegment}
\alias{msHMMsegment}
\title{msHMMsegment}
\usage{
msHMMsegment(
  ematrix,
  annotation,
  tran.namelist = NULL,
  param = NULL,
  autosomes = NULL,
  maxiter = 50,
  tolerance = 0.1 * ncol(ematrix),
  getparam = FALSE,
  verbose = TRUE
)
}
\arguments{
\item{ematrix}{The input expression matrix with rows as genes and columns as
cells.}

\item{annotation}{A data.frame object which should contain the chromosome
name, the start position, and the end position of each gene.}

\item{tran.namelist}{A list containing the cell indexes. The cell indexes in
each index of \code{tran.namelist} can only be ether numeric (integer)
variables or the names of cells. The cells in the same index of
\code{tran.namelist} will share a same transition probability matrix.}

\item{param}{The parameter used in HMM. See details in function
\code{\link[msHMM]{generateParam}}}

\item{autosomes}{Array of LOGICAL values corresponding to the 'chr' argument
where an element is TRUE if the chromosome is an autosome, otherwise FALSE.
If not provided, will automatically set the following chromosomes to false:
"X", "Y", "23", "24", "chrX", chrY", "M", "MT", "chrM". (Copied from HMMcopy)}

\item{maxiter}{The maximum number of iteration for EM algorithm.}

\item{getparam}{If \code{TRUE}, return with the parameters.}

\item{verbose}{If \code{TRUE}, print the messages.}
}
\value{
The segments
}
\description{
Find segmentations across sample.
}
