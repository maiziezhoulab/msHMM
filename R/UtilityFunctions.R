# Author: Daniel Lai
# https://bioconductor.org/packages/release/bioc/html/HMMcopy.html
# Normalize a given array to sum to 1
normalize <- function(A) {
  M <- A / (sum(A) + (sum(A) == 0))
  return(M);
}


# Dirichlet probability density function, returns the probability of vector
# x under the Dirichlet distribution with parameter vector alpha
# Author: David Ross http://www.cs.toronto.edu/~dross/code/dirichletpdf.m
dirichletpdf <- function(x, alpha) {
  if (any(x < 0)) {
    return(0);
  }
  if (abs(sum(x) - 1) > 1e-3) {
    stop("Dirichlet PDF: probabilities do not sum to 1")
    return(0);
  }
  p <- exp(lgamma(sum(alpha)) - sum(lgamma(alpha))) * prod(x ^ (alpha - 1))
  return(p)
}


# Author: Daniel Lai
# https://bioconductor.org/packages/release/bioc/html/HMMcopy.html
tdistPDF <- function(x, mu, lambda, nu) {
  p <- (gamma(nu / 2 + 0.5)/gamma(nu / 2)) * ((lambda / (pi * nu)) ^ (0.5)) *
    (1 + (lambda * (x - mu) ^ 2) / nu) ^ (-0.5 * nu - 0.5)
  p[is.na(p)] <- 1
  return(p)
}


# Original Author: Akdes Serin Harmanci
# https://github.com/akdess/CaSpER
# Modified to perform our multisample RNA-sequencing HMM in CaSpER method
# (Same transition probability for tumors or controls)
PerformSegmentationWithHMM <- function(object, cnv.scale, removeCentromere = T, cytoband) {

  ematrix <- object@control.normalized[[cnv.scale]]
  annotation <- object@annotation.filt[]

  sample.ids = colnames(ematrix)
  control.sample.ids = object@control.sample.ids
  tumor.sample.ids = setdiff(sample.ids, control.sample.ids)

  if (removeCentromere) {
    isCentromer = annotation$isCentromer == "no"
    ematrix = ematrix[isCentromer, ]
    annotation = annotation[isCentromer, ]
  }

  annotation = data.frame(chr = as.factor(annotation$cytoband),
                          start = annotation$start,
                          end = annotation$end)

  param = generateParam(ematrix)

  segments <- NULL
  hmm.segments.list = MRHMMsegment(ematrix, annotation, tumor.sample.ids,
                                   sample.ids, param = param, autosomes = NULL,
                                   maxiter = 50, getparam = FALSE,
                                   verbose = TRUE)
  for (i in 1:dim(ematrix)[2]) {
    hmm.segments <- hmm.segments.list[[i]]
    segments <- rbind(segments, data.frame(ID = colnames(data)[i], hmm.segments$segs))
  }

  arms <- paste(cytoband$V1, cytoband$V4, sep = "")
  arm_sizes <- cytoband$V3 - cytoband$V2

  object@segments <- segments
  object@segments$event_scale <- rep("", nrow(object@segments))

  for (i in 1:dim(object@segments)[1]) {
    ind <- which(annotation$Position >= object@segments$start[i] & annotation$cytoband == object@segments[i, "chr"] & annotation$Position <=
                   object@segments$end[i])
    object@segments$size[i] <- object@segments$end[i] - object@segments$start[i]
    object@segments$num.marks[i] <- length(ind)
    pair_arm_sizes <- arm_sizes[match(as.character(object@segments[i, "chr"]), arms)]
    object@segments$arm.size.perc[i] <- object@segments$size[i]/pair_arm_sizes
    if (object@segments$size[i] > pair_arm_sizes * (1/3))
      object@segments$event_scale[i] <- "large_scale"
    if ((object@segments$size[i] > pair_arm_sizes * (1/10)) & (object@segments$size[i] < pair_arm_sizes * (1/3)))
      object@segments$event_scale[i] <- "focal"

  }

  return(object)
}

