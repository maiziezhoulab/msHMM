# This page Original Author: Akdes Serin Harmanci
# https://github.com/akdess/CaSpER
# Modified to perform our multisample RNA-sequencing HMM in CaSpER
# (Same transition probability for tumors or controls)

#' @title runCaSpER_MRHMM()
#'
#' @description  Main casper function that performs a pairwise comparison of all scales from BAF and expression signals to ensure a coherent set of CNV calls.
#'
#' @param object casper object
#'
#' @param removeCentromere boolean values determining if centromere regions should be removed from the analysis
#'
#' @param cytoband cytoband information downloaded from UCSC hg19: http://hgdownload.cse.ucsc.edu/goldenpath/hg19/database/cytoBand.txt.gz hg38:http://hgdownload.cse.ucsc.edu/goldenpath/hg38/database/cytoBand.txt.gz
#'
#' @param method iterative or fixed method. Fixed performs CNV calls on desired baf and expression scale whereas iterative performs pairwise comparison of all expression and baf scale pairs. Iterative method is recommendend. (default: iterative)
#'
#' @return list of objects
#'
#' @export
#'
#'
runCaSpER_MRHMM <- function(object, removeCentromere = T, cytoband = object@cytoband, method = "iterative", maxiter = 100) {
  final.objects <- list()

  if (method == "iterative") {
    loh.list <- list()
    cnv.list <- list()

    message("Performing recursive median filtering...")

    for (i in 1:object@loh.scale) {
      loh.list[[i]] <- lohCallMedianFilterByChr(object, loh.scale = i)
    }

    message("Performing HMM segmentation...")

    for (i in 1:object@cnv.scale) {
      cnv.list[[i]] <- MRHMM::PerformSegmentationWithMRHMM(object, cnv.scale = i, removeCentromere = T, cytoband = cytoband, maxiter)
    }

    combin <- expand.grid(1:object@cnv.scale, 1:object@loh.scale)
    list.names <- apply(combin, 1, function(x) paste(x[1], x[2], sep = "_vs_"))

    for (i in 1:nrow(combin)) {
      loh.scale <- combin[i, 2]
      cnv.scale <- combin[i, 1]
      message("Processing cnv.scale:", cnv.scale, " loh.scale:", loh.scale, "...")
      object <- cnv.list[[cnv.scale]]
      object@loh.median.filtered.data <- loh.list[[loh.scale]]@loh.median.filtered.data
      object <- calculateLOHShiftsForEachSegment(object)
      object <- assignStates(object)
      final.objects[[i]] <- generateLargeScaleEvents(object)
    }
    names(final.objects) <- list.names
  } else if (method == "fixed") {
    object <- MRHMM::PerformSegmentationWithMRHMM(object, cnv.scale = object@cnv.scale, removeCentromere = T, cytoband = cytoband, maxiter)
    object <- lohCallMedianFilterByChr(object, loh.scale = object@loh.scale)
    object <- calculateLOHShiftsForEachSegment(object)
    object <- assignStates(object)
    final.objects[[1]] <- generateLargeScaleEvents(object)
  }
  return(final.objects)
}


PerformSegmentationWithMRHMM <- function(object, cnv.scale, removeCentromere = T, cytoband, maxiter) {

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

  param = MRHMM::generateParam(ematrix)

  segments <- NULL
  hmm.segments.list = MRHMM::MRHMMsegment(ematrix, annotation, tumor.sample.ids,
                                          sample.ids, param = param, autosomes = NULL,
                                          maxiter = maxiter, getparam = FALSE,
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


