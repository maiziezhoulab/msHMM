# Same transition probability for tumor
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
  hmm.segments.list = MRHMMcopy(ematrix, sample.ids, tumor.sample.ids, annotation, param = param,
                        autosomes = NULL, maxiter = 50,
                        getparam = FALSE, verbose = TRUE)
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
generateParam <- function(ematrix) {
  # paste("Parameter missing, ensure all parameters exist as columns in",
  #       "data frame: mu, lambda, nu, kappa, eta, gamma")
  param <- data.frame(strength = 1e+07,
                      e = 0.9999999,
                      mu = quantile(ematrix,
                                    na.rm = TRUE,
                                    prob = c(0.01, 0.05, 0.5, 0.95, 0.99)),
                      lambda = 20,
                      nu = 2.1,
                      kappa = c(0.05, 0.05, 0.8, 0.05, 0.05) * 1000,
                      m = 0,
                      eta = c(5, 5,50, 5, 5) * 10000,
                      gamma = 3,
                      S = 0)
  param$m <- param$mu
  param$S <- ((sd(2^ematrix, na.rm = TRUE)/sqrt(nrow(param)))^2)
  rownames(param) <- seq(1, 5)
  return(param)
}

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

  output.list = manualSegment(ematrix, sample.ids, tumor.sample.ids, chr, autosomes, param, maxiter,
                               verbose)
  for(m in 1:dim(ematrix)[2]){
    output.list[[m]]$segs = processSegments(output.list[[m]]$segs, chr,
                                             annotation$start, annotation$end,
                                             ematrix[, m])
  }
  return(output.list)
}

#
# output = .Call("forward_backward", P.pi, P.T, pyc[m,,char.ind],
#                PACKAGE="HMMcopy")
# maxiter = 50
# chrom = annotation$chr
# verbose = T
# autosomes <- (chrom != "X" & chrom != "Y" & chrom != "23" & chrom != "24"
#               & chrom != "chrX" & chrom != "chrY" & chrom != "M" & chrom != "MT"
#               & chrom != "chrM")

manualSegment <- function(ematrix, sample.ids, tumor.sample.ids, chrom, autosomes, param, maxiter,
                          verbose = TRUE){
  K = dim(param)[1]               # number of states
  # ematrix N*M
  N = nrow(ematrix)               # number of genes
  M = ncol(ematrix)               # number of cells
  # rho = matrix(0, K, N)
  P.gamma = array(0, c(M, K, N))
  # py = matrix(0, K, N)          # Local evidence
  pyc = array(0, c(M, K, N))      # Local evidence
  mus = matrix(0, K, maxiter)     # State means
  # lambdas = matrix(0, K, maxiter) # State Variances (inverse) -- sigma^(-2)
  sigmam2 = matrix(0, K, maxiter) # State Variances (inverse) -- sigma^(-2)
  P.pi = param$kappa               # Initial state distribution
  converged = FALSE               # Flag for convergence
  # Z = rep(0, N)
  # Zc = matrix(0, M, N)
  # Zcounts = matrix(0, K, K)
  P.xi = array(0, c(M, K, K))     # sum the genes of xi
  loglik = rep(0, maxiter)

  # SET UP
  # Set up the chromosome indices and make cell array of chromosome indicies
  char.uniq = levels(chrom)              # WARNING, gets resorted here...

  # initialise the chromosome index and the init state distributions
  char.ind.list = lapply(1:length(char.uniq),
                         function(i) which(chrom == char.uniq[i]))

  # INITIALIZATION
  if (verbose) { message("Initialization") }
  i = 1
  mus[, 1] = param$mu # First column
  sigmam2[, 1] = param$lambda
  # Recalculate the likelihood
  for(m in 1:M){
    for (k in 1:K) {
      # py[k, ] = tdistPDF(copy, mus[k, i], lambdas[k, i], param$nu[k])
      pyc[m, k, ] = tdistPDF(ematrix[, m], mus[k, i], sigmam2[k, i], param$nu[k])
    }
  }
  # Initialize transition matrix to the prior
  # P.T = matrix(0, K, K)
  # Same transition probability for tumor
  num.P.T = 2          # tumor: 1 and control: 2
  P.T = array(0, c(num.P.T, K, K))
  A_prior = array(0, c(num.P.T, K, K))
  dirPrior = array(0, c(num.P.T, K, K))
  for(mPT in 1:num.P.T){
    for (k in 1:K) {
      P.T[mPT, k, ] = (1 - param$e[1]) / (K - 1)
      P.T[mPT, k, k] = param$e[1]
    }
    A_prior[mPT, , ] = P.T[mPT, , ]
    dirPrior[mPT, , ] = P.T[mPT, , ] * param$strength[1]
  }
  loglik[i] = -Inf

  while(!converged && (i < maxiter)) {
    if (verbose) { message(paste("EM iteration:", i,
                                 "Log likelihood:", loglik[i])) }
    i = i + 1

    # E-step
    ############################################################################
    if (verbose) { message("Expectation") }
    for(m in 1:M){
      for (j in 1:length(char.ind.list)) {
        char.ind = char.ind.list[[j]]
        # output = .Call("forward_backward",P.pi,P.T,py[, char.ind],PACKAGE = "HMMcopy")
        # Same transition probability for tumor
        if(sample.ids[m] %in% tumor.sample.ids){
          output = .Call("forward_backward", P.pi, P.T[1, , ], pyc[m, , char.ind],
                         PACKAGE="HMMcopy")
        }else{
          output = .Call("forward_backward", P.pi, P.T[2, , ], pyc[m, , char.ind],
                         PACKAGE="HMMcopy")
        }

        # rho[, char.ind] = output$rho
        P.gamma[m, , char.ind] = output$rho
        loglik[i] = loglik[i] + output$loglik
        # accumulate P.xi for each cell
        P.xi[m, ,] = P.xi[m, ,] + t(colSums(aperm(output$xi, c(3, 2, 1))))
      }
    }

    # M-step
    ############################################################################
    # Update the noise hyperparams

    if (verbose) { message("Maximization") }
    # mu_i = mus[, i - 1]
    # lambda_i = lambdas[, i - 1]
    # output <- estimateTNoiseParamsMap(copy[autosomes], mu_i, lambda_i, param$nu,
    #                                   rho[, autosomes], param$eta, param$m, param$gamma, param$S, param$kappa)
    output = estimate_emi_pi(ematrix[autosomes, ], mus[, i - 1], sigmam2[, i - 1],
                             P.gamma[, , autosomes], param)
    mus[, i] = output$mu.hat
    sigmam2[, i] = output$sigmam2.hat
    P.pi = output$P.pi

    # Recalculate the likelihood
    for(m in 1:M){
      for (k in 1:K) {
        # py[k, ] = tdistPDF(copy, mus[k, i], lambdas[k, i], param$nu[k])
        pyc[m, k, ] = tdistPDF(ematrix[, m], mus[k, i], sigmam2[k, i], param$nu[k])
      }
    }

    # Update transition matrix P.T
    # Same transition probability for tumor
    priorA = rep(0,num.P.T)
    P.T.temp = array(0,dim(P.T))
    for(m in 1:M){
      if(sample.ids[m] %in% tumor.sample.ids){
        P.T.temp[1, , ] = P.xi[m, ,] + P.T.temp[1, , ]
      }else{
        P.T.temp[2, , ] = P.xi[m, ,] + P.T.temp[2, , ]
      }
    }
    for(mPT in 1:num.P.T){
      for (k in 1:K) {
        P.T[mPT, k, ] = P.T.temp[mPT, k, ] + dirPrior[mPT, k, ]
        P.T[mPT, k, ] = normalize(P.T[mPT, k, ])
        priorA[mPT] = priorA[mPT] +
          log(dirichletpdf(A_prior[mPT, k, ], P.T[mPT, k, ]))
      }
    }

    # Compute log likelihood and check convergence
    # Same transition probability for tumor
    for(mPT in 1:num.P.T){
      priorMu = c()
      for(k in 1:K) {
        priorMu[k] = log(dnorm(mus[k, i], param$mu[k], 1))
      }
      loglik[i] = loglik[i] + priorA[mPT] + sum(priorMu)
    }
    if (abs(loglik[i] - loglik[i - 1]) < 1e-1 || loglik[i] < loglik[i - 1]) {
      converged = T
    }
  }

  if(converged){
    i = i - 1
    # Perform one last round of E-step to get latest responsibilities
    #E-step
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if (verbose) { message("Re-calculating latest responsibilties for output") }
    for(m in 1:M){
      for (j in 1:length(char.ind.list)) {
        char.ind = char.ind.list[[j]]
        # output = .Call("forward_backward",P.pi,P.T,py[, char.ind],PACKAGE = "HMMcopy")
        # Same transition probability for tumor
        if(sample.ids[m] %in% tumor.sample.ids){
          output = .Call("forward_backward", P.pi, P.T[1, , ], pyc[m, , char.ind],
                         PACKAGE="HMMcopy")
        }else{
          output = .Call("forward_backward", P.pi, P.T[2, , ], pyc[m, , char.ind],
                         PACKAGE="HMMcopy")
        }

        # rho[, char.ind] = output$rho
        P.gamma[m, , char.ind] = output$rho
        loglik[i] = loglik[i] + output$loglik
        # accumulate P.xi for each cell
        P.xi[m, ,] = P.xi[m, ,] + t(colSums(aperm(output$xi, c(3, 2, 1))))
      }
    }
  }
  if (verbose) {
    message("Optimal parameters found, segmenting and classifying")
  }
  output.list = list()
  for(m in 1:M){
    segs <- list()
    Z = rep(0, N)
    for (j in 1:length(char.ind.list)) {
      char.ind = char.ind.list[[j]]
      if(sample.ids[m] %in% tumor.sample.ids){
        output = .Call("viterbi", log(P.pi), log(P.T[1, , ]), log(pyc[m, , char.ind]),
                       PACKAGE = "HMMcopy")
      }else{
        output = .Call("viterbi", log(P.pi), log(P.T[2, , ]), log(pyc[m, , char.ind]),
                       PACKAGE = "HMMcopy")
      }
      Z[char.ind] <- output$path
      segs[[j]] <- output$seg
    }

    # mus = mus[, 1:i]
    # sigmam2 = sigmam2[, 1:i]
    # loglik = loglik[1:i]
    #
    # output <- vector('list', 0);
    # output$state <- Z
    output$segs <- segs
    # output$mus <- mus
    # output$sigmam2 <- sigmam2
    # output$pi <- P.pi
    # output$loglik <- loglik
    # output$P.gamma <- P.gamma
    output.list[[m]] = output
  }
  return(output.list)
}
# segs.cell = list()
# for(m in 1:M){
#   segs = vector('list', length(chrs))
#   for (j in 1:length(char.ind.list)) {
#     char.ind = char.ind.list[[j]]
#     # output <- .Call("viterbi", log(kappa), log(A), log(py[, I]),
#     #                 PACKAGE = "CaSpER")
#     if(sample.ids[m] %in% tumor.sample.ids){
#       output = .Call("viterbi", log(P.pi), log(P.T[1, , ]), log(pyc[m, , char.ind]),
#                      PACKAGE = "HMMcopy")
#     }else{
#       output = .Call("viterbi", log(P.pi), log(P.T[2, , ]), log(pyc[m, , char.ind]),
#                      PACKAGE = "HMMcopy")
#     }
#
#     # Z[I] = output$path
#     Zc[m, I] = output$path
#     segs[[j]] = output$seg
#   }
#   segs.cell[[m]] = segs
# }
# mus = mus[, 1:i]
# lambdas = lambdas[, 1:i]
# loglik = loglik[1:i]
#
# output = vector('list', 0);
# # output$state = Z
# output$state = Zc
# # output$segs = segs
# output$segs = segs.cell
# output$mus = mus
# output$lambdas = lambdas
# output$pi = kappa
# output$loglik = loglik
# # output$rho = rho
# output$rho = rhoc


# mu = mus[, i - 1]
# sigmam2 = sigmam2[, i - 1]
# P.gamma = P.gamma[, , autosomes]
estimate_emi_pi <- function(ematrix, mu, sigmam2, P.gamma, param) {
  nu = param$nu
  eta = param$eta
  m = param$m
  gamma = param$gamma
  S = param$S
  P.pi.ini = param$kappa
  ematrix.reshape = c(ematrix)
  P.gamma.reshape = matrix(aperm(P.gamma,c(2,3,1)), dim(P.gamma)[2])
  yr = t(matrix(ematrix.reshape, length(ematrix.reshape), length(mu)))        # Vectorize parameter
  u = (1 + nu) / (((yr - mu) ^ 2) * sigmam2 + nu); # scale parameter

  # Calculate the mean
  mu.hat = (rowSums(P.gamma.reshape * u * yr, na.rm = TRUE) + (eta * m)) /
    (rowSums(P.gamma.reshape * u, na.rm = TRUE) + eta)

  # Calculate the precision
  sigmam2.hat = (rowSums(P.gamma.reshape, na.rm = TRUE) + gamma + 1) /
    (rowSums(P.gamma.reshape * u * ((yr - mu.hat) ^ 2), na.rm = TRUE) +
       (eta * (mu.hat - m) ^ 2 + S))

  # Calculate the stationary distribution
  P.pi.hat = (rowSums(P.gamma.reshape, na.rm = TRUE) + P.pi.ini - 1) /
    (length(ematrix.reshape) + sum(P.pi.ini, na.rm = TRUE) + length(mu))

  out = vector('list', 0)
  out$mu.hat = mu.hat
  out$sigmam2.hat = sigmam2.hat
  out$P.pi.hat = P.pi.hat
  return(out)
}


processSegments <- function(seg, chr, start, end, copy) {
  segment <- data.frame()

  if (!is.factor(chr)) {
    warning("chr is not a factor, converting to factor")
    chr <- as.factor(chr)
  }

  chromosomes <- levels(chr)
  for (i in 1:length(chromosomes)) {
    seg_length = dim(seg[[i]])[1]
    chr_name <- rep(chromosomes[i], seg_length)
    chr_index <- which(chr == chromosomes[i])
    chr_start <- start[chr_index][seg[[i]][, 1]]
    chr_stop <- end[chr_index][seg[[i]][, 2]]
    chr_state <- seg[[i]][, 3]
    chr_median <- rep(0, seg_length)
    for(j in 1:seg_length) {
      chr_median[j] <-
        median(na.rm = TRUE, copy[chr_index][seg[[i]][j, 1]:seg[[i]][j, 2]])
    }
    segment <- rbind(segment, cbind(chr = chr_name,
                                    start = as.numeric(chr_start), end = chr_stop, state = chr_state,
                                    median = chr_median))
  }
  segment <- transform(segment, start = as.numeric(as.character(start)),
                       end = as.numeric(as.character(end)),
                       median = as.numeric(as.character(median)))
  return(segment)
}
