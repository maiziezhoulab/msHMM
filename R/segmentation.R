#' @title generateParam
#' @description Generate parameters for HMM
#' @param ematrix The input expression matrix.
#'
#' @return Parameters formatted in a data.frame.
#' @details
#' \code{e}: Probability of extending a segment, increase to lengthen segments,
#' decrase to shorten segments. Range: (0, 1)
#'
#' \code{strength}: Strength of initial e suggestion, reducing allows e to
#' change, increasing makes e undefiable. Range: [0, Inf)
#'
#' \code{mu}: Suggested median for copy numbers in state, change to readjust
#' classification of states. Range: (-Inf, Inf)
#'
#' \code{lambda}: Suggested precision (inversed variance) for copy numbers in
#' state, increase to reduce overlap between states. Range: [0, Inf)
#'
#' \code{nu}: Suggested degree of freedom between states, increase to reduce
#' overlap between states. Range: [0, Inf)
#'
#' \code{kappa}: Suggested distribution of states. Should sum to 1.
#'
#' \code{m}: Optimal value for mu, difference from corresponding mu value
#' determines elasticity of the mu value. i.e. Set to identical value as mu if you don't want mu to move much.
#'
#' \code{eta}: Mobility of mu, increase to allow more movement. Range: [0, Inf)
#'
#' \code{gamma}: Prior shape on lambda, gamma distribution. Effects flexibility
#' of lambda.
#'
#' \code{S}: Prior scale on lambda, gamma distribution. Effects flexibility of
#' lambda.

#' @export
#'
#' @examples
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

EMSegment <- function(ematrix, tumor.sample.ids, sample.ids, chrom, autosomes,
                      param, maxiter, verbose = TRUE){
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

  while(!converged && (i <= maxiter)) {
    if (verbose) { message(paste("EM iteration:", i,
                                 "Log likelihood:", loglik[i])) }
    i = i + 1

    # E-step
    ############################################################################
    if (verbose) { message("Expectation") }
    for(m in 1:M){
      for (j in 1:length(char.ind.list)) {
        char.ind = char.ind.list[[j]]
        # output = .Call("forward_backward",P.pi,P.T,py[, char.ind],PACKAGE = "MRHMM")
        # Same transition probability for tumor
        if(sample.ids[m] %in% tumor.sample.ids){
          output = .Call("forward_backward", P.pi, P.T[1, , ], pyc[m, , char.ind],
                         PACKAGE="MRHMM")
        }else{
          output = .Call("forward_backward", P.pi, P.T[2, , ], pyc[m, , char.ind],
                         PACKAGE="MRHMM")
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
    if (verbose) {
      message("Optimal parameters found, segmenting and classifying")
    }

    i = i - 1
    # Perform one last round of E-step to get latest responsibilities
    #E-step
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if (verbose) { message("Re-calculating latest responsibilties for output") }
    for(m in 1:M){
      for (j in 1:length(char.ind.list)) {
        char.ind = char.ind.list[[j]]
        # output = .Call("forward_backward",P.pi,P.T,py[, char.ind],PACKAGE = "MRHMM")
        # Same transition probability for tumor
        if(sample.ids[m] %in% tumor.sample.ids){
          output = .Call("forward_backward", P.pi, P.T[1, , ], pyc[m, , char.ind],
                         PACKAGE="MRHMM")
        }else{
          output = .Call("forward_backward", P.pi, P.T[2, , ], pyc[m, , char.ind],
                         PACKAGE="MRHMM")
        }

        # rho[, char.ind] = output$rho
        P.gamma[m, , char.ind] = output$rho
        loglik[i] = loglik[i] + output$loglik
        # accumulate P.xi for each cell
        P.xi[m, ,] = P.xi[m, ,] + t(colSums(aperm(output$xi, c(3, 2, 1))))
      }
    }
  }else{
    if (verbose)
      message("Reach maxiter, segmenting and classifying")
  }

  mus = mus[, 1:i]
  sigmam2 = sigmam2[, 1:i]
  loglik = loglik[1:i]

  output.list = list()
  for(m in 1:M){
    segs <- list()
    Z = rep(0, N)
    for (j in 1:length(char.ind.list)) {
      char.ind = char.ind.list[[j]]
      if(sample.ids[m] %in% tumor.sample.ids){
        output = .Call("viterbi", log(P.pi), log(P.T[1, , ]), log(pyc[m, , char.ind]),
                       PACKAGE = "MRHMM")
      }else{
        output = .Call("viterbi", log(P.pi), log(P.T[2, , ]), log(pyc[m, , char.ind]),
                       PACKAGE = "MRHMM")
      }
      Z[char.ind] <- output$path
      segs[[j]] <- output$seg
    }

    output <- vector('list', 0);
    output$state <- Z
    output$segs <- segs
    output$mus <- mus
    output$sigmam2 <- sigmam2
    output$pi <- P.pi
    output$loglik <- loglik
    output$P.gamma <- P.gamma
    output.list[[m]] = output
  }
  return(output.list)
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

  # Calculate the precision (the inverse of the variance)
  sigmam2.hat = (rowSums(P.gamma.reshape, na.rm = TRUE) + gamma + 1) /
    (rowSums(P.gamma.reshape * u * ((yr - mu.hat) ^ 2), na.rm = TRUE) +
       (eta * (mu.hat - m) ^ 2 + S))

  # Calculate the start distribution
  P.pi.hat = (rowSums(P.gamma.reshape, na.rm = TRUE) + P.pi.ini - 1) /
    (length(ematrix.reshape) + sum(P.pi.ini, na.rm = TRUE) + length(mu))

  out = vector('list', 0)
  out$mu.hat = mu.hat
  out$sigmam2.hat = sigmam2.hat
  out$P.pi.hat = P.pi.hat
  return(out)
}

