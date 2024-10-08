################################################################################

#' @useDynLib pcadapt, .registration = TRUE
#' @importFrom Rcpp sourceCpp
NULL

#' Principal Component Analysis for outlier detection
#'
#' \code{pcadapt} performs principal component analysis and computes p-values to
#' test for outliers. The test for outliers is based on the correlations between
#' genetic variation and the first \code{K} principal components. \code{pcadapt}
#' also handles Pool-seq data for which the statistical analysis is performed on
#' the genetic markers frequencies. Returns an object of class \code{pcadapt}.
#'
#' @details First, a principal component analysis is performed on the scaled and 
#' centered genotype data. Depending on the specified \code{method}, different 
#' test  statistics can be used.
#'
#' \code{mahalanobis} (default): the robust Mahalanobis distance is computed for 
#' each genetic marker using a robust estimate of both mean and covariance 
#' matrix between the \code{K} vectors of z-scores.
#'
#' \code{communality}: the communality statistic measures the proportion of 
#' variance explained by the first \code{K} PCs. Deprecated in version 4.0.0.
#'
#' \code{componentwise}: returns a matrix of z-scores.
#'
#' To compute p-values, test statistics (\code{stat}) are divided by a genomic 
#' inflation factor (\code{gif}) when \code{method="mahalanobis"}. When using 
#' \code{method="mahalanobis"}, the scaled statistics 
#' (\code{chi2_stat}) should follow a chi-squared distribution with \code{K} 
#' degrees of freedom. When using \code{method="componentwise"}, the z-scores 
#' should follow a chi-squared distribution with \code{1} degree of freedom. For 
#' Pool-seq data, \code{pcadapt} provides p-values based on the Mahalanobis 
#' distance for each SNP.
#'
#' @param input The output of function \code{read.pcadapt}.
#' @param K an integer specifying the number of principal components to retain.
#' @param method a character string specifying the method to be used to compute
#'   the p-values. Two statistics are currently available, \code{"mahalanobis"},
#'   and \code{"componentwise"}.
#' @param min.maf Threshold of minor allele frequencies above which p-values are 
#'   computed. Default is \code{0.05}.
#' @param LD.clumping Default is \code{NULL} and doesn't use any SNP thinning.
#'   If you want to use SNP thinning, provide a named list with parameters 
#'   \code{$size} and \code{$thr} which corresponds respectively to the window 
#'   radius and the squared correlation threshold. A good default value would 
#'   be \code{list(size = 500, thr = 0.1)}.
#' @param pca.only a logical value indicating whether PCA results should be 
#'   returned (before computing any statistic).
#' @param ploidy Number of trials, parameter of the binomial distribution. 
#'   Default is 2, which corresponds to diploidy, such as for the human genome.
#' @param tol Convergence criterion of \code{RSpectra::svds()}. 
#'   Default is \code{1e-4}.
#' 
#' @return The returned value is an object of class \code{pcadapt}.
#' 
#' @useDynLib pcadapt
#' 
#' @name pcadapt
#'
#' @export
#'
pcadapt <- function(input, 
                    K = 2, 
                    method = "mahalanobis", 
                    min.maf = 0.05, 
                    ploidy = 2,
                    LD.clumping = NULL,
                    pca.only = FALSE,
                    tol = 1e-4) {
  
  if (missing(input)) {
    appDir <- system.file("shiny-examples/app-pcadapt", package = "pcadapt")
    if (appDir == "") {
      stop("Could not find Shiny app in pcadapt.", call. = FALSE)
    }
    shiny::runApp(appDir, display.mode = "normal")  
  } else {
    if (any(grepl("^pcadapt_", class(input)))) {
      UseMethod("pcadapt")
    } else {
      stop("Remember to always use read.pcadapt() before using pcadapt().")
    }
  }
}

#' @rdname pcadapt
#' @export
pcadapt.pcadapt_matrix <- function(input, 
                                   K = 2, 
                                   method = c("mahalanobis", "componentwise"), 
                                   min.maf = 0.05, 
                                   ploidy = 2,
                                   LD.clumping = NULL,
                                   pca.only = FALSE,
                                   tol = 1e-4) {
  
  pcadapt0(input, K, match.arg(method), min.maf, ploidy, LD.clumping, pca.only, tol)
}

#' @rdname pcadapt
#' @export
pcadapt.pcadapt_bed <- function(input, 
                                K = 2, 
                                method = c("mahalanobis", "componentwise"), 
                                min.maf = 0.05, 
                                ploidy = 2,
                                LD.clumping = NULL,
                                pca.only = FALSE, 
                                tol = 1e-4) {
  
  # File mapping
  n <- attr(input, "n")
  p <- attr(input, "p")
  xptr <- structure(bedXPtr(input, n, p), n = n, p = p, class = "xptr_bed")
  
  pcadapt0(xptr, K, match.arg(method), min.maf, ploidy, LD.clumping, pca.only, tol)
}

#' @rdname pcadapt
#' @export
pcadapt.pcadapt_pool <- function(input, 
                                 K = (nrow(input) - 1),
                                 method = "mahalanobis",
                                 min.maf = 0.05,
                                 ploidy = NULL,
                                 LD.clumping = NULL,
                                 pca.only = FALSE,
                                 tol) {
  
  w <- matrix(NA_real_, nrow = ncol(input), ncol = K)
  
  tmat <- scale(input, center = TRUE, scale = FALSE) 
  tmat[is.na(tmat)] <- 0 # mean imputation
  
  mean_freq <- attr(tmat, "scaled:center")
  mean_freq <- pmin(mean_freq, 1 - mean_freq)

  pass <- mean_freq > min.maf
  
  if (nrow(input) == 2) {
    obj.pca <- list(
      u = matrix(NA_real_, nrow = 2, ncol = 1),
      v = t(tmat[1, pass, drop = FALSE]),
      d = 1
    )
  } else {
    obj.pca <- svd(tmat[, pass, drop = FALSE])
    #obj.pca <- RSpectra::svds(tmat[, pass, drop = FALSE], k = K)
  }
  
  w[pass, ] <- obj.pca$v[, 1:K, drop = FALSE]
  res <- get_statistics(w, method = method, pass = pass)

  structure(
    list(
      scores = obj.pca$u[, 1:K, drop = FALSE],
      singular.values = sqrt(obj.pca$d[1:K]^2 / sum(obj.pca$d^2)) ,
      loadings = w,
      zscores = w,
      af = attr(tmat, "scaled:center"),
      maf = mean_freq,
      chi2.stat = res$chi2.stat,
      stat = res$stat,
      gif = res$gif,
      pvalues = res$pvalues,
      pass = which(pass)
    ),
    K = K, method = method, min.maf = min.maf, class = "pcadapt"
  )
}

################################################################################

#' pcadapt statistics
#'
#' \code{get_statistics} returns chi-squared distributed statistics. 
#'
#' @param zscores a numeric matrix containing the z-scores.
#' @param method a character string specifying the method to be used to compute
#' the p-values. Two statistics are currently available, \code{"mahalanobis"},
#' and \code{"componentwise"}.
#' @param pass a boolean vector.
#' 
#' @return The returned value is a list containing the test statistics and the 
#' associated p-values.
#' 
#' @importFrom stats median na.omit pchisq qchisq
#' 
#' @keywords internal
#'
get_statistics <- function(zscores, method, pass) {
  
  nSNP <- nrow(zscores)
  K <- ncol(zscores)
  if (method == "mahalanobis") {
    res <- rep(NA, nSNP)
    if (K == 1) {
      res[pass] <- (zscores[pass] - median(zscores[pass]))^2 
    } else if (K > 1) {
      res[pass] <- bigutilsr::dist_ogk(zscores[pass, ])
    }
    gif <- median(res, na.rm = TRUE) / qchisq(0.5, df = K)
    res.gif <- res / gif
    pval <- as.numeric(pchisq(res.gif, df = K, lower.tail = FALSE))
  } else if (method == "communality") {
    # res <- sapply(1:nSNP, FUN = function(h) {sum(zscores[h, ]^2 * values^2 / nSNP)})
    # c <- sum(values^2) / K
    # gif <- median(res * nSNP / c, na.rm = TRUE) / qchisq(0.5, df = K)
    # res.gif <- res * nSNP / (c * gif)
    # pval <- pchisq(res.gif, df = K, lower.tail = FALSE)
  } else if (method == "componentwise") {
    res <- apply(zscores, MARGIN = 2, FUN = function(h) {h^2})
    gif <- sapply(1:K, FUN = function(h) {
      median(zscores[, h]^2, na.rm = TRUE) / qchisq(0.5, df = 1)
    })
    res.gif = res / gif
    pval <- NULL
    for (k in 1:K) {
      pval <- cbind(pval, pchisq(res.gif[, k], df = 1, lower.tail = FALSE))
    }
  }
  return(list(stat = res, gif = gif, chi2.stat = res.gif, pvalues = pval))
}

################################################################################

pcadapt0 <- function(input, K, method, min.maf, ploidy, LD.clumping, pca.only, tol) {
  
  # Test arguments and init
  if (!(class(K) %in% c("numeric", "integer")) || K <= 0)
    stop("K has to be a positive integer.")
  
  if (!is.numeric(min.maf) || min.maf < 0 || min.maf > 0.45) 
    stop("min.maf has to be a real number between 0 and 0.45.")
  
  # Compute PCs and z-scores    
  obj.pca <- iram_and_reg(input, K = K, 
                          min.maf = min.maf, 
                          ploidy = ploidy,
                          LD.clumping = LD.clumping, 
                          tol = tol)
  if (pca.only) return(obj.pca)
  
  res <- get_statistics(obj.pca$zscores, 
                        method = method, 
                        pass = obj.pca$pass)
  
  structure(
    list(
      scores = obj.pca$u,
      singular.values = obj.pca$d / sqrt(obj.pca$total_var),
      loadings = obj.pca$v,
      zscores = obj.pca$zscores,
      af = obj.pca$af,
      maf = pmin(obj.pca$af, 1 - obj.pca$af),
      chi2.stat = res$chi2.stat,
      stat = res$stat,
      gif = res$gif,
      pvalues = res$pvalues,
      pass = obj.pca$pass
    ),
    K = K, method = method, min.maf = min.maf, class = "pcadapt"
  )
}

################################################################################

#' Shiny app
#'
#' \code{pcadapt} comes with a Shiny interface.
#'
#' @export
run.pcadapt <- function() {
  appDir <- system.file("shiny-examples", "app-pcadapt", package = "pcadapt")
  if (appDir == "") {
    stop("Could not find example directory. Try re-installing `pcadapt`.", 
         call. = FALSE)
  }
  shiny::runApp(appDir, launch.browser = TRUE)
}

################################################################################
