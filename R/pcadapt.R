#' Principal Component Analysis for outlier detection
#' 
#' \code{pcadapt} performs principal component analysis and computes p-values to test for outliers. The test for
#' outliers is based on the correlations between genetic variation and the first \code{K} principal components. 
#' \code{pcadapt} also allows the user to read outputs from the software PCAdapt, which is implemented in \code{C}. Using the \code{C} software
#' might be useful for very large datasets. \code{pcadapt} also handles Pool-seq data for which the statistical analysis is
#' performed on the genetic markers frequencies. Returns an object of class \code{pcadapt}. 
#'
#' @details First, a principal component analysis is performed on the scaled and centered genotype data. To account for missing
#' data, the correlation matrix between individuals is computed using only the markers available for each
#' pair of individuals. The scores and the loadings (correlations between PCs and genetic markers) are then found using
#' the \code{\link{eigen}} function. Depending on the specified \code{method}, different test statistics can be used. 
#' 
#' \code{mahalanobis} (default): the Mahalanobis distance is computed for each genetic marker using a robust
#' estimate of the mean and of the covariance matrix between the \code{K} vectors of loadings. 
#' 
#' \code{communality}: the communality statistic measures the proportion of variance explained by the first \code{K} PCs.
#' 
#' \code{euclidean}: the Euclidean distance between the \code{K} scaled loadings of each genetic marker and the mean of the \code{K} vectors of scaled loadings is computed.
#' Scaled loadings correspond to loadings divided by a robust estimate of their standard deviation.
#' 
#' \code{componentwise}: returns a matrix of scaled loadings. Scaled loadings correspond to loadings divided by a robust estimate of their standard deviation.
#' 
#' To compute p-values, test statistics (\code{stat}) are divided by a genomic inflation factor (\code{gif}) when \code{method="mahalanobis","euclidean"}. When \code{method="communality"}, the test
#' statistic is first multiplied by \code{K} and divided by the percentage of variance explained by the first \code{K} PCs before accounting for genomic inflation factor. When using \code{method="mahalanobis","communality","euclidean"}, the scaled statistics (\code{chi2_stat}) should follow a chi-squared 
#' distribution with \code{K} degrees of freedom. When using \code{method="componentwise"}, the rescaled loadings should follow a standard normal distribution.
#' For Pool-seq data, \code{pcadapt} provides p-values based on the Mahalanobis distance for each SNP.
#'
#' @param data a data matrix or a data frame if `PCAdapt = FALSE`. The name of the file generated with the \code{C} software PCAdapt (with no extension) if \code{data.type="PCAdapt"}.
#' @param K an integer specifying the number of principal components to retain. In the case where \code{data.type="PCAdapt"}, it is not necessary to set a value for \code{K}, but specifying \code{K}
#' will reduce the number of principal components taken into account.
#' @param method a character string that specifies the test statistic to compute the p-values. Four statistics are currently available,
#' \code{"mahalanobis"}, \code{"communality"}, \code{"euclidean"} and \code{"componentwise"}.
#' @param data.type a character string that specifies the type of data being read, either a \code{genotype} matrix (\code{data.type="genotype"}), or a matrix of allele frequencies (\code{data.type="pool"}, or outputs from the software PCAdapt (\code{data.type="PCAdapt"}).
#' @param min.maf a value between \code{0} and \code{0.5} specifying the threshold of minor allele frequencies above which p-values are computed.
#' @param ploidy an integer specifying the ploidy of the individuals.
#'
#' @return  The returned value \code{x} is an object of class \code{pcadapt}.
#' The different fields can be viewed using the dollar sign (example: \code{x$pvalues}).
#' The returned value contains the following components, depending on the choice of method:
#' \item{stat}{is a vector containing the test statistics associated with the chosen method for each genetic marker. \code{NULL} if \code{method="componentwise"}. \code{method} default value set to \code{mahalanobis}.}
#' \item{pvalues}{is a data frame containing p-values.}
#' \item{maf}{is a vector containing minor allele frequencies.}
#' \item{chi2_stat}{is a vector containing the scaled statistics equal to the values contained in \code{stats} divided by \code{gif} (\code{method}="mahalanobis","euclidean"). It should follow a chi-squared distribution with K degrees of freedom.}
#' \item{gif}{is a numerical value corresponding to the genomic inflation factor estimated from \code{stat}}
#' \item{scores}{is a matrix corresponding to the projections of the individuals onto each PC.}
#' \item{loadings}{is a matrix containing the correlations between each genetic marker and each PC.}
#' \item{singular_values}{contains the ordered squared root of the proportion of variance explained by each PC.}
#' 
#' @examples
#' data <- read4pcadapt("geno3pops",option="example")
#' x <- pcadapt(data,K=10)
#' 
#' ## Screeplot
#' plot(x,option="screeplot")
#' 
#' ## PCA
#' plot(x,option="scores")
#' 
#' ## Neutral SNPs distribution
#' plot(x,option="stat.distribution",K=2)
#' 
#' ## Manhattan Plot
#' plot(x,option="manhattan")
#' 
#' ## Q-Q Plot
#' plot(x,option="qqplot")
#' 
#'@importFrom utils read.table
#'
#' @export
pcadapt = function(data=NULL,K,method="mahalanobis",data.type="genotype",min.maf=0.05,ploidy=2){
  if (data.type == "genotype"){
    res <- corpca(data,K,ploidy=ploidy) 
    nSNP <- dim(data)[2]    
    freq <- apply(data,2,sum)/(ploidy*dim(data)[1])
    res$maf <- pmin(freq,1-freq)
    res$loadings[res$maf<min.maf,] <- NA
  } else if (data.type == "PCAdapt"){
    res <- NULL
    ldgs <- read.table(paste0(data,".loadings"),header=FALSE)
    nSNP <- dim(ldgs)[1]
    res$singular_values <- as.numeric(read.table(paste0(data,".sigma"),header=FALSE))   
    if (missing(K)){
      K <- length(res$singular_values)
    } else {
      if (K > length(res$singular_values)){
        warning("K has to be less than or equal to the number of principal components retained in PCAdapt.")
      }
    }    
    res$loadings <- as.matrix(ldgs[,1:K]*sqrt(nSNP))
    res$scores <- t(read.table(paste0(data,".scores"),header=FALSE))[,1:K]
    out <- read.table(data,header=TRUE)  
    res$maf <- out$mAF 
    res$loadings[res$maf<min.maf,] <- NA
    res$singular_values <- res$singular_values[1:K]  
  } else if (data.type == "pool"){
    nPOP <- dim(data)[1]
    nSNP <- dim(data)[2]
    if (missing(K)){
      K <- nPOP-1
      res <- corpca(data,K,scale=FALSE)
    } else {
      res <- corpca(data,K,scale=FALSE)
    }
    freq <- apply(data,2,sum)/nPOP
    res$maf <- as.vector(pmin(freq,1-freq))
    res$maf[which(is.na(res$maf))] <- 0
  }
  aux <- computeStats(res,method,nSNP,K)
  res$pvalues <- pval(aux,method,K)  
  if (method == "componentwise"){
    res$pvalues[res$maf<min.maf,] <- NA
    resf <- list(pvalues=res$pvalues,maf=res$maf,gif=aux$gif,scores=res$scores,loadings=res$loadings,singular_values=res$singular_values)
  } else {
    res$pvalues[res$maf<min.maf] <- NA
    resf <- list(stat=aux$stat,pvalues=res$pvalues,maf=res$maf,gif=aux$gif,chi2_stat=aux$chi2_stat,scores=res$scores,loadings=res$loadings,singular_values=res$singular_values)
  }
  class(resf) <- 'pcadapt'  
  attr(resf,"method") <- method
  attr(resf,"min.maf") <- min.maf
  attr(resf,"K") <- K
  attr(resf,"data.type") <- data.type
  if (!is.null(attr(data,"pop"))){
    attr(resf,"pop") <- attr(data,"pop")
    resf$names <- colnames(data)
  }
  return(resf)
}

#' Principal Component Analysis based on the correlation matrix
#'
#' \code{corpca} is an auxiliary function that performs principal components analysis on a dataset. It returns an object \code{x} 
#' which contains the loadings, the scores and the singular values of the \code{K} first principal components. 
#' It handles missing values in a dataset and actually computes the eigen elements of the \code{n x n} 
#' covariance matrix, where \code{n} is the number of individuals. 
#' 
#' @param data a data matrix or a data frame. 
#' @param K an integer specifying the number of principal components that are retained.
#' @param scale a logical value indicating whether the data has to be scaled.
#' @param ploidy an integer specifying the ploidy of the individuals.
#' 
#' @importFrom stats cov
#' 
#' @examples
#' data <- read4pcadapt("geno3pops",option="example")
#' x <- corpca(data,K=3)
#' 
#' @keywords internal
#' 
#' @export
corpca = function(data,K,scale=TRUE,ploidy=2){
  n <- dim(data)[1]
  p <- dim(data)[2]
  if (scale==TRUE){
    if (ploidy==2){
      pp <- apply(data,FUN=mean,MARGIN=2)/2
      data_aux <- scale(data,scale=sqrt(2*(pp*(1-pp))))*sqrt(p/(n-1))
    } else {
      pp <- apply(data,FUN=mean,MARGIN=2)
      data_aux <- scale(data,scale=sqrt(pp*(1-pp)))*sqrt(p/(n-1))
    }
  } else {
    data_aux <- scale(data,scale=FALSE)*sqrt(p/(n-1))
  }
  covmat <- cov(t(data_aux),use="pairwise.complete.obs")
  res <- NULL
  aux <- eigen(covmat,symmetric=TRUE)
  sdev <- aux$values[1:K]
  res$scores <- aux$vectors[,1:K]
  aux_ldgs <- t(aux$vectors)%*%data_aux
  res$loadings <- array(0,dim=c(p,K))
  res$loadings[,] <- t((1/(sqrt(sdev)))*aux_ldgs[1:K,])
  res$singular_values <- sqrt(abs(sdev))
  return(res)
}

#' Principal Components p-values
#'
#' \code{pval} computes p-values for each genetic marker, depending on the test statistics distribution.
#' 
#' @param x an output of \code{pcadapt} containing the test statistics.
#' @param method a character string that specifies the method used to compute the p-values.
#' 
#' @examples
#' data <- read4pcadapt("geno3pops",option="example")
#' x <- pcadapt(data,K=3)
#' pval(x,method="mahalanobis",K=3)
#' 
#' @keywords internal
#' 
#' @importFrom stats pnorm pchisq qchisq mad
#' 
#' @export
pval = function(x,method,K){  
  if (method == "componentwise"){ 
    pvalues <- NULL
    column_names <- NULL
    for (k in 1:K){
      sigma <- apply(x$loadings,FUN=function(xx){mad(xx,na.rm=TRUE)},MARGIN=2)
      pvalues <- cbind(pvalues,2*pnorm(abs(x$loadings[,k]/x$gif[k]),mean=0,sd=sigma[k],lower.tail=FALSE))      
      column_names <- c(column_names,paste0("p",k))      
    }
    q <- as.data.frame(pvalues)
    colnames(q) <- column_names
  } else{
    q <- as.numeric(pchisq(x$chi2_stat,df=K,lower.tail=FALSE))
  }  
  return(q)
}

#' Test Statistics
#'
#' \code{computeStats} computes the test statistics for each genetic marker.
#' 
#' @param res a list of quantities among which are the loadings.
#' @param method a character string that specifies the method used to compute the p-values.
#' @param nSNP an integer specifying the number of genetic markers.
#' @param K an integer specifying the number of principal components to retain.
#' 
#' @examples
#' ## see ?pcadapt for examples
#' 
#' @keywords internal
#' 
#' @importFrom stats median mad na.omit
#' @importFrom robust covRob
#' @importFrom MASS cov.rob
#' 
#' @export
computeStats = function(res,method,nSNP,K){
  aux <- NULL
  aux$loadings <- res$loadings
  if (method == "mahalanobis"){
    if (K>1){
      nanlist <- which(!is.na(apply(abs(aux$loadings),1,sum)))
      u <- as.vector(robust::covRob(res$loadings,na.action=na.omit,estim="pairwiseGK")$dist)        
      completed_stat <- array(NA,dim=nSNP)
      completed_stat[nanlist] <- u
      aux$stat <- completed_stat
      aux$gif <- median(aux$stat,na.rm=TRUE)/qchisq(0.5,df=K)
      aux$chi2_stat <- aux$stat/aux$gif
    } else if (K==1){
      nanlist <- which(!is.na(res$loadings[,1]))
      onedcov <- as.vector(MASS::cov.rob(res$loadings[nanlist,1]))
      aux$stat <- (res$loadings[,1]-onedcov$center)^2
      aux$gif <- onedcov$cov[1]
      aux$chi2_stat <- aux$stat/aux$gif
    }
  } else if (method == "euclidean"){
    loadmad <- apply(aux$loadings,FUN=function(xx){mad(xx,na.rm=TRUE)},MARGIN=2)
    aux$stat <- sapply(1:nSNP,FUN=function(xx){sum(res$loadings[xx,1:K]^2/(loadmad^2))})  
    aux$gif <- median(aux$stat,na.rm=TRUE)/qchisq(0.5,df=K)
    aux$chi2_stat <- aux$stat/aux$gif
  } else if (method == "communality"){
    aux$stat <- sapply(1:nSNP,FUN=function(xx){sum(res$loadings[xx,1:K]^2*res$singular_values[1:K]^2/(nSNP))})
    c <- sum(res$singular_values[1:K]^2)/K
    aux$gif <- median(aux$stat*nSNP/c,na.rm=TRUE)/qchisq(0.5,df=K)      
    aux$chi2_stat <- aux$stat*nSNP/(c*aux$gif)
  } else if (method == "componentwise"){
    aux$gif <- sapply(1:K,FUN=function(xx){median(res$loadings[,xx]^2,na.rm=TRUE)/qchisq(0.5,df=1)})
  } 
  return(aux)
}

