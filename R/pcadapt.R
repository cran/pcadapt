#' Auxiliary function
#' 
#' \code{fdrList} returns a list of false discovery rates for different thresholds.
#' 
#' @param snps a vector of indices.
#' @param adjpval a vector of adjusted p-values.
#' @param list_selection a vector of indices containing the list of all the SNPs under selection.
#' @param threholds a vector containing the control thresholds. Default value set to \code{seq(0.01,0.50,length=100)}.
#'   
#' @examples
#' x <- read4pcadapt("geno3pops",option="example")
#' y <- pcadapt(x,K=3)
#' P <- pval(y)
#' ranked_snps <- order(P[,1],decreasing=FALSE)
#' adjpval <- p.adjust(sort(P[,1]),method="BH",n=length(P[,1]))
#' fdrList(ranked_snps,adjpval,1:250,thresholds=seq(0.05,0.5,length=10))
#' 
#' @keywords internal
#' 
#' @export
fdrList = function(snps,adjpval,list_selection,thresholds=seq(0.01,0.50,length=100)){
  fdr <- NULL  
  for (a in thresholds)
  {
    aux <- snps[adjpval<=a]
    fdr <- c(fdr,fdrCalc(aux,list_selection))
  }  
  return(fdr)
}

#' Auxiliary function
#' 
#' \code{sdevList} returns a list of estimated standard deviations using the kurtosis-based method.
#' 
#' @param V a numeric matrix containing the loadings of the Principal Components Analysis.
#' @param K an integer specifying the number of retained principal components.
#'   
#' @examples
#' x <- replicate(1000, rnorm(100))
#' x <- floor(abs(x))
#' freq <- apply(x,2,sum)/(2*dim(x)[1])
#' minfreq <- pmin(freq,1-freq)
#' respc <- prcomp(x)
#' V <- respc$rotation*sqrt(dim(x)[2])
#' sdevList(V,K=3,maf=minfreq,minmaf=0.05)
#' 
#' @keywords internal
#' 
#' @export
sdevList = function(V,K,maf,minmaf){
  out <- NULL
  if (K==1){
    idxmaf <- maf>=minmaf
    aux <- sdevEstimation(V[idxmaf])
    out$neutral_sdev <- c(aux$neutral_sdev)
    out$prop_removed <- c(aux$prop_removed)
    out$q <- cbind(out$q,aux$resq)
  } else {
    idxmaf <- maf>=minmaf
    for (k in 1:K){
      aux <- sdevEstimation(V[idxmaf,k])
      out$neutral_sdev <- c(out$neutral_sdev,aux$neutral_sdev)
      out$q <- cbind(out$q,aux$resq)
      out$prop_removed <- c(out$prop_removed,aux$prop_removed)
    }
  }
  out$p <- aux$vecp
  out$q <- data.frame(out$q)
  return(out)
}

#' False Discovery Rate
#'
#' \code{fdrCalc} computes the false discovery rate (FDR) for a given list of discoveries,
#' knowing the list of the Single Nucleotide Polymorphisms (SNPs) under selection.
#' 
#' @param x a vector of indices corresponding to the discoveries.
#' @param list_selection a vector of indices containing the list of all the SNPs under selection.
#' 
#' @examples
#' x <- c(1,4,5,8,9)
#' list_selection <- c(4,5,7,10,11)
#' fdrCalc(x,list_selection)
#'
#' @export
fdrCalc = function(x,list_selection){
  if (length(x)==0){
    fdr <- 0
  } else {
    fdr <- sum(1-(x %in% list_selection))/length(x)
  }
  return(fdr)  
}

#' Kurtosis
#'
#' \code{kurtosis} computes the kurtosis of a given distribution \code{x}.
#' 
#' @param x a numeric vector. 
#' @param prop a value between \code{0} and \code{1} specifying the proportion of the distribution that is removed when computing the kurtosis. 
#' Default value set to \code{0}.
#' 
#' @examples
#' x <- rnorm(10000)
#' kurtosis(x, prop = 0.2)
#' 
#' @importFrom stats quantile
#'
#' @export
kurtosis = function(x,prop=0){
  x_rm <- x[which(!is.na(x))]
  q <- quantile(abs(x_rm),1-prop)
  aux <- x_rm[abs(x_rm)<q]
  mu <- mean(aux)
  kurt <- mean((aux-mu)^4)/mean((aux-mu)^2)^2 - 3 # + abs(mean((aux-mu)^3)/(mean((aux-mu)^2))^(3/2))
  return(kurt)
}

#' False Discovery Rate Control Plotting
#'
#' \code{fdrControl} plots the false discovery rate as a function of thresholds under which the FDR are expected to be.
#' 
#' @param x an output from \code{outlier} containing the false discovery rates for each threshold that is specified in \code{alpha_thresholds}.
#' @param alpha_threholds a vector containing the control thresholds. Default value set to \code{seq(0.01,0.50,length=100)}.
#' 
#' @examples
#' x <- NULL
#' x$fdr <- seq(0.0,0.4,length=100)
#' fdrControl(x)
#' 
#' @keywords internal
#' 
#' @importFrom graphics abline axis plot.default title
#'
#' @export
fdrControl = function(x,alpha_thresholds=seq(0.01,0.50,length=100)){
  plot.default(alpha_thresholds,x$fdr,col="red",cex=0.2,xlim=c(0,0.5),ylim=c(0,0.5),xlab="thresholds",ylab="FDR",lwd=10,type="l",xaxt="n",yaxt="n")
  tick_marks <- seq(0,0.5,by=0.05)
  axis(1,at=tick_marks)
  axis(2,at=tick_marks)
  abline(v=tick_marks,col="black",lty="dotted")
  abline(h=tick_marks,col="black",lty="dotted")
  abline(0,1)
  title("FDR Control")
}

#' Principal Components Analysis
#' 
#' \code{pcadapt} performs principal component analysis and compute P-values to test for selection as indicated
#' by significant correlations between genetic variation and principal components. \code{pcadapt} also allows the user to read PCAdapt outputs 
#' for further analysis in \code{R}. Returns an object of class \code{pcadapt}.
#'
#' @details First, a principal component analysis is performed on the scaled and centered genotype data. To account for missing
#' data, the correlation matrix between individuals is computed using only the markers available for each
#' pair of individuals. The scores and the loadings (correlations between PCs and genetic markers) are then found using
#' the \code{\link{eigen}} function. The p-values are then computed based on the matrix of loadings. The loadings of the neutral markers
#' are assumed to follow a centered Gaussian distribution. The standard deviation of the Gaussian distribution
#' is estimated after removing a proportion of genetic markers with the largest loadings (in absolute values).
#' The removal proportion is the smallest percentage such that the kurtosis of the truncated distribution of the loadings matches the kurtosis 
#' of a Gaussian distribution, which is equal to 3. The standard deviation of the loadings is finaly estimated using the maximum likelihood of a truncated Gaussian distribution.
#'
#' @param data a data matrix or a data frame. 
#' @param file the name of the file which the data are to be read from. Basically, the file generated with PCAdapt which has no extension. 
#' @param K an integer specifying the number of principal components to retain.
#' @param communality_test a logical value indicating whether a communality test should be performed. Default value set to \code{FALSE}.
#' @param ploidy an integer specifying the ploidy of the individuals.
#' @param minmaf a value between \code{0} and \code{0.5} specifying the threshold under which the frequencies are considered minor allele frequencies.
#'
#' @return  The returned value \code{x} is an object of class \code{pcadapt}.
#' The different fields can be viewed using the dollar sign (example: \code{x$neutral_sdev}).
#' The returned value contains the following components:
#' \item{loadings}{is a matrix containing the correlations between each genetic marker and each PC.}
#' \item{scores}{is a matrix corresponding to the projections of the individuals onto each PC.}
#' \item{singular_values}{contains the ordered squared root of the proportion of variance explained by each PC.}
#' \item{pvalues}{is a data frame containing the p-values for the \code{K} first principal components.}
#' \item{communality}{contains the communality for each PC which corresponds to the proportion of variance explained by the first \code{K} PCs.}
#' \item{p}{is a data frame with \code{K} columns. Gives the proportions removed from the loadings distributions in order to estimate the standard deviation of the neutral markers.}
#' \item{q}{is a data frame with \code{K} columns. Each column of \code{q} represents the kurtosis evaluated on the distribution of the loadings for each cut-off provided by \code{p}.}
#' \item{proportion_removed}{is a list of size \code{K} corresponding to the proportions of markers to remove from the loading distributions to match the kurtosis expected for a Gaussian distribution.}
#'
#' @examples
#' x <- read4pcadapt("geno3pops",option="example")
#' x <- floor(abs(x))
#' y <- pcadapt(x,K=10)
#' 
#' ## Screeplot
#' plot(y,option="screeplot")
#' 
#' ## PCA
#' plot(y,option="scores")
#' 
#' ## Neutral SNPs distribution
#' plot(y,option="neutral",K=1)
#' 
#' ## Manhattan Plot
#' plot(y,option="manhattan",K=1)
#' 
#' ## Q-Q Plot
#' plot(y,option="qqplot",K=1)
#' 
#'@importFrom utils read.table
#'
#' @export
pcadapt = function(data=NULL,file=NULL,K,communality_test=FALSE,ploidy=2,minmaf=0.05){
  res <- NULL
  if (is.null(file)){
    respc <- corpca(data,K)
    V <- respc$loadings
    freq <- apply(data,2,sum)/(ploidy*dim(data)[1])
    minfreq <- pmin(freq,1-freq)
    aux <- sdevList(V,K,maf=minfreq,minmaf=minmaf)
    res$neutral_sdev <- aux$neutral_sdev
    res$loadings <- V
    nSNP <- dim(data)[2]
    res$correlations <- V/sqrt(nSNP)
    res$scores <- respc$scores
    res$singular_values <- sqrt(abs(respc$sdev))
    res$maf <- minfreq
    if (K > 1){
      res$communality <- sapply(1:nSNP,FUN=function(xx){sum(res$loadings[xx,1:K]^2*res$singular_values[1:K]^2/nSNP)})
    } else {
      res$communality <- res$singular_values[1]^2*V^2/dim(data)[2]
    }  
    res$p <- aux$p
    res$q <- aux$q
    res$proportion_removed <- aux$prop_removed
    auxpval <- pval(res,communality_test=communality_test)
    auxpval[minfreq<minmaf,] <- NA
    res$pvalues <- auxpval
  } else {
    ldgs <- read.table(paste0(file,".loadings"),header=FALSE)
    res$loadings <- ldgs*sqrt(dim(ldgs)[1])
    res$correlations <- (res$loadings)/sqrt(dim(ldgs)[1])
    res$scores <- t(read.table(paste0(file,".scores"),header=FALSE))
    res$singular_values <- read.table(paste0(file,".sigma"),header=FALSE) 
    out <- read.table(file,header=TRUE)
    res$maf <- out$mAF
    res$communality <- out$h
    aux <- sdevList(res$loadings,dim(res$loadings)[2],maf=res$maf,minmaf=minmaf)
    res$neutral_sdev <- aux$neutral_sdev
    res$p <- aux$p
    res$q <- aux$q
    res$proportion_removed <- aux$prop_removed
    auxpval <- pval(res,communality_test=communality_test)
    auxpval[res$maf<minmaf,] <- NA
    res$pvalues <- auxpval
  }
  resf <- list(scores=res$scores,loadings=res$loadings,correlations=res$correlations,singular_values=res$singular_values,pvalues=res$pvalues,maf=res$maf,communality=res$communality,
              neutral_sdev=res$neutral_sdev,p=res$p,q=res$q,proportion_removed=res$proportion_removed)
  class(resf) <- 'pcadapt'  
  attr(resf,"communality_test") <- communality_test
  attr(resf,"minmaf") <- minmaf
  return(resf)
}

#' Principal Components Analysis
#'
#' \code{corpca} performs a principal components analysis on a dataset, and returns an object \code{x} 
#' which contains the loadings, the scores and the singular values of the \code{K} first principal components. 
#' It handles missing values in a dataset and actually computes the eigen elements of the rather small \code{n x n} 
#' covariance matrix. All these variables are accesible using
#' the dollar sign (example : \code{stds <- x$neutral_sdev}).
#' 
#' @param data a data matrix or a data frame. 
#' @param K an integer specifying the number of principal components that are retained.
#' 
#' @importFrom stats cov
#' 
#' @examples
#' x <- read4pcadapt("geno3pops",option="example")
#' y <- pcadapt(x,K=3)
#' 
#' @export
corpca = function(data,K){
  n <- dim(data)[1]
  p <- dim(data)[2]
  data_aux <- scale(data,scale=TRUE)*sqrt(p/(n-1))
  covmat <- cov(t(data_aux),use="pairwise.complete.obs")
  res <- NULL
  aux <- eigen(covmat,symmetric=TRUE)
  res$sdev <- aux$values[1:K]
  res$scores <- aux$vectors[,1:K]
  aux_ldgs <- t(aux$vectors)%*%data_aux
  res$loadings <- t((1/(sqrt(res$sdev)))*aux_ldgs[1:K,])
  return(res)
}

#' Principal Components p-values
#'
#' \code{pval} computes p-values for each principal component or for a set of principal components in the meta-analysis case.
#' 
#' @param x an output of \code{pcadapt} containing the loadgins and the estimated standard deviations.
#' @param communality_test a logical value indicating whether a communality test should be performed. Default value set to \code{FALSE}.
#' 
#' @examples
#' x <- replicate(1000, rnorm(100))
#' x <- floor(abs(x))
#' y <- pcadapt(x,K=3)
#' pval(y)
#' 
#' @keywords internal
#' 
#' @importFrom stats pnorm median pchisq qchisq
#' 
#' @export
pval = function(x,communality_test=FALSE){
  pvalues <- NULL
  z <- x$loadings
  sigma <- x$neutral_sdev
  K <- length(sigma)
  column_names <- NULL
  if (K==1){
    pvalues <- 2*pnorm(abs(z),sd=sigma[1],lower.tail=FALSE) 
    column_names <- "p1"
    q <- as.data.frame(t(pvalues))
    colnames(q) <- column_names
  } else {
    for (k in 1:K){
      pvalues <- cbind(pvalues,2*pnorm(abs(z[,k]),sd=sigma[k],lower.tail=FALSE))  
      column_names <- c(column_names,paste0("p",k))
      q <- as.data.frame(pvalues)
      colnames(q) <- column_names
    }
  }
  if (communality_test == TRUE){
    h2 <- x$communality
    nSNP <- length(h2)
    df <- K
    c <- sum(x$singular_values[1:df]^2)/df
    gif <- median(h2*nSNP/c,na.rm=TRUE)/qchisq(0.5,df=df)
    pvalues <- pchisq(h2*nSNP/(c*gif),df=df,lower.tail=FALSE)
    column_names <- "p1"
    q <- as.data.frame(pvalues)
    colnames(q) <- column_names
  }
  return(q)
}

#' Analysis of the computed p-values
#'
#' \code{outlier} returns an object \code{x} containing a list of indices (\code{x$outliers}) corresponding to the SNPs ranked by the chosen p-values 
#' associated with \code{K}-th principal component or obtained with the meta-analysis (\code{x$pvalues}), as well as a list of false discovery rates
#' computed for different thresholds (which are set in the built-in function \code{fdrList}) in the case of a simulation. 
#' 
#' @param y an object of class \code{pcadapt} generated with \code{pcadapt}.
#' @param K an integer specifying the principal component the user is interested in. Default value set to \code{1}. NB : this field can be left empty in the
#' meta-analysis case. 
#' @param threshold a real value between \code{0} and \code{1} indicating the proportion of false discoveries to be tolerated.
#' @param list_selection a vector of indices containing the list of all the SNPs under selection in the case of 
#' a simulation. When this argument is provided, a list of false discovery rates is calculated for
#' different thresholds. Default value is set to \code{NULL}.
#' 
#' @examples
#' x <- read4pcadapt("geno3pops",option="example")
#' y <- pcadapt(x,K=3)
#' outlier(y,K=1)
#' 
#' @importFrom stats p.adjust
#'
#' @export
outlier = function(y,K=NULL,threshold=0.10,list_selection = NULL){
  x <- NULL
  minmaf <- attr(y,"minmaf")
  idxmaf <- y$maf>=minmaf
  if (attr(y,"communality_test")==TRUE){
    pval <- y$pvalues[idxmaf,1]
    x$rankedSNPs <- order(pval)
    x$adjusted_pval <- p.adjust(sort(pval),method="BH",n=length(pval))
    if (!is.null(list_selection)){
      x$fdr <- fdrList(x$rankedSNPs,x$adjusted_pval,list_selection)
    }
    num_out <- sum(x$adjusted_pval<threshold)
    if (num_out > 0){
      x$outliers <- (x$rankedSNPs)[1:num_out]
    } else  {
      x$outliers <- 0
    }
    x$pvalues <- y$pvalues
    class(x) <- 'outlier'
  } else {
    if (is.null(K)){
      warning("K has to be specified")
    } else {
      pval <- y$pvalues[idxmaf,K]
      x$rankedSNPs <- order(pval)
      x$adjusted_pval <- p.adjust(sort(pval),method="BH",n=length(pval))
      if (!is.null(list_selection)){
        x$fdr <- fdrList(x$rankedSNPs,x$adjusted_pval,list_selection)
      }
      num_out <- sum(x$adjusted_pval<threshold)
      if (num_out > 0){
        x$outliers <- (x$rankedSNPs)[1:num_out]
      } else  {
        x$outliers <- 0
      }
      x$pvalues <- y$pvalues
      class(x) <- 'outlier'
    }
  }  
  return(x)
}

#' Kurtosis-based Standard Deviation Estimation
#'
#' \code{sdevEstimation} returns an estimate of the standard deviation of a set of values \code{x} which is supposed to follow a normal distribution. 
#' 
#' @param x a numeric vector whose distribution is supposed to be normal.
#' @param method a character string that specifies which method is used to estimate the standard deviations. Default value set to \code{kurtosis.null}.
#' An alternative option is \code{kurtosis.truncated}.
#' 
#' @examples
#' x <- rnorm(10000,sd=0.9)
#' sdevEstimation(x)
#' 
#' @importFrom stats rnorm quantile dnorm pnorm
#' 
#' @export
sdevEstimation = function(x,method="kurtosis.null"){
  x_rm <- x[which(!is.na(x))]
  vecp <- 10^seq(-log10(length(x_rm)),log10(0.2),length.out=100)
  resq <- NULL
  if (method=="kurtosis.null"){
    for (j in vecp){
      resq <- c(resq,kurtosis(x_rm,j))
    }
  } else if (method=="kurtosis.truncated"){
    for (j in vecp){
      gauss <- rnorm(100000)
      resq <- c(resq,kurtosis(x_rm,j)-kurtosis(gauss,j))
    }
  }    
  index_min_kurt <- which.min(abs(resq))
  best_p <- vecp[index_min_kurt]
  best_q <- quantile(abs(x_rm),1-best_p,na.rm=TRUE)
  aux <- x_rm[abs(x_rm)<best_q]  
  n = length(aux)
  sigma_likelihood_fn = function(x){
    sum(dnorm(aux,sd=x,log=TRUE)) - n*(log(pnorm(best_q,sd=x) - pnorm(-best_q,sd=x)))
  }
  resfn <- sapply(sigma<-seq(0.5,3,length=100),sigma_likelihood_fn)
  out <- NULL
  out$vecp <- vecp
  out$resq <- resq
  out$neutral_sdev <- sigma[which.max(resfn)]
  out$prop_removed <- best_p
  return(out)
}

