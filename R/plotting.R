#' pcadapt visualization tool
#'
#' \code{plot.pcadapt} is a method designed for objects of class \code{pcadapt}.
#' It provides a plotting utile for quick visualization of a \code{pcadapt} object.
#' Different options are available : \code{"screeplot"}, \code{"scores"}, \code{"neutral"}, \code{kurtosis},
#' \code{"manhattan"} and \code{"qqplot"}.
#' \code{"screeplot"} shows the decay of the singular values of the genotype matrix and provides
#' a figure which may help in the choice of \code{K}.
#' \code{"scores"} plots the projection of the individuals onto the first two principal components.
#' \code{"neutral"} displays the loadings histogram of the principal component of interest, as well as
#' the estimated distribution of the loadings of the neutral SNPs.
#' \code{"kurtosis"} plots the kurtosis as a function of the proportion of observations removed from the 
#' distribution of interest. 
#' \code{"manhattan"} draws the Manhattan plot of the p-values associated with the principal component
#' of interest.
#' \code{"qqplot"} draws a Q-Q plot of the p-values associated with the principal component
#' of interest.
#' 
#' @param x an object of class "pcadapt" generated with \code{pcadapt}. 
#' @param ... \dots
#' @param option a character string specifying the figures to be displayed. If \code{NULL} (the default), all three plots are printed.
#' @param K an integer specifying the principal component of interest. \code{K} has to be specified uniquely when using the \code{neutral} option.
#' @param i an integer indicating onto which principal component the individuals are projected when the "scores" option is chosen. 
#' Default value is set to \code{1}.
#' @param j an integer indicating onto which principal component the individuals are projected when the "scores" option is chosen. 
#' Default value is set to \code{2}.
#' @param subpop a list of lists of indices specifying the subpopulations. This attributes specific colors to each subpopulations.
#' @param subcol a list of colors for each subpopulation. If not specified, default colors will be used.
#' @param threshold a parameter which depends on the chosen option. When used with \code{"scores"}, it's the threshold which the user expects the false discovery rate 
#' to be under. For the \code{"qqplot"} option, it displays an addiotional bar which separates the SNPs with the highest p-values from the ones with
#' lower p-values.
#' @param num_pc an integer specifying the number of components to take into account in the scree plot. \code{num_axis} should be lower than \code{K}.
#' By default, \code{num_axis = K}.
#' 
#' @examples
#' ## see ?pcadapt for examples
#'
#' @method plot pcadapt
#' 
#' @importFrom graphics plot title axis abline legend 
#' 
#' @export
plot.pcadapt = function(x,...,option=NULL,K=NULL,i=1,j=2,subpop=NULL,subcol=NULL,threshold=NULL,num_pc=NULL){
  if (is.null(option)){
    screePlot(x,num_pc = num_pc)
    scoresPlot(x,i,j,subpop)
  } else if (!is.null(option)){
    if (option == "screeplot"){
      screePlot(x,num_pc = num_pc)
    } else if (option == "scores"){
      scoresPlot(x,i,j,subpop,subcol)
    } else if (option == "neutral"){
      if (is.null(K)){
        warning("K has to be specified")
      } else{
        neutralDistribution(x,K)
      }
    } else if (option == "kurtosis"){
      if (is.null(K)){
        warning("K has to be specified")
      } else{
        plot(x$p,x$q[,K],xlab="Proportion removed from the distribution",ylab="Excess kurtosis",xaxt='n')
        title(paste0("Excess kurtosis - PC",K))
        tick_marks_x <- seq(0,0.2,by=0.025)
        tick_marks_y <- seq(floor(min(x$q[,K])),floor(max(x$q[,K])),by=0.5)
        axis(1,at=tick_marks_x)
        axis(2,at=tick_marks_y)
        abline(v=tick_marks_x,col="black",lty="dotted")
        abline(h=tick_marks_y,col="black",lty="dotted")
        abline(v=x$proportion_removed[K],col="red",lwd=2)
        legend("topright",legend=paste0("p = ",round(x$proportion_removed[K]*100,digits = 2),"%"),col="red",lty=1,lwd=3)
      }
    } else if (option == "manhattan"){
      if (attr(x,"communality_test")==FALSE){
        if (is.null(K)){
          warning("K has to be specified")
        } else{
          manhattanPlot(x,K,alpha=threshold)
        }
      } else {
        manhattanPlot(x,K=1,alpha=threshold)        
      }
    } else if (option == "qqplot"){
      if (attr(x,"communality_test")==FALSE){
        if (is.null(K)){
          warning("K has to be specified")
        } else{
          pvalqqPlot(x,K,threshold=threshold)
        }
      } else {
        pvalqqPlot(x,K=1,threshold=threshold)        
      }
    }
  }
}

#' Outlier visualization tool
#'
#' \code{plot.outlier} is a method designed for objects of class \code{outlier}. For the moment,
#' the only option available is "fdr", which displays the observed false discovery rate as a function
#' of the expected false discovery rate.
#' 
#' @param x an object of class "outlier" generated with \code{outlier} containing the p-values of interest. 
#' @param ... \dots
#' @param option a character string specifying the figures to be displayed. If \code{NULL} (the default), all available plots are printed.
#' 
#' @examples
#' ## see ?outlier for examples
#'
#' @export
plot.outlier = function(x,...,option=NULL){
  if (is.null(x$fdr)){
    warning("x has no attribute fdr. This option is only available for known simulations.")  
  } else {
    fdrControl(x)
  }
}

#' p-values Q-Q Plot
#'
#' \code{pvalqqPlot} plots a Q-Q plot of the p-values computed ).
#' 
#' @param x an output from \code{outlier} containing the p-values of interest. 
#' @param alpha a real number between \code{0} et \code{1}.
#' 
#' @examples
#' ## see ?pcadapt for examples
#' 
#' @keywords internal
#' 
#' @importFrom graphics plot.default abline legend title
#'
#' @export
pvalqqPlot = function(x,K,threshold = NULL){
  minmaf <- attr(x,"minmaf")
  idxmaf <- x$maf>=minmaf
  sorted_pval <- sort(x$pvalues[idxmaf,K])
  p <- length(sorted_pval) 
  expected_p <- 1:p/p
  plot.default(-log10(expected_p),-log10(sorted_pval),xlab="Expected -log10(p-values)",ylab="Observed -log10(p-values)",cex=0.3,pch=19)
  abline(0,1,col="red",lwd=3)
  if (!is.null(threshold)){
    q <- floor(threshold*p)
    pval_thresh <- expected_p[q]
    abline(v=-log10(pval_thresh),col="blue",lty=3,lwd=3)
    legend("topleft",col="blue",lty=3,lwd=3,legend=paste0("Top ",100*threshold,"%"),cex=0.75)
  }
  title("p-values Q-Q plot")
}

#' Principal Components Analysis Scores Plot
#'
#' \code{"scoresPlot"} plots the projection of the individuals onto the first two principal components.
#' 
#' @param x an output from \code{pcadapt} containing the scores.
#' @param i an integer indicating onto which principal component the individuals are projected when the "scores" option is chosen. 
#' Default value is set to \code{1}.
#' @param j an integer indicating onto which principal component the individuals are projected when the "scores" option is chosen. 
#' Default value is set to \code{2}.
#' @param subpop a list of lists of indices specifying the subpopulations. This attributes specific colors to each subpopulations.
#' @param subcol a list of colors for each subpopulation. If not specified, default colors will be used.
#' 
#' @examples
#' ## see ?pcadapt for examples
#'
#' @keywords internal
#'
#' @importFrom graphics plot.default title points
#' @importFrom grDevices rainbow
#'
#' @export
scoresPlot = function(x,i,j,subpop=NULL,subcol=NULL){
  if (dim(x$loadings)[1]<2){
    warning("K=1, option not available since two principal components have to be computed at least")
  } else {
  if (is.null(subpop)){
    plot.default(x$scores[,i],x$scores[,j],cex=0.5,pch=19,xlab=paste0("PC",i),ylab=paste0("PC",j))  
    title(paste0("Projection onto PC",i," and PC",j))
  } else {
    if (is.null(subcol)){
      palette <- rainbow(length(subpop))
    } else {
      palette <- subcol
    }
    x0<-min(x$scores[,i])
    x1<-max(x$scores[,i])
    y0<-min(x$scores[,j])
    y1<-max(x$scores[,j])
    plot.default(x$scores[unlist(subpop[1]),i],x$scores[unlist(subpop[1]),j],cex=0.5,pch=19,xlab=paste0("PC",i),ylab=paste0("PC",j),col=palette[1],xlim=c(x0,x1),ylim=c(y0,y1))
    for (k in 2:length(subpop)){
      points(x$scores[unlist(subpop[k]),i],x$scores[unlist(subpop[k]),j],cex=0.5,pch=19,col=palette[k])
    }
    title(paste0("Projections onto PC",i," and PC",j))
  }
  }
}

#' Principal Components Analysis Scree Plot
#'
#' \code{screePlot} plots the scee plot associated with the principal components analysis performed on the dataset. 
#' NB : \code{pcadapt} has to be run on the dataset in order to get an output readable by \code{plot.screePlot}
#' 
#' @param x an output from \code{pcadapt} containing the singular values. 
#' @param num_pc an integer specifying the number of components to take into account in the scree plot. \code{num_axis} should be lower than \code{K}.
#' By default, \code{num_axis = K}.
#' 
#' @examples
#' ## see ?pcadapt for examples
#'
#' @keywords internal
#'
#' @importFrom graphics plot.default axis title
#'
#' @export
screePlot = function(x,num_pc=NULL){
  if (is.null(num_pc)){
    K <- length(x$singular_values)
  } else {
    K <- num_pc
  }
  
  if (K<2){
    warning("K = 1, the scree plot is thus composed of a unique point")
  }
  plot.default(1:K,(x$singular_values[1:K])^2/length(x$communality),xlab="PCs",ylab="Proportion of explained variance",type="b",pch=19,lwd=3,col="red",xaxt="n")
  axis(1,at=seq(1,K,by=1))
  title(paste("Scree Plot - K =",K))
}

#' Manhattan Plot
#'
#' \code{manhattanPlot} displays a Manhattan plot which represents the p-values for each SNP for a particular principal component (or set of
#' principal components in the communality test case).
#' 
#' @param x an object of class "outlier" generated with \code{outlier} containing the p-values of interest. 
#' @param alpha a real number between \code{0} et \code{1}.
#' 
#' @examples
#' ## see ?pcadapt for examples
#' 
#' @keywords internal
#' 
#' @importFrom stats p.adjust 
#' @importFrom graphics plot.default title abline legend
#'
#' @export
manhattanPlot = function(x,K,alpha){ 
  minmaf <- attr(x,"minmaf")
  pval_K <- x$pvalues[x$maf>=minmaf,K]
  p <- length(pval_K)  
  plot.default(1:p,-log10(pval_K),xlab="SNP",ylab="-log10(p-values)",cex=0.3,pch=19)
  title("Manhattan plot")
  if (!is.null(alpha)){
    adjpval <- p.adjust(sort(pval_K),method="BH",n=p)
    idx_pmax <- sum(adjpval<=alpha)
    if (idx_pmax > 0){
      pval_alpha = adjpval[idx_pmax]*idx_pmax/p
      abline(h = -log10(pval_alpha),col="red",lwd=3)
      legend("topright",legend=paste("E[FDR]<",alpha),col="red",lty=1,lwd=3)
    }
  }
}

#' Neutral Distribution Estimation
#'
#' \code{neutralDistribution} plots the histogram of the loadings associated with the \code{K}-th principal component, as well
#' as the normal distribution with the estimated standard deviation.
#' 
#' @param data a numeric matrix which provides the data for the principal components analysis (PCA). 
#' @param K an integer indicating which principal component the histogram will be associated with.
#' @param sigma a numeric vector containing the estimated standard deviations for each principal component.
#' @param n_breaks an integer indicating the number of bars the histogram will be composed of. Default value set to \code{100}.
#' 
#' @examples
#' ## see ?pcadapt for examples
#' 
#' @keywords internal
#' 
#' @importFrom stats dnorm
#' @importFrom graphics hist lines title legend 
#'
#' @export
neutralDistribution = function(x,K,n_breaks=100){
  minmaf <- attr(x,"minmaf")
  idxmaf <- x$maf>minmaf
  if (dim(x$loadings)[1]>1){
    z <- x$loadings[idxmaf,K]
  } else {
    z <- x$loadings[idxmaf]
  }
  sigma <- x$neutral_sdev[K]
  h <- hist(z,breaks=n_breaks,freq=FALSE,xlab="Loadings",main=NULL)
  x <- seq(floor(min(z[which(!is.na(z))])),floor(max(z[which(!is.na(z))])+1),length=200)
  y <- dnorm(x,mean=0,sd=sigma)
  lines(x,y,type="l",lwd=4,col="red")
  title(paste0("Loadings distribution - PC",K))
  legend("topright",legend=paste("estimated sdev =",(round(sigma,digits=5))),cex=0.75)
}

