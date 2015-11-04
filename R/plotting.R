#' pcadapt visualization tool
#'
#' \code{plot.pcadapt} is a method designed for objects of class \code{pcadapt}.
#' It provides a plotting utile for quick visualization of a \code{pcadapt} object.
#' Different options are available : \code{"screeplot"}, \code{"scores"}, \code{"stat.distribution"},
#' \code{"manhattan"} and \code{"qqplot"}.
#' \code{"screeplot"} shows the decay of the singular values of the genotype matrix and provides
#' a figure to guide in the choice of \code{K}.
#' \code{"scores"} plots the projection of the individuals onto the first two principal components.
#' \code{"stat.distribution"} displays the histogram of the selected test statistics, as well as
#' the estimated distribution for the neutral SNPs.
#' \code{"manhattan"} draws the Manhattan plot of the p-values associated with the principal component
#' of interest.
#' \code{"qqplot"} draws a Q-Q plot of the p-values associated with the principal component
#' of interest.
#' 
#' @param x an object of class "pcadapt" generated with \code{pcadapt}. 
#' @param ... \dots
#' @param option a character string specifying the figures to be displayed. If \code{NULL} (the default), all three plots are printed.
#' @param K an integer specifying the principal component of interest. \code{K} has to be specified only when using the \code{loadings} option.
#' @param i an integer indicating onto which principal component the individuals are projected when the "scores" option is chosen. 
#' Default value is set to \code{1}.
#' @param j an integer indicating onto which principal component the individuals are projected when the "scores" option is chosen. 
#' Default value is set to \code{2}.
#' @param pop a list of lists of indices specifying the subpopulations. This attributes specific colors to each subpopulations.
#' @param subcol a list of colors for each subpopulation. If not specified, default colors will be used.
#' @param threshold for the \code{"qqplot"} option, it displays an additional bar which shows the \code{threshold} percent of SNPs with smallest p-valuesseparates the SNPs with the highest p-values.
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
plot.pcadapt = function(x,...,option=NULL,K=NULL,i=1,j=2,pop=attr(x,"pop"),subcol=NULL,threshold=NULL,num_pc=NULL){
  if (!(option %in% c("screeplot","scores","manhattan","qqplot","stat.distribution"))){
    warning(paste("Plotting option",option,"not valid, options currently available are: screeplot, scores, manhattan, qqplot, loadings."))
  } else {
    if (option == "screeplot"){
      screePlot(x,num_pc = num_pc)   
    } else if (option == "scores"){
      if (attr(x,"data.type")!="pool"){
        scoresPlot(x,i,j,pop,subcol)
      } else {
        scoresPlot(x,i,j,pop=1:dim(x$scores)[1],subcol)
      }   
    } else if (option == "stat.distribution"){
      if ((attr(x,"method") %in% c("mahalanobis","euclidean","communality")) == FALSE){
        if (is.null(K)){
          warning("K has to be specified")
        } else {
          neutralDistribution(x,K)
        }
      } else {
        neutralDistribution(x,1)
      }
    } else if (option == "manhattan"){
      if ((attr(x,"method") %in% c("mahalanobis","euclidean","communality")) == FALSE){
        if (is.null(K)){
          warning("K has to be specified")
        } else {
          manhattanPlot(x,K)
        }
      } else {
        manhattanPlot(x,K=1)        
      }
    } else if (option == "qqplot"){
      if ((attr(x,"method") %in% c("mahalanobis","euclidean","communality")) == FALSE){
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
  if (attr(x,"method")=="componentwise"){
    sorted_pval <- sort(x$pvalues[x$maf>=attr(x,"minmaf"),K])
  } else {
    sorted_pval <- sort(x$pvalues[x$maf>=attr(x,"minmaf")])
  }
  p <- length(sorted_pval) 
  expected_p <- 1:p/p
  plot.default(-log10(expected_p),-log10(sorted_pval),xlab="Expected -log10(p-values)",ylab="Observed -log10(p-values)",cex=0.3,pch=19)
  abline(0,1,col="red",lwd=3)
  if (!is.null(threshold)){
    q <- floor(threshold*p)
    pval_thresh <- expected_p[q]
    abline(v=-log10(pval_thresh),col="grey",lty=3,lwd=3)
    legend("topleft",col="grey",lty=3,lwd=3,legend=paste0("Top ",100*threshold,"%"),cex=0.75)
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
scoresPlot = function(x,i,j,pop,subcol=NULL){
  if (dim(x$loadings)[1]<2){
    warning("K=1, option not available since two principal components have to be computed at least")
  } else {
  if (is.null(pop)){
    plot.default(x$scores[,i],x$scores[,j],cex=0.5,pch=19,xlab=paste0("PC",i),ylab=paste0("PC",j))   
    title(paste0("Projection onto PC",i," and PC",j))    
  } else {
    pop_aux <- getPopColors(pop)
    if (is.null(subcol)){
      palette <- rainbow(length(pop_aux))
    } else {
      palette <- subcol
    }
    x0<-min(x$scores[,i])
    x1<-max(x$scores[,i])
    y0<-min(x$scores[,j])
    y1<-max(x$scores[,j])
    plot.default(x$scores[unlist(pop_aux[1]),i],x$scores[unlist(pop_aux[1]),j],cex=0.5,pch=19,xlab=paste0("PC",i),ylab=paste0("PC",j),col=palette[1],xlim=c(x0,x1),ylim=c(y0,y1))
    for (k in 2:length(pop_aux)){
      points(x$scores[unlist(pop_aux[k]),i],x$scores[unlist(pop_aux[k]),j],cex=0.5,pch=19,col=palette[k])
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
  if (is.null(num_pc)){K <- length(x$singular_values)} 
  else {K <- num_pc}  
  if (K<2){warning("K = 1, the scree plot is thus composed of a unique point")}
  if (attr(x,"method")=="componentwise"){
    nSNP <- length(x$pvalues[,1])
  } else {
    nSNP <- length(x$pvalues)
  }
  plot.default(1:K,(x$singular_values[1:K])^2/nSNP,xlab="PCs",ylab="Proportion of explained variance",type="b",pch=19,lwd=3,col="red",xaxt="n")
  axis(1,at=seq(1,K,by=1))
  title(paste("Scree Plot - K =",K))
}

#' Manhattan Plot
#'
#' \code{manhattanPlot} displays a Manhattan plot which represents the p-values for each SNP for a particular test statistic.
#' 
#' @param x an object of class "pcadapt" generated with \code{pcadapt} containing the p-values of interest.
#' 
#' @examples
#' ## see ?pcadapt for examples
#' 
#' @keywords internal
#' 
#' @importFrom graphics plot.default title
#'
#' @export
manhattanPlot = function(x,K){ 
  if (attr(x,"method")=="componentwise"){
    pval_K <- x$pvalues[x$maf>=attr(x,"minmaf"),K]
  } else {
    pval_K <- x$pvalues[x$maf>=attr(x,"minmaf")]
  }
  p <- length(pval_K)  
  plot.default(1:p,-log10(pval_K),xlab="SNP",ylab="-log10(p-values)",cex=0.3,pch=19)
  title("Manhattan plot")
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
#' @importFrom stats dnorm dchisq mad
#' @importFrom graphics hist lines title legend 
#'
#' @export
neutralDistribution = function(x,K,n_breaks=100){
  idxmaf <- x$maf>attr(x,"minmaf")
  if (attr(x,"method")=="componentwise"){
    z <- x$loadings[idxmaf,K]/x$gif[K]
    sigma <- apply(x$loadings,2,mad)
    h <- hist(z,breaks=n_breaks,freq=FALSE,xlab="Loadings",main=NULL)
    x1 <- seq(floor(min(z[which(!is.na(z))])),floor(max(z[which(!is.na(z))])+1),length=200)
    y <- dnorm(x1,mean=0,sd=sigma)
    lines(x1,y,type="l",lwd=4,col="red")
    title(paste0("Loadings distribution - PC",K))
    legend("topright",legend="null distribution",lwd=2,col="red")
  } else if (attr(x,"method")=="mahalanobis" && attr(x,"data.type")!="pool"){
    z <- x$stat[idxmaf]/x$gif
    df <- attr(x,"K")
    h <- hist(z,breaks=n_breaks,freq=FALSE,xlab="Scaled Mahalanobis distance",main=NULL)
    x1 <- seq(floor(min(z[which(!is.na(z))])),floor(max(z[which(!is.na(z))])+1),length=200)
    y <- dchisq(x1,df=df)
    lines(x1,y,type="l",lwd=4,col="red")
    title("Distance distribution")
    legend("topright",legend="null distribution",lwd=2,col="red")
  } else if (attr(x,"method")=="communality"){
    df <- attr(x,"K")
    c <- sum(x$singular_values[1:df]^2)/df
    z <- x$stat[idxmaf]*length(x$stat)/(x$gif*c)
    h <- hist(z,breaks=n_breaks,freq=FALSE,xlab="Scaled Communality",main=NULL)
    x1 <- seq(floor(min(z[which(!is.na(z))])),floor(max(z[which(!is.na(z))])+1),length=200)
    y <- dchisq(x1,df=df)
    lines(x1,y,type="l",lwd=4,col="red")
    title("Communality distribution")
    legend("topright",legend="null distribution",lwd=2,col="red")
  } else if (attr(x,"method")=="euclidean"){
    df <- attr(x,"K")
    z <- x$stat[idxmaf]/x$gif
    h <- hist(z,breaks=n_breaks,freq=FALSE,xlab="Scaled Euclidean distance",main=NULL)
    x1 <- seq(floor(min(z[which(!is.na(z))])),floor(max(z[which(!is.na(z))])+1),length=200)
    y <- dchisq(x1,df=df)
    lines(x1,y,type="l",lwd=4,col="red")
    title("Euclidean distance distribution")
    legend("topright",legend="null distribution",lwd=2,col="red")
  } else if (attr(x,"method")=="mahalanobis" && attr(x,"data.type")=="pool"){
    z <- x$stat/x$gif
    df <- attr(x,"K")
    h <- hist(z,breaks=n_breaks,freq=FALSE,xlab="Scaled Mahalanobis distance",main=NULL)
    x1 <- seq(floor(min(z[which(!is.na(z))])),floor(max(z[which(!is.na(z))])+1),length=200)
    y <- dchisq(x1,df=df)
    lines(x1,y,type="l",lwd=4,col="red")
    title("Distance distribution")
    legend("topright",legend="null distribution",lwd=2,col="red")
  }
}

