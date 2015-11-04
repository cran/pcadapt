#' Fst
#' 
#' \code{fstCalc} returns the list of Fst associated with each SNP using Weir and Cockerham's formula. 
#' 
#' @param data a data matrix or a data frame. 
#' @param POPs a list of integers specifying which subpopulation the individuals belong to.
#' @param nPOP an integer specifying the number of populations in the dataset.
#' @param PopSizes a list of length \code{nPOP} specifying the number of individuals for each subpopulation.
#' @param ploidy an integer specifying the ploidy of the individuals.
#' 
#' @references Weir, B. S., & Cockerham, C. C. (1984). Estimating F-statistics for the analysis of population structure. Evolution, 1358-1370.
#' 
#' @examples
#' x <- read4pcadapt("geno3pops",option="example")
#' popsizes <- c(50,50,50) 
#' pops <- c(rep(1,times=50),rep(2,times=50),rep(3,times=50))
#' nPOP <- 3
#' fst <- fstCalc(x,pops,nPOP,popsizes)
#'  
#' @export
fstCalc = function(data,POPs,nPOP,PopSizes,ploidy=2){
  nIND <- dim(data)[1]
  nSNP <- dim(data)[2]
  FstWC <- array(0,dim=nSNP)
  freq <- apply(data,2,sum)/(ploidy*dim(data)[1])
  minfreq <- pmin(freq,1-freq)
  for (i in 1:nSNP){
    af <- freq[i]
    if ((af>0) && (af<1)){
      Geno <- data[,i]
      FstWC[i] <- fstGeno(Geno,POPs,nIND,nPOP,PopSizes)
      if (FstWC[i] > 1.5 | FstWC[i] < -.5){FstWC[i]<-0}
    } 
  } 
  return(FstWC)
}

#' Fst auxiliary function
#' 
#' \code{fstGeno} returns the Fst associated with a particular SNP using Weir and Cokram's formula 
#' 
#' @param Geno a vector of \code{0,1} or \code{2} of length \code{nIND}.
#' @param POPs a list of integers specifying which subpopulation the individuals belong to.
#' @param nPOP an integer specifying the number of SNPs present in the dataset.
#' @param PopSizes a list of length \code{nPOP} specifying the number of individuals for each subpopulation.
#' 
#' @examples
#' ## see also ?fstCalc for examples
#' 
#' @keywords internal
#'  
#' @export
fstGeno = function(Geno,POPs,nIND,nPOP,PopSizes){
  n_bar <- 0
  p_bar <- 0
  h_bar <- 0
  s2 <- 0
  nIND <- length(Geno)
  PopFreqs <- array(0,dim=nPOP)
  PopHet <- array(0,dim=nPOP)
  for (p in 1:nPOP){
    n_bar <- n_bar + PopSizes[p]/nPOP
    for (i in 1:nIND){
      if (POPs[i] == p){
        PopFreqs[p] <- PopFreqs[p] + Geno[i]/2        
      }
      if ((POPs[i] == p) && (Geno[i] == 1)){
        PopHet[p] <- PopHet[p] + 1
      }
    }
    PopFreqs[p] <- PopFreqs[p]/PopSizes[p]
    PopHet[p] <- PopHet[p]/PopSizes[p]
  }  
  
  for (p in 1:nPOP){
    p_bar <- p_bar + PopSizes[p]*PopFreqs[p]/(n_bar*nPOP)
    h_bar <- h_bar + PopSizes[p]*PopHet[p]/(n_bar*nPOP)
  }
  
  n_c <- nPOP*n_bar
  for (p in 1:nPOP){
    n_c = n_c - PopSizes[p]^2/(nPOP*n_bar)
  }
  n_c <- n_c/(nPOP-1)
  
  for (p in 1:nPOP){
    s2 <- s2 + PopSizes[p]*(PopFreqs[p]-p_bar)^2/(n_bar*(nPOP-1))
  }
  
  h_bar <- 0
  a <- (n_bar/n_c)*(s2 - 1.0/(n_bar - 1.0)*(p_bar*(1.0 - p_bar) - s2*(nPOP - 1)/(nPOP) - h_bar/4.0))
  b <- (n_bar/(n_bar - 1.0))*(p_bar*(1.0 - p_bar) - ((nPOP - 1.0)/nPOP)*s2 - h_bar*(2*n_bar - 1)/(4*n_bar))
  c <- h_bar/2
  
  return(a/(a+b+c))
}