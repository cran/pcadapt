#' Read files for pcadapt
#'
#' \code{read.genotype} returns the genotype matrix for different types of files.
#' 
#' @param file the name of the file without the extention which the data are to be read from.
#' @param option a character string that specifies which type of data the user is willing to read.
#' 
#' @keywords internal
#' 
#' @importFrom utils read.table
#' 
#' @export
read.genotype = function(file,option=NULL){
  if (option=="ped"){
    dat <- read.table(paste0(file,".ped"),sep=" ",header=FALSE)
    rs <- read.table(paste0(file,".map"),sep="\t",header=FALSE)
    genotype <- dat[,-(1:6)]
    n <- dim(genotype)[1]
    p <- dim(genotype)[2]/2
    
    c<-array(0,dim=p)
    a<-array(0,dim=p)
    t<-array(0,dim=p)
    g<-array(0,dim=p)
    geno <- array(0,c(n,p))
    for (k in 1:p){
      c[k]<-sum(genotype[,2*k-1]=="C")+sum(genotype[,2*k]=="C")
      t[k]<-sum(genotype[,2*k-1]=="T")+sum(genotype[,2*k]=="T")
      a[k]<-sum(genotype[,2*k-1]=="A")+sum(genotype[,2*k]=="A")
      g[k]<-sum(genotype[,2*k-1]=="G")+sum(genotype[,2*k]=="G")
      if (c[k]>0){
        if (c[k]>t[k]){
          geno[,k] <- (genotype[,2*k-1]=="T") + (genotype[2*k]=="T")
        } else {
          geno[,k] <- (genotype[,2*k-1]=="C") + (genotype[2*k]=="C")
        }
      }
      if (a[k]>0){
        if (a[k]>g[k]){
          geno[,k] <- (genotype[,2*k-1]=="G") + (genotype[2*k]=="G")
        } else {
          geno[,k] <- (genotype[,2*k-1]=="A") + (genotype[2*k]=="A")
        }
      }
    }
    output <- as.data.frame(geno)
    colnames(output) <- as.character(rs[,2])
    return(output)
  }
}

#' Multi-purpose read function
#'
#' \code{read4pcadatp} allows the user to display individuals of the same pre-defined subpopulation with the same color when using the option 
#' \code{"scores"} in \code{pcadapt}. Creates a global variable storing the lists of indices for each subpopulation.
#' 
#' @param x can be a file path or a file name, depending on the chosen option. 
#' @param option a character string specifying the type of data going to be read.
#' 
#' @examples
#' ## see also ?pcadapt for examples
#' 
#' @importFrom utils read.table
#' 
#' @export
read4pcadapt = function(x=NULL,option="genotype"){
  if (option == "example"){
    out <- read.table(system.file("extdata",x,package="pcadapt"),header=FALSE)
  } else if (option == "genotype"){
    out <- read.table(x,header=FALSE,na.strings=c("-9","NA","NaN"))
  } else if (option == "ped"){
    out <- read.genotype(file=x,option="ped")
  } 
  return(out)
}

#' Subpopulation colorization
#'
#' \code{getPopColors} allows the user to display individuals of the same pre-defined subpopulation with the same color when using the option 
#' \code{"scores"} in \code{pcadapt}.
#' 
#' @param pop a list of integers or strings specifying which subpopulation the individuals belong to.
#' 
#' @examples
#' ## see also ?pcadapt for examples
#' 
#' @keywords internal
#' 
#' @export
getPopColors = function(pop){
  subpop <- list()
  listref <- unlist(pop)
  idx <- 1
  while (length(listref)>0){
    col <- listref[1]    
    subpop[[idx]] <- which(pop==col)
    idx <- idx+1
    listref <- listref[listref != col]
  }
  return(subpop)
}