#' Read files for pcadapt
#'
#' \code{read.genotype} returns the genotype matrix for different types of files.
#' 
#' @param file the name of the file without the extention which the data are to be read from.
#' @param option a character string that specifies which type of data the user is willing to read.
#' @param header a logical value indicating whether the file contains the names of the variables as its first line.
#' @param sep the field separator character.
#' @param mapheader a logical value indicating whether the map file contains the names of the variables as its first line.
#' @param mapsep the field separator character for the map file.
#' 
#' @keywords internal
#' 
#' @importFrom utils read.table
#' 
#' @export
read.genotype = function(file,option=NULL,header=FALSE,sep=" ",mapheader,mapsep=NULL){
  if (option=="ped"){
    dat <- read.table(paste0(file,".ped"),sep=sep,header=header)
    rs <- read.table(paste0(file,".map"),sep=mapsep,header=mapheader)
    genotype <- dat[,-(1:6)]
    n <- dim(genotype)[1]
    p <- dim(genotype)[2]/2
    c <- array(0,dim=p)
    a <- array(0,dim=p)
    t <- array(0,dim=p)
    g <- array(0,dim=p)
    geno <- array(0,c(n,p))
    for (k in 1:p){
      c[k] <- sum(genotype[,2*k-1]=="C") + sum(genotype[,2*k]=="C")
      t[k] <- sum(genotype[,2*k-1]=="T") + sum(genotype[,2*k]=="T")
      a[k] <- sum(genotype[,2*k-1]=="A") + sum(genotype[,2*k]=="A")
      g[k] <- sum(genotype[,2*k-1]=="G") + sum(genotype[,2*k]=="G")
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
    freq <- apply(geno,2,sum)/(2*n)
    minfreq <- pmin(freq,1-freq)
    geno <- geno[,minfreq>0]
    output <- as.data.frame(geno)
    colnames(output) <- as.character(rs[minfreq>0,2])
    attr(output,"pop") <- as.character(dat[,1])
    return(output)
  }
}

#' pcadapt read function
#' 
#' \code{read4pcadapt} is a function which output may be read by the main function \code{pcadapt}. It currently supports two formats.
#' The first one is a text file containing the genotype matrix with individuals in rows and genotype markers in columns (same format as for the \code{R} package LEA).
#' The second supported format is the .ped format.
#' 
#' @param x can be a file path or a file name, depending on the chosen option. 
#' @param option a character string that specifies the type of data to be read: "genotype" (default), "ped" or "pool".
#' @param header a logical value indicating whether the file contains the names of the variables as its first line.
#' @param sep the field separator character.
#' @param transpose a logical value indicating whether the pooled data has the right format, which is
#' n rows and p columns where n is the number of populations and p is the number of genetic markers.
#' @param mapheader a logical value indicating whether the map file contains the names of the variables as its first line.
#' @param mapsep a character string that specifies the type of separator used when reading the data.
#' 
#' @examples
#' ## see also ?pcadapt for examples
#' 
#' @importFrom utils read.table
#' 
#' @export
read4pcadapt = function(x=NULL,option="genotype",header=FALSE,sep="",transpose=FALSE,mapheader=FALSE,mapsep="\t"){
  if (option == "example"){
    out <- read.table(system.file("extdata",x,package="pcadapt"),header=FALSE)
  } else if (option == "genotype" | option == "pool"){
    out <- read.table(x,header=header,sep=sep,na.strings=c("-9","NA","NaN"))
    if (transpose==TRUE){
      out <- t(out)
    }
  } else if (option == "ped"){
    out <- read.genotype(file=x,option="ped",header=header,sep=sep,mapheader=mapheader,mapsep=mapsep)
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

