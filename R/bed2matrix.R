#' Convert a bed to a matrix
#'
#' @param bedfile Path to a bed file. 
#' @param n Number of samples. Default reads it from corresponding fam file.
#' @param p Number of SNPs. Default reads it from corresponding bim file.
#'
#' @return An integer matrix.
#' @export
#'
#' @examples
#' bedfile <- system.file("extdata", "geno3pops.bed", package = "pcadapt")
#' mat <- bed2matrix(bedfile)
#' dim(mat)
#' table(mat)
bed2matrix <- function(bedfile, n = NULL, p = NULL) {
  
  if (inherits(bedfile, "pcadapt_bed")) {
    n <- attr(bedfile, "n")
    p <- attr(bedfile, "p")
  }
  
  if (is.null(n) || is.null(p)) {
    bed <- read.pcadapt(bedfile, type = "bed")
    n <- attr(bed, "n")
    p <- attr(bed, "p")
  } 
  
  xptr <- bedXPtr(bedfile, n, p)
  bed2mat(xptr)
}
