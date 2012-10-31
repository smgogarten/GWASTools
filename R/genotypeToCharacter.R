genotypeToCharacter <- function(geno, alleleA=NULL, alleleB=NULL, sort=TRUE) {
  if (length(dim(geno)) < 2) {
    geno <- matrix(geno, ncol=1)
  }
  if (is.null(alleleA) | is.null(alleleB)) {
    geno[geno %in% 0] <- "B/B"
    geno[geno %in% 1] <- "A/B"
    geno[geno %in% 2] <- "A/A"                
  } else {
    alleleA <- as.character(alleleA)
    alleleB <- as.character(alleleB)
    aa <- paste(alleleA, alleleA, sep="/")
    bb <- paste(alleleB, alleleB, sep="/")
    if (sort) {
      ab <- paste(pmin(alleleA, alleleB), pmax(alleleA, alleleB), sep="/")
    } else {
      ab <- paste(alleleA, alleleB, sep="/")
    }
    for (k in 1:ncol(geno)) {
      geno[geno[,k] %in% 0, k] <- bb[geno[,k] %in% 0]
      geno[geno[,k] %in% 1, k] <- ab[geno[,k] %in% 1]
      geno[geno[,k] %in% 2, k] <- aa[geno[,k] %in% 2]
    }
  }
  if (ncol(geno) == 1) {
    geno <- as.vector(geno)
  }
  return(geno)
}
