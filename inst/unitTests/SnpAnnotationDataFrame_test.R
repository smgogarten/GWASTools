test_SnpAnnotationDataFrame <- function() {
  # test object
  snpID <- 1:10
  chrom <- c(rep(1L,5), 23:27)
  pos <- 101:110
  rsID <- paste("rs", 1:10, sep="")
  A <- rep("A",10)
  B <- rep("B",10)
  x <- data.frame(snpID=snpID, chromosome=chrom, position=pos, rsID=rsID,
                  alleleA=A, alleleB=B, stringsAsFactors=FALSE)
  obj <- SnpAnnotationDataFrame(x)
  checkIdentical(x, as(obj, "data.frame"))

  # retrieve named columns
  checkIdentical(snpID, getSnpID(obj))
  checkIdentical(chrom, getChromosome(obj))
  checkIdentical(c(rep("1",5),"X","XY","Y","M","U"),
                 getChromosome(obj, char=TRUE))
  checkIdentical(pos, getPosition(obj))
  checkIdentical(A, getAlleleA(obj))
  checkIdentical(B, getAlleleB(obj))

  # other columns
  checkTrue(hasVariable(obj, "rsID"))
  checkIdentical(rsID, getVariable(obj,"rsID"))
  checkTrue(!hasVariable(obj, "foo"))
  checkIdentical(NULL, getVariable(obj, "foo"))
  vars <- c("snpID", "rsID")
  checkIdentical(x[,vars], getVariable(obj, vars))

  # check indexing
  sel <- 1:5 # numeric
  checkIdentical(snpID[sel], getSnpID(obj, sel))
  sel <- snpID > 5 # logical
  checkIdentical(snpID[sel], getSnpID(obj, sel))

  # missing required columns
  x <- data.frame(snpID=snpID, chromosome=chrom)
  checkException(SnpAnnotationDataFrame(x))
  x <- data.frame(snpID=snpID, position=pos)
  checkException(SnpAnnotationDataFrame(x))
  x <- data.frame(chromosome=chrom, position=pos)
  checkException(SnpAnnotationDataFrame(x))

  # bad column format
  x <- data.frame(snpID=as.character(snpID), chromosome=chrom, position=pos)
  checkException(SnpAnnotationDataFrame(x))
  x <- data.frame(snpID=snpID, chromosome=as.character(chrom), position=pos)
  checkException(SnpAnnotationDataFrame(x))
  x <- data.frame(snpID=snpID, chromosome=chrom, position=as.character(pos))
  checkException(SnpAnnotationDataFrame(x))

  # snpID not unique
  x <- data.frame(snpID=rep(1L,10), chromosome=chrom, position=pos)
  checkException(SnpAnnotationDataFrame(x))

  # get and set methods
  x <- data.frame(snpID=snpID, chromosome=chrom, position=pos)
  pData(obj) <- x
  checkIdentical(x, getAnnotation(obj))
  meta <- getMetadata(obj)
  meta["snpID", "labelDescription"] <- "id"
  varMetadata(obj) <- meta
  checkIdentical(meta, getMetadata(obj))

  # check alternate chromosome codes
  x <- data.frame(snpID=snpID, chromosome=chrom, position=pos)
  obj <- SnpAnnotationDataFrame(x, YchromCode=24L, XYchromCode=25L)
  checkIdentical(c(rep("1",5),"X","Y","XY","M","U"),
                 getChromosome(obj, char=TRUE))
}
