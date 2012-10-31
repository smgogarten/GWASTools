test_SnpAnnotationSQLite <- function() {
  # create new database
  dbpath <- tempfile()  
  obj <- SnpAnnotationSQLite(dbpath)
  
  # test object
  snpID <- 1:10
  chrom <- c(rep(1L,5), 23:27)
  pos <- 101:110
  x <- data.frame(snpID=snpID, chromosome=chrom, position=pos)
  writeAnnotation(obj, x)
  checkIdentical(x, getAnnotation(obj))
  checkEquals(10, nsnp(obj))
  metadata <- data.frame(varname=c("snpID", "chromosome", "position"),
                         description=c("id", "chrom", "pos"),
                         stringsAsFactors=FALSE)
  writeMetadata(obj, metadata)
  checkIdentical(metadata, getMetadata(obj))

  # test errors
  xx <- data.frame(snpID=snpID, chromosome=chrom)
  checkException(writeAnnotation(obj, xx))
  xmeta <- data.frame(var="snpID", descr="id")
  checkException(writeMetadata(obj, xmeta))

  # add columns
  alleleA <- rep("A",10)
  alleleB <- rep("B",10)
  rsID <- paste("rs", 1:10, sep="")
  x <- cbind(x, alleleA, alleleB, rsID, stringsAsFactors=FALSE)
  writeAnnotation(obj, x)
  checkIdentical(x, getAnnotation(obj))
  checkException(validObject(obj))
  newmeta <- data.frame(varname=c("alleleA", "alleleB", "rsID"),
                        description=c("A", "B", "rs id"),
                        stringsAsFactors=FALSE)
  writeMetadata(obj, newmeta, append=TRUE)
  checkIdentical(rbind(metadata, newmeta), getMetadata(obj))

  # retrieve named columns
  checkIdentical(snpID, getSnpID(obj))
  checkIdentical(chrom, getChromosome(obj))
  checkIdentical(c(rep("1",5),"X","XY","Y","M","U"),
                 getChromosome(obj, char=TRUE))
  checkIdentical(pos, getPosition(obj))
  checkIdentical(alleleA, getAlleleA(obj))
  checkIdentical(alleleB, getAlleleB(obj))
  
  # other columns
  checkTrue(hasVariable(obj, "rsID"))
  checkIdentical(rsID, getVariable(obj,"rsID"))
  checkTrue(!hasVariable(obj, "foo"))
  checkIdentical(NULL, getVariable(obj, "foo"))
  vars <- c("snpID", "rsID")
  checkIdentical(x[,vars], getVariable(obj, vars))

  # check indexing
  sel <- 1:5 # numeric
  checkIdentical(snpID[sel], getSnpID(obj, index=sel))
  sel <- snpID > 5 # logical
  checkIdentical(snpID[sel], getSnpID(obj, index=sel))

  # check condition
  cond <- "LIMIT 10"
  checkIdentical(snpID[1:10], getSnpID(obj, condition=cond))
  cond <- "WHERE chromosome=1"
  checkIdentical(snpID[chrom == 1], getSnpID(obj, condition=cond))

  close(obj)

  # check alternate chromosome codes
  obj <- SnpAnnotationSQLite(dbpath, YchromCode=24L, XYchromCode=25L)
  checkIdentical(c(rep("1",5),"X","Y","XY","M","U"),
                 getChromosome(obj, char=TRUE))
  
  close(obj)
  file.remove(dbpath)
}
