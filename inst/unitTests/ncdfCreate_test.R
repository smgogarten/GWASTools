test_ncdfCreate <- function() {
  snpID <- 1:10
  chrom <- c(rep(1L,5), 23:27)
  pos <- c(101:109, NA)
  annot <- data.frame(snpID=snpID, chromosome=chrom, position=pos)
  file <- tempfile()
  vars <- c("X","Y")
  nsamp <- 5
  arrname <- "Array"
  build <- "Build"

  ncdfCreate(annot, file, variables=vars, n.samples=nsamp,
                  precision="single", array.name=arrname,
                  genome.build=build)
  
  nc <- NcdfIntensityReader(file)
  checkEquals(nsamp, nscan(nc))
  checkEquals(length(snpID), nsnp(nc))
  checkIdentical(c("sampleID", "position", "chromosome", vars),
                 getVariableNames(nc))
  checkIdentical(snpID, getSnpID(nc))
  checkIdentical(chrom, getChromosome(nc))
  checkIdentical(pos, getPosition(nc))
  checkIdentical(arrname, getAttribute(nc, "array_name"))
  checkIdentical(build, getAttribute(nc, "genome_build"))
  close(nc)
  
  # errors
  checkException(ncdfCreate(annot, file, variables="foo")) # wrong var
  annot <- data.frame(snp=snpID, c=chrom, p=pos) # wrong names
  checkException(ncdfCreate(annot, file))
  snpID <- paste("rs", 1:10, sep="") # snpID not an integer
  annot <- data.frame(snpID=snpID, chromosome=chrom, position=pos,
                      stringsAsFactors=FALSE)
  checkException(ncdfCreate(annot, file))
  snpID <- rep(1L, 10) # snpID not unique
  annot <- data.frame(snpID=snpID, chromosome=chrom, position=pos)
  checkException(ncdfCreate(annot, file))
  snpID <- 10:1 # snpID not sorted
  annot <- data.frame(snpID=snpID, chromosome=chrom, position=pos)
  checkException(ncdfCreate(annot, file))
  
  unlink(file)
}
