test_gds <- function() {
  # simulated data 
  ncfile <- tempfile()
  simulateGenotypeMatrix(n.snps=100, n.chromosomes=3,
                         n.samples=20, ncdf.filename=ncfile)
  gdsfile <- tempfile()
  convertNcdfGds(ncfile, gdsfile)
  gds <- GdsGenotypeReader(gdsfile)
  sampID <- getScanID(gds)
  sampsel <- c(1,3,6,7,8,9,10,13,16,20)
  geno <- getGenotype(gds)
  close(gds)

  regions <- data.frame("scanID"=c(1,1,7,19), "chromosome"=c(1,2,3,1),
                     "left.base"=c(NA,10,30,70), "right.base"=c(NA,50,40,100),
                     "whole.chrom"=c(TRUE, FALSE, FALSE, FALSE))
  geno.na <- geno
  geno.na[1:100,1] <- NA
  geno.na[110:150,1] <- NA
  geno.na[230:240,7] <- NA
  geno.na[70:100,19] <- NA

  # test with empty data frame
  reg.null <- regions[0,]
  newfile <- tempfile()
  setMissingGenotypes(gdsfile, newfile, reg.null, file.type="gds", sample.include=NULL)
  ncfile2 <- tempfile()
  convertGdsNcdf(newfile, ncfile2)
  ncdfSubsetCheck(ncfile, ncfile2)
  file.remove(c(newfile, ncfile2))
  
  newfile <- tempfile()
  setMissingGenotypes(gdsfile, newfile, reg.null, file.type="gds", sample.include=sampID[sampsel])
  ncfile2 <- tempfile()
  convertGdsNcdf(newfile, ncfile2)
  ncdfSubsetCheck(ncfile, ncfile2, sample.include=sampID[sampsel])
  file.remove(c(newfile, ncfile2))

  # test with regions
  newfile <- tempfile()
  setMissingGenotypes(gdsfile, newfile, regions, file.type="gds", sample.include=NULL)
  gds.new <- GdsGenotypeReader(newfile)
  checkEquals(geno.na, getGenotype(gds.new))
  close(gds.new)
  file.remove(newfile)
  
  newfile <- tempfile()
  setMissingGenotypes(gdsfile, newfile, regions, file.type="gds", sample.include=sampID[sampsel])
  gds.new <- GdsGenotypeReader(newfile)
  checkEquals(geno.na[,sampsel], getGenotype(gds.new))
  close(gds.new)
  file.remove(newfile)
  
  file.remove(gdsfile)
}

test_ncdf <- function() {
  # simulated data 
  ncfile <- tempfile()
  simulateGenotypeMatrix(n.snps=100, n.chromosomes=3,
                         n.samples=20, ncdf.filename=ncfile)
  nc <- NcdfGenotypeReader(ncfile)
  sampID <- getScanID(nc)
  sampsel <- c(1,3,6,7,8,9,10,13,16,20)
  geno <- getGenotype(nc)
  close(nc)

  regions <- data.frame("scanID"=c(1,1,7,19), "chromosome"=c(1,2,3,1),
                     "left.base"=c(NA,10,30,70), "right.base"=c(NA,50,40,100),
                     "whole.chrom"=c(TRUE, FALSE, FALSE, FALSE))
  geno.na <- geno
  geno.na[1:100,1] <- NA
  geno.na[110:150,1] <- NA
  geno.na[230:240,7] <- NA
  geno.na[70:100,19] <- NA

  # test with empty data frame
  reg.null <- regions[0,]
  newfile <- tempfile()
  setMissingGenotypes(ncfile, newfile, reg.null, file.type="ncdf", sample.include=NULL)
  ncdfSubsetCheck(ncfile, newfile)
  file.remove(newfile)
  
  newfile <- tempfile()
  setMissingGenotypes(ncfile, newfile, reg.null, file.type="ncdf", sample.include=sampID[sampsel])
  ncdfSubsetCheck(ncfile, newfile, sample.include=sampID[sampsel])
  file.remove(newfile)

  # test with regions
  newfile <- tempfile()
  setMissingGenotypes(ncfile, newfile, regions, file.type="ncdf", sample.include=NULL)
  nc.new <- NcdfGenotypeReader(newfile)
  checkEquals(geno.na, getGenotype(nc.new))
  close(nc.new)
  file.remove(newfile)
  
  newfile <- tempfile()
  setMissingGenotypes(ncfile, newfile, regions, file.type="ncdf", sample.include=sampID[sampsel])
  nc.new <- NcdfGenotypeReader(newfile)
  checkEquals(geno.na[,sampsel], getGenotype(nc.new))
  close(nc.new)
  file.remove(newfile)
  
  file.remove(ncfile)
}
