test_getPlinkGenotype <- function() {
  g1 <- c(0,1,2,NA,0,1,2,NA)
  snpID <- 1:8
  chromosome <- 1:8
  position <- 1:8
  alleleA <- c("A", "A", "A", "A", "G", "G", "G", "G")
  alleleB <- c("T", "T", "T", "T", "C", "C", "C", "C")
  mgr <- MatrixGenotypeReader(matrix(g1, ncol=1), snpID, chromosome, position, 1L)
  snpdf <- SnpAnnotationDataFrame(data.frame(snpID, chromosome, position, alleleA, alleleB,
                                             stringsAsFactors=FALSE))
  gd <- GenotypeData(mgr, snpdf)
  exp <- c("T T", "A T", "A A", "0 0", "C C", "C G", "G G", "0 0")
  pg <- GWASTools:::.getPlinkGenotype(gd, 1, 1,)
  checkIdentical(exp, pg)

  g2 <- rep(NA, 8)
  mgr <- MatrixGenotypeReader(matrix(c(g1, g2), ncol=2), snpID, chromosome, position, 1:2)
  gd <- GenotypeData(mgr, snpdf)
  exp2 <- matrix(c(exp, rep("0 0", 8)), ncol=2)
  pg <- GWASTools:::.getPlinkGenotype(gd, 1, 2)
  checkIdentical(exp2, pg)

  # what happens if we have some missing values in the alleles?
  alleleA <- c("A", "A", "A", "A", NA, "G", "G", NA)
  alleleB <- c("T", "T", NA, "T", "C", "C", "C", NA)
  snpdf <- SnpAnnotationDataFrame(data.frame(snpID, chromosome, position, alleleA, alleleB,
                                             stringsAsFactors=FALSE))
  gd <- GenotypeData(mgr, snpdf)
  pg <- GWASTools:::.getPlinkGenotype(gd, 1, 1)
  checkIdentical(exp, pg)
}

test_plinkWrite <- function() {
  A <- "A"; C <- "C"; G <- "G"; T <- "T"
  scandf <- ScanAnnotationDataFrame(
              data.frame(family=rep(1,3), scanID=21:23, father=rep(0,3),
                         mother=rep(0,3), sex=c("M","F","M"),
                         stringsAsFactors=FALSE))
  snpdf <- SnpAnnotationDataFrame(
             data.frame(snpID=1:7, chromosome=21:27, position=rep(10L,7),
                        rsID=paste("rs", 1:7, sep=""),
                        alleleA=c(A,A,A,A,C,0,C), alleleB=c(C,C,T,G,T,0,0),
                        stringsAsFactors=FALSE))
  geno <- matrix(c(1,2,1,
                   2,1,0,
                   NA,1,2,
                   0,0,2,
                   0,1,0,
                   NA,NA,NA,
                   2,2,2), ncol=3, byrow=TRUE)
  ped <- matrix(c(1, 21, 0, 0, 1, -9, A, C, A, A, 0, 0, G, G, T, T, 0, 0, C, C,
                  1, 22, 0, 0, 2, -9, A, A, A, C, A, T, G, G, C, T, 0, 0, C, C,
                  1, 23, 0, 0, 1, -9, A, C, C, C, A, A, A, A, T, T, 0, 0, C, C),
                nrow=3, byrow=TRUE)
  map <- data.frame(chromosome=c(21,22,23,25,24,26,0),
                    rsID=paste("rs", 1:7, sep=""),
                    mapdist=rep(0, 7), position=rep(10, 7),
                    stringsAsFactors=FALSE)

  mgr <- MatrixGenotypeReader(geno, snpdf$snpID, snpdf$chromosome, snpdf$position,
                              scandf$scanID)
  genoData <- GenotypeData(mgr, scanAnnot=scandf, snpAnnot=snpdf)
  pedfile <- tempfile()
  plinkWrite(genoData, pedfile)
  ped.in <- as.matrix(read.table(paste(pedfile, "ped", sep="."),
                       colClasses=rep("character", ncol(ped)),
                       as.is=TRUE, header=FALSE))
  checkTrue(allequal(ped, ped.in))
  map.in <- read.table(paste(pedfile, "map", sep="."), as.is=TRUE, header=FALSE)
  checkTrue(allequal(map, map.in))

  ## test with all but one sample excluded
  plinkWrite(genoData, pedfile, scan.exclude=22:23)
  ped.in <- as.matrix(read.table(paste(pedfile, "ped", sep="."),
                       colClasses=rep("character", ncol(ped)),
                       as.is=TRUE, header=FALSE))
  checkTrue(allequal(ped[1,,drop=FALSE], ped.in))
  
  unlink(paste(pedfile, "*", sep=""))
}

test_plinkCheck <- function() {
  A <- "A"; C <- "C"; G <- "G"; T <- "T"
  scandf <- ScanAnnotationDataFrame(
              data.frame(family=rep(1,3), scanID=21:23, father=rep(0,3),
                         mother=rep(0,3), sex=c("M","F","M"),
                         stringsAsFactors=FALSE))
  snpdf <- SnpAnnotationDataFrame(
             data.frame(snpID=1:7, chromosome=21:27, position=rep(10L,7),
                        rsID=paste("rs", 1:7, sep=""),
                        alleleA=c(A,A,A,A,C,0,C), alleleB=c(C,C,T,G,T,0,0),
                        stringsAsFactors=FALSE))
  geno <- matrix(c(1,2,1,
                   2,1,0,
                   NA,1,2,
                   0,0,2,
                   0,1,0,
                   NA,NA,NA,
                   2,2,2), ncol=3, byrow=TRUE)
  ped <- matrix(c(1, 21, 0, 0, 1, -9, C, C, A, C, A, A, 0, 0, T, T, G, G, 0, 0,
                  1, 22, 0, 0, 2, -9, C, C, A, A, A, C, A, T, C, T, G, G, 0, 0,
                  1, 23, 0, 0, 1, -9, C, C, A, C, C, C, A, A, T, T, A, A, 0, 0),
                nrow=3, byrow=TRUE)
  map <- data.frame(chromosome=c(0,21,22,23,24,25,26),
                    rsID=paste("rs", c(7,1,2,3,5,4,6), sep=""),
                    mapdist=rep(0, 7), position=rep(10, 7),
                    stringsAsFactors=FALSE)

  mgr <- MatrixGenotypeReader(geno, snpdf$snpID, snpdf$chromosome, snpdf$position,
                              scandf$scanID)
  genoData <- GenotypeData(mgr, scanAnnot=scandf, snpAnnot=snpdf)
  pedfile <- tempfile()
  write.table(ped, file=paste(pedfile, "ped", sep="."), quote=FALSE,
              row.names=FALSE, col.names=FALSE)
  write.table(map, file=paste(pedfile, "map", sep="."), quote=FALSE,
              row.names=FALSE, col.names=FALSE)
  
  log <- tempfile()
  checkTrue(plinkCheck(genoData, pedfile, log))

  # change ped file
  ped2 <- ped
  ped2[1,1] <- 3
  write.table(ped2, file=paste(pedfile, "ped", sep="."), quote=FALSE,
              row.names=FALSE, col.names=FALSE)
  checkTrue(!plinkCheck(genoData, pedfile, log))
  ped2 <- ped
  ped2[2,2] <- 3
  write.table(ped2, file=paste(pedfile, "ped", sep="."), quote=FALSE,
              row.names=FALSE, col.names=FALSE)
  checkTrue(!plinkCheck(genoData, pedfile, log))
  # change map file
  write.table(ped, file=paste(pedfile, "ped", sep="."), quote=FALSE,
              row.names=FALSE, col.names=FALSE)
  map2 <- map
  map2$rsID[1] <- "foo"
  write.table(map2, file=paste(pedfile, "map", sep="."), quote=FALSE,
              row.names=FALSE, col.names=FALSE)
  checkTrue(!plinkCheck(genoData, pedfile, log))
  
  unlink(paste(pedfile, "*", sep=""))
  unlink(log)
}

test_plinkCheck2 <- function() {
  ncfile <- tempfile()
  simulateGenotypeMatrix(n.snps=10, n.chromosomes=26,
                         n.samples=20, filename=ncfile, file.type="ncdf")
  nc.orig <- NcdfGenotypeReader(ncfile)
  
  scanID <- getScanID(nc.orig)
  sex <- c(rep("M", 10), rep("F", 10))
  family <- rep(0, 20)
  father <- rep(0, 20)
  mother <- rep(0, 20)
  scandf <- data.frame(family, scanID, father, mother, sex, stringsAsFactors=FALSE)
  snpID <- getSnpID(nc.orig)
  chromosome <- getChromosome(nc.orig)
  position <- getPosition(nc.orig)
  rsID <- paste("rs", snpID, sep="")
  allele.A <- sample(c("A","T"), nsnp(nc.orig), replace=TRUE)
  allele.B <- sample(c("C","G"), nsnp(nc.orig), replace=TRUE)
  snpdf <- data.frame(snpID, chromosome, position, rsID, allele.A,
                      allele.B, stringsAsFactors=FALSE)
  geno <- getGenotype(nc.orig)
  genoData <- GenotypeData(nc.orig, scanAnnot=ScanAnnotationDataFrame(scandf),
                           snpAnnot=SnpAnnotationDataFrame(snpdf))
  
  pedfile <- tempfile()
  plinkWrite(genoData, pedfile)
  logfile <- tempfile()
  checkTrue(plinkCheck(genoData, pedfile, logfile))
  close(genoData)

  # change ncdf file
  nc.new <- nc_open(ncfile, write=TRUE)
  geno <- ncvar_get(nc.new, "genotype")
  geno[geno == 1] <- 2
  ncvar_put(nc.new, "genotype", geno)
  nc_close(nc.new)
  nc.new <- NcdfGenotypeReader(ncfile)
  genoData <- GenotypeData(nc.new, scanAnnot=ScanAnnotationDataFrame(scandf),
                           snpAnnot=SnpAnnotationDataFrame(snpdf))
  checkTrue(!plinkCheck(genoData, pedfile, logfile))
  close(genoData)
  
  unlink(c(ncfile, logfile, paste(pedfile, "*", sep=".")))
}

test_plinkCheck_map <- function() {
  ncfile <- tempfile()
  simulateGenotypeMatrix(n.snps=10, n.chromosomes=26,
                         n.samples=20, filename=ncfile, file.type="ncdf")
  nc.orig <- NcdfGenotypeReader(ncfile)
  
  scanID <- getScanID(nc.orig)
  sex <- c(rep("M", 10), rep("F", 10))
  family <- rep(0, 20)
  father <- rep(0, 20)
  mother <- rep(0, 20)
  scandf <- data.frame(family, scanID, father, mother, sex, stringsAsFactors=FALSE)
  snpID <- getSnpID(nc.orig)
  chromosome <- getChromosome(nc.orig)
  position <- getPosition(nc.orig)
  rsID <- paste("rs", snpID, sep="")
  allele.A <- sample(c("A","T"), nsnp(nc.orig), replace=TRUE)
  allele.B <- sample(c("C","G"), nsnp(nc.orig), replace=TRUE)
  snpdf <- data.frame(snpID, chromosome, position, rsID, allele.A,
                      allele.B, stringsAsFactors=FALSE)
  geno <- getGenotype(nc.orig)
  genoData <- GenotypeData(nc.orig, scanAnnot=ScanAnnotationDataFrame(scandf),
                           snpAnnot=SnpAnnotationDataFrame(snpdf))
  
  pedfile <- tempfile()
  plinkWrite(genoData, pedfile)
  close(genoData)

  # change chromosomes in netCDF file
  nc.new <- nc_open(ncfile, write=TRUE)
  chr.new <- snpdf$chromosome
  chr.new[1:5] <- 28L
  ncvar_put(nc.new, "chromosome", chr.new)
  nc_close(nc.new)
  snpdf$chromosome <- chr.new
  nc.new <- NcdfGenotypeReader(ncfile)
  genoData <- GenotypeData(nc.new, scanAnnot=ScanAnnotationDataFrame(scandf),
                           snpAnnot=SnpAnnotationDataFrame(snpdf))
  
  # map to describe how to map to PLINK
  map <- snpdf[,c("snpID", "chromosome","rsID","position")]
  # reorder map to test matching on snpID
  map <- map[order(map$chromosome),]
  # specify PLINK conversion
  chr.new <- map$chromosome
  # we know that chr=28 was chr=1 in PLINK
  map$chromosome[chr.new == 28] <- 1
  # need to specify Y/XY swap since this will not happen automatically
  map$chromosome[chr.new == 24] <- 25
  map$chromosome[chr.new == 25] <- 24
  
  logfile <- tempfile()
  checkTrue(plinkCheck(genoData, pedfile, logfile, map.alt=map))
  close(genoData)
  
  unlink(c(ncfile, logfile, paste(pedfile, "*", sep=".")))
}

