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
                         n.samples=20, ncdf.filename=ncfile)
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
  nc.new <- open.ncdf(ncfile, write=TRUE)
  geno <- get.var.ncdf(nc.new, "genotype")
  geno[geno == 1] <- 2
  put.var.ncdf(nc.new, "genotype", geno)
  close.ncdf(nc.new)
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
                         n.samples=20, ncdf.filename=ncfile)
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
  nc.new <- open.ncdf(ncfile, write=TRUE)
  chr.new <- snpdf$chromosome
  chr.new[1:5] <- 28L
  put.var.ncdf(nc.new, "chromosome", chr.new)
  close.ncdf(nc.new)
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


test_plinkToNcdf <- function() {
  A <- "A"; C <- "C"; G <- "G"; T <- "T"
  ped <- matrix(c(1, 21, 0, 0, 1, -9, C, A, A, A, 0, 0, G, G, T, T, 0, 0, C, C,
                  1, 22, 0, 0, 2, -9, A, A, A, C, A, T, G, G, C, T, 0, 0, C, C,
                  1, 23, 0, 0, 1, -9, A, C, C, C, A, A, A, A, T, T, 0, 0, C, C),
                nrow=3, byrow=TRUE)
  map <- data.frame(chromosome=1:7, rsID=paste("rs", 1:7, sep=""),
                    mapdist=rep(0, 7), position=rep(10, 7),
                    stringsAsFactors=FALSE)
  bim <- cbind(map, alleleA=c(A,A,A,A,C,0,C), alleleB=c(C,C,T,G,T,0,0),
               stringsAsFactors=FALSE)
  # what we expect to see if using alleles provided
  geno <- matrix(c(1,2,1,
                   2,1,0,
                   NA,1,2,
                   0,0,2,
                   0,1,0,
                   NA,NA,NA,
                   2,2,2), ncol=3, byrow=TRUE)
  # what we expect to see if deducing alleles from data
  alleleA.map <- c(C,A,A,G,T,NA,C)
  alleleB.map <- c(A,C,T,A,C,NA,NA)
  geno.map <- matrix(c(1,0,1,
                       2,1,0,
                       NA,1,2,
                       2,2,0,
                       2,1,2,
                       NA,NA,NA,
                       2,2,2), ncol=3, byrow=TRUE)

  # input files
  prefix <- tempfile()
  pedfile <- paste(prefix, "ped", sep=".")
  write.table(ped, file=pedfile, quote=FALSE, row.names=FALSE, col.names=FALSE)
  mapfile <- paste(prefix, "map", sep=".")
  write.table(map, file=mapfile, quote=FALSE, row.names=FALSE, col.names=FALSE)
  bimfile <- paste(prefix, "bim", sep=".")
  write.table(bim, file=bimfile, quote=FALSE, row.names=FALSE, col.names=FALSE)

  # output files
  ncfile <- tempfile()
  scanfile <- tempfile()
  snpfile <- tempfile()
  log <- tempfile()
  
  # test with map file
  plinkToNcdf(pedfile, mapfile, nSamples=nrow(ped), ncdfFile=ncfile,
              snpAnnotFile=snpfile, scanAnnotFile=scanfile)
  snpAnnot <- getobj(snpfile)
  checkEquals(alleleA.map, snpAnnot$alleleA)
  checkEquals(alleleB.map, snpAnnot$alleleB)
  scanAnnot <- getobj(scanfile)
  # individ should be NULL if it was unique integer (now called scanID)
  checkIdentical(NULL, scanAnnot$individ)
  nc <- NcdfGenotypeReader(ncfile)
  genoData <- GenotypeData(nc, scanAnnot=scanAnnot, snpAnnot=snpAnnot)
  checkEquals(geno.map, getGenotype(genoData))
  checkTrue(plinkCheck(genoData, prefix, log))
  close(genoData)
  
  # test with bim file
  plinkToNcdf(pedfile, bimfile, nSamples=nrow(ped), ncdfFile=ncfile,
              snpAnnotFile=snpfile, scanAnnotFile=scanfile)
  snpAnnot <- getobj(snpfile)
  checkEquals(bim$alleleA[bim$alleleA != 0], snpAnnot$alleleA[!is.na(snpAnnot$alleleA)])
  checkEquals(bim$alleleB[bim$alleleB != 0], snpAnnot$alleleB[!is.na(snpAnnot$alleleB)])
  checkEquals(bim$alleleA == 0, is.na(snpAnnot$alleleA))
  checkEquals(bim$alleleB == 0, is.na(snpAnnot$alleleB))
  scanAnnot <- getobj(scanfile)
  nc <- NcdfGenotypeReader(ncfile)
  genoData <- GenotypeData(nc, scanAnnot=scanAnnot, snpAnnot=snpAnnot)
  checkEquals(geno, getGenotype(genoData))
  checkTrue(plinkCheck(genoData, prefix, log))
  close(genoData)

  # test with non-unique individual ID
  ped[,1] <- 1:3
  ped[,2] <- rep(1,3)
  write.table(ped, file=pedfile, quote=FALSE, row.names=FALSE, col.names=FALSE)
  plinkToNcdf(pedfile, mapfile, nSamples=nrow(ped), ncdfFile=ncfile,
              snpAnnotFile=snpfile, scanAnnotFile=scanfile)
  snpAnnot <- getobj(snpfile)
  scanAnnot <- getobj(scanfile)
  checkEquals(scanAnnot$scanID, 1:3)
  checkEquals(scanAnnot$individ, ped[,2])
  nc <- NcdfGenotypeReader(ncfile)
  genoData <- GenotypeData(nc, scanAnnot=scanAnnot, snpAnnot=snpAnnot)
  close(genoData)
  
  # test with character individual ID
  ped[,1] <- rep(1,3)
  ped[,2] <- c("A1", "B1", "C1")
  write.table(ped, file=pedfile, quote=FALSE, row.names=FALSE, col.names=FALSE)
  plinkToNcdf(pedfile, mapfile, nSamples=nrow(ped), ncdfFile=ncfile,
              snpAnnotFile=snpfile, scanAnnotFile=scanfile)
  snpAnnot <- getobj(snpfile)
  scanAnnot <- getobj(scanfile)
  checkEquals(scanAnnot$scanID, 1:3)
  nc <- NcdfGenotypeReader(ncfile)
  genoData <- GenotypeData(nc, scanAnnot=scanAnnot, snpAnnot=snpAnnot)
  checkTrue(plinkCheck(genoData, prefix, log, individual.col="individ"))
  close(genoData)
  
  # test with character chromosomes
  map$chromosome <- c("X", "Y", "XY", "MT", "0", "X", "Y")
  write.table(map, file=mapfile, quote=FALSE, row.names=FALSE, col.names=FALSE)
  plinkToNcdf(pedfile, mapfile, nSamples=nrow(ped), ncdfFile=ncfile,
              snpAnnotFile=snpfile, scanAnnotFile=scanfile)
  snpAnnot <- getobj(snpfile)
  scanAnnot <- getobj(scanfile)
  nc <- NcdfGenotypeReader(ncfile)
  genoData <- GenotypeData(nc, scanAnnot=scanAnnot, snpAnnot=snpAnnot)
  checkTrue(plinkCheck(genoData, prefix, log, individual.col="individ"))
  checkEquals(c(23,23,24,25,25,26,27), getChromosome(genoData))
  checkEquals(c("X", "X", "XY", "Y", "Y", "M", "U"), getChromosome(genoData, char=TRUE))
  close(genoData)

  # bimfile
  bim$chromosome <- c("X", "Y", "XY", "MT", "0", "X", "Y")
  write.table(bim, file=bimfile, quote=FALSE, row.names=FALSE, col.names=FALSE)
  plinkToNcdf(pedfile, bimfile, nSamples=nrow(ped), ncdfFile=ncfile,
              snpAnnotFile=snpfile, scanAnnotFile=scanfile)
  snpAnnot <- getobj(snpfile)
  scanAnnot <- getobj(scanfile)
  nc <- NcdfGenotypeReader(ncfile)
  genoData <- GenotypeData(nc, scanAnnot=scanAnnot, snpAnnot=snpAnnot)
  checkTrue(plinkCheck(genoData, prefix, log, individual.col="individ"))
  checkEquals(c(23,23,24,25,25,26,27), getChromosome(genoData))
  checkEquals(c("X", "X", "XY", "Y", "Y", "M", "U"), getChromosome(genoData, char=TRUE))
  close(genoData)

  # PLINK integer codes
  map$chromosome <- c(0, 23, 24, 24, 25, 25, 26)
  write.table(map, file=mapfile, quote=FALSE, row.names=FALSE, col.names=FALSE)
  plinkToNcdf(pedfile, mapfile, nSamples=nrow(ped), ncdfFile=ncfile,
              snpAnnotFile=snpfile, scanAnnotFile=scanfile)
  snpAnnot <- getobj(snpfile)
  scanAnnot <- getobj(scanfile)
  nc <- NcdfGenotypeReader(ncfile)
  genoData <- GenotypeData(nc, scanAnnot=scanAnnot, snpAnnot=snpAnnot)
  checkTrue(plinkCheck(genoData, prefix, log, individual.col="individ"))
  checkEquals(c(23,24,24,25,25,26,27), getChromosome(genoData))
  checkEquals(c("X", "XY", "XY", "Y", "Y", "M", "U"), getChromosome(genoData, char=TRUE))
  close(genoData)
  
  # test with different coding
  map$chromosome <- c("X", "Y", "XY", "MT", "0", "X", "Y")
  write.table(map, file=mapfile, quote=FALSE, row.names=FALSE, col.names=FALSE)
  plinkToNcdf(pedfile, mapfile, nSamples=nrow(ped), ncdfFile=ncfile,
              snpAnnotFile=snpfile, scanAnnotFile=scanfile,
              ncdfXchromCode=15, ncdfXYchromCode=14, ncdfYchromCode=13,
              ncdfMchromCode=12, ncdfUchromCode=0)
  snpAnnot <- getobj(snpfile)
  scanAnnot <- getobj(scanfile)
  nc <- NcdfGenotypeReader(ncfile, XchromCode=15L, XYchromCode=14L,
                           YchromCode=13L, MchromCode=12L)
  genoData <- GenotypeData(nc, scanAnnot=scanAnnot, snpAnnot=snpAnnot)
  checkTrue(plinkCheck(genoData, prefix, log, individual.col="individ"))
  checkEquals(c(0,12,13,13,14,15,15), getChromosome(genoData))
  checkEquals(c("U", "M", "Y", "Y", "XY", "X", "X"), getChromosome(genoData, char=TRUE))
  close(genoData)
  
  unlink(paste(prefix, "*", sep=""))
  unlink(c(ncfile, scanfile, snpfile))
  unlink(log)
}


test_plinkToNcdf2<- function() {
  ncfile <- tempfile()
  simulateGenotypeMatrix(n.snps=1, n.chromosomes=26,
                         n.samples=10, ncdf.filename=ncfile)
  nc.orig <- NcdfGenotypeReader(ncfile)
  
  scanID <- getScanID(nc.orig)
  sex <- c(rep("M", 5), rep("F", 5))
  family <- rep(0, 10)
  father <- rep(0, 10)
  mother <- rep(0, 10)
  scandf <- data.frame(family, scanID, father, mother, sex, stringsAsFactors=FALSE)
  snpID <- getSnpID(nc.orig)
  chromosome <- getChromosome(nc.orig)
  position <- getPosition(nc.orig)
  rsID <- paste("rs", snpID, sep="")
  alleleA <- sample(c("A","T"), nsnp(nc.orig), replace=TRUE)
  alleleB <- sample(c("C","G"), nsnp(nc.orig), replace=TRUE)
  snpdf <- data.frame(snpID, chromosome, position, rsID, alleleA,
                      alleleB, stringsAsFactors=FALSE)
  geno <- getGenotype(nc.orig)
  genoData <- GenotypeData(nc.orig, scanAnnot=ScanAnnotationDataFrame(scandf),
                           snpAnnot=SnpAnnotationDataFrame(snpdf))
  
  pedfile <- tempfile()
  plinkWrite(genoData, pedfile)
  close(genoData)

  # test without alleles
  pedFile <- paste(pedfile, "ped", sep=".")
  mapFile <- paste(pedfile, "map", sep=".")
  ncfile2 <- tempfile()
  scanfile <- tempfile()
  snpfile <- tempfile()
  plinkToNcdf(pedFile, mapFile, nSamples=10, ncdfFile=ncfile2,
              snpAnnotFile=snpfile, scanAnnotFile=scanfile)
  snpAnnot <- getobj(snpfile)
  scanAnnot <- getobj(scanfile)
  nc <- NcdfGenotypeReader(ncfile2)
  genoData2 <- GenotypeData(nc, scanAnnot=scanAnnot, snpAnnot=snpAnnot)
  scancols <- c("scanID", "family", "father", "mother", "sex")
  checkTrue(allequal(pData(scanAnnot)[,scancols], scandf[,scancols]))
  snpcols <- c("snpID", "chromosome", "position", "rsID")
  checkTrue(allequal(pData(snpAnnot)[,snpcols], snpdf[,snpcols]))
  #cannot compare genotypes directly since definition of A and B alleles
  # has changed  
  log <- tempfile()
  checkTrue(plinkCheck(genoData2, pedfile, log))
  close(genoData2)

  # test with alleles
  map <- read.table(mapFile, as.is=TRUE, header=FALSE)
  map <- cbind(map, alleleA, alleleB)
  write.table(map, file=mapFile, quote=FALSE, row.names=FALSE, col.names=FALSE)
  plinkToNcdf(pedFile, mapFile, nSamples=10, ncdfFile=ncfile2,
              snpAnnotFile=snpfile, scanAnnotFile=scanfile)
  snpAnnot <- getobj(snpfile)
  scanAnnot <- getobj(scanfile)
  nc <- NcdfGenotypeReader(ncfile2)
  genoData2 <- GenotypeData(nc, scanAnnot=scanAnnot, snpAnnot=snpAnnot)
  scancols <- c("scanID", "family", "father", "mother", "sex")
  checkTrue(allequal(pData(scanAnnot)[,scancols], scandf[,scancols]))
  snpcols <- c("snpID", "chromosome", "position", "rsID", "alleleA", "alleleB")
  checkTrue(allequal(pData(snpAnnot)[,snpcols], snpdf[,snpcols]))
  checkEquals(getGenotype(genoData2), geno)  
  checkTrue(plinkCheck(genoData2, pedfile, log)) 
  close(genoData2)
  
  unlink(paste(pedfile, "*", sep=""))
  unlink(c(ncfile, ncfile2, scanfile, snpfile))
  unlink(log)
}
