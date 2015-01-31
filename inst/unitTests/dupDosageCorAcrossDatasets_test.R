
test_dosCorSelectGenotype <- function() {
  # make GenotypeData object
  # using test Illumina dataset
  gds.file <- system.file("extdata", "illumina_geno.gds", package="GWASdata")
  gds <- GdsGenotypeReader(gds.file)
  # select random 20 SNPs
  snp.sel <- sample(getSnpID(gds),20)
  # select random 10 scans
  samp.sel <- sample(getScanID(gds),10)
  close(gds)

  # subset to 20 SNPs and 10 samples
  gds.sub <- tempfile()
  gdsSubset(gds.file, gds.sub,
            snp.include=snp.sel,
            sample.include=samp.sel)

  gds <- GdsGenotypeReader(gds.sub)

  # create snpAnnot
  snpAnnot <- SnpAnnotationDataFrame(data.frame(snpID=getSnpID(gds),
                                                chromosome=getChromosome(gds),
                                                position=getPosition(gds),
                                                alleleA=getAlleleA(gds),
                                                alleleB=getAlleleB(gds)))

  # create scanAnnot
  scanAnnot <- ScanAnnotationDataFrame(data.frame(scanID=getScanID(gds)))

  genoData <- GenotypeData(gds, snpAnnot, scanAnnot)

  ## correct way to source?
  # shuffle scan order - exclude 1 SNP
  scan.order <- getScanID(gds)[c(5:10,1:4)]
  snp.order <- getSnpID(gds)[c(1:8,10:20)]
  #  geno.sel <- GWASTools:::.dosCorSelectGenotype( )
  geno.sel <- GWASTools:::.dosCorSelectGenotype(genoData,
                                    scanIDs=scan.order,
                                    snpIDs=snp.order)

  # now check selection
  rownames(geno.sel) <- NULL # remove snpIDs as row names  
  # loop over scans (columns)
  for (i in getScanID(gds)){
    geno.allsnps <- getGenotype(gds,snp=c(1,-1),scan=c(which(getScanID(gds)==i),1))
    geno.chk <- geno.allsnps[c(1:8,10:20)]
    checkIdentical(geno.sel[,as.character(i)],geno.chk)
  }

  # checkException - check error conditions. question - should these only happen for checks that are explicit in the function itself?

  # non integer scanIDs
  checkException(GWASTools:::.dosCorSelectGenotype(genoData, scanIDs=letters[1:3], snpIDs=snp.order))  
  # checkException(.dosCorSelectGenotype(genoData, scanIDs=letters[1:3], snpIDs=snp.order))

  # non integer snpIDs
  checkException(GWASTools:::.dosCorSelectGenotype(genoData, scanIDs=scan.order, snpIDs=letters[1:3]))
  # checkException(.dosCorSelectGenotype(genoData, scanIDs=scan.order, snpIDs=letters[1:3]))
  
  # unlink tempfile
  unlink(gds.sub)

}


test_dupDosageCorAcrossDatasets <- function() {
  # make first GenotypeData object 
  # using test Illumina dataset
  gds.file <- system.file("extdata", "illumina_geno.gds", package="GWASdata")
  gds <- GdsGenotypeReader(gds.file)
  # select random 20 SNPs
  snp.sel <- sample(getSnpID(gds),20)
  # select random 10 scans
  samp.sel <- sample(getScanID(gds),10)
  close(gds)

  # subset to 20 SNPs and 10 scans
  gds.sub <- tempfile()
  gdsSubset(gds.file, gds.sub,
            snp.include=snp.sel,
            sample.include=samp.sel)
  gds <- GdsGenotypeReader(gds.sub)

  # create snpAnnot
  snpAnnot <- SnpAnnotationDataFrame(data.frame(snpID=getSnpID(gds),
                                                chromosome=getChromosome(gds),
                                                position=getPosition(gds),
                                                alleleA=getAlleleA(gds),
                                                alleleB=getAlleleB(gds),
                                                stringsAsFactors=FALSE))

  # create scanAnnot
  scanAnnot <- ScanAnnotationDataFrame(data.frame(scanID=getScanID(gds),
                                                  subjectID=letters[1:nscan(gds)],
                                                  sex=sample(c("M","F"),size=nscan(gds),replace=TRUE),
                                                  stringsAsFactors=FALSE))          
  
                                       
  genoData1 <- GenotypeData(gds, snpAnnot, scanAnnot)

  # make second GenotypeData as variation on genoData1
  geno1.matx <- getGenotype(gds,snp=c(1,-1), scan=c(1,-1))
  geno2.matx <- abs(geno1.matx-0.1)
  mgr <- MatrixGenotypeReader(genotype=geno2.matx,snpID=getSnpID(gds),
                                         chromosome=getChromosome(gds),
                                         position=getPosition(gds),
                                         scanID=1:nscan(gds))
  # change scanIDs in second dataset
  scan2 <- pData(scanAnnot)
  scan2$scanID <- 1:nrow(scan2)
  # remove one of the overlapping samples
  scan2$subjectID[4] <- "z"
  # change allele A and alleleB mapping for one SNP
  snp1 <- pData(snpAnnot)
  snp2 <- snp1
  snp2$alleleA[3] <- snp1$alleleB[3]
  snp2$alleleB[3] <- snp1$alleleA[3]
  
  genoData2 <- GenotypeData(mgr, SnpAnnotationDataFrame(snp2),
                            scanAnnot=ScanAnnotationDataFrame(scan2))

  ## test genoData1 vs genoData1 --
  out.equal <- dupDosageCorAcrossDatasets(genoData1=genoData1,
                                          genoData2=genoData1,
                                          snp.block.size=19,
                                          scan.block.size=9)
  snp.equal <- out.equal$snp
  samp.equal <- out.equal$samps

  # dosage r2 should be 1 or NA
  checkEquals(snp.equal$dosageR2[!is.na(snp.equal$dosageR2)],
              rep(1,sum(!is.na(snp.equal$dosageR2))))
  checkTrue(max(snp.equal$nsamp.dosageR2)<=nrow(samp.equal))

  checkEquals(samp.equal$dosageR2[!is.na(samp.equal$dosageR2)],
              rep(1,sum(!is.na(samp.equal$dosageR2))))
  checkTrue(max(samp.equal$nsnp.dosageR2)<=nrow(snp.equal))  

  # check exception - no matching samples
  scan3 <- pData(scanAnnot)
  scan3$subjectID <- paste0(scan3$subjectID,"x")
  genoData3 <-  GenotypeData(gds, snpAnnot, ScanAnnotationDataFrame(scan3))
  checkIdentical(NULL, dupDosageCorAcrossDatasets(genoData1=genoData1, genoData2=genoData3))
  # or check warnings - convert warnings to errors
  # options(warn=2)
  # checkException(dupDosageCorAcrossDatasets(genoData1=genoData1, genoData2=genoData3))
  # options(warn=0)  
  

  # check > 1 dup pair per sample - should still work, but select 1 dup pair per sample
  scan4 <- pData(scanAnnot)
  scan4$subjectID[5] <-  scan4$subjectID[4]
  genoData4 <-  GenotypeData(gds, snpAnnot, ScanAnnotationDataFrame(scan4))

  out.twopair <- dupDosageCorAcrossDatasets(genoData1=genoData1,
                                            genoData2=genoData4)
  snp.twopair <- out.twopair$snp
  samp.twopair <- out.twopair$samps
  npair <- length(intersect(unique(scan4$subjectID), getVariable(scanAnnot,"subjectID")))
  checkEquals(npair, nrow(samp.twopair))

  # test genoData1 vs genoData2
  out.diff <- dupDosageCorAcrossDatasets(genoData1=genoData1,
                                        genoData2=genoData2)
  snp.diff <- out.diff$snp
  samp.diff <- out.diff$samp

  # check there are 20 matching snps
  checkEquals(nrow(snp.diff),20)

  # check there are 9 sample dup pairs
  checkEquals(nrow(samp.diff), 9)  
  
  # check output r2 values to manual calculations
  colnames(geno1.matx) <- getScanVariable(genoData1, "subjectID")
  nms <- getSnpVariable(genoData1, c("chromosome","position"))
  rownames(geno1.matx) <- paste(nms$chromosome,nms$position)
  colnames(geno2.matx) <- getScanVariable(genoData2, "subjectID")
  nms <- getSnpVariable(genoData2, c("chromosome","position"))
  rownames(geno2.matx) <- paste(nms$chromosome,nms$position)

  snp.diff$rownames <- paste(snp.diff$chromosome, snp.diff$position)

  # order the matrics to match the function output
  geno1.tmp <- geno1.matx[,is.element(colnames(geno1.matx),samp.diff$subjectID)]
  geno2.tmp <- geno2.matx[,is.element(colnames(geno2.matx),samp.diff$subjectID)]  
  geno1.srt <- geno1.tmp[match(rownames(geno1.tmp),snp.diff$rownames),
                          match(colnames(geno1.tmp),samp.diff$subjectID)]

  geno2.srt <- geno2.tmp[match(rownames(geno2.tmp),snp.diff$rownames),
                          match(colnames(geno2.tmp),samp.diff$subjectID)]

  stopifnot(allequal(rownames(geno1.srt),rownames(geno2.srt)))
  stopifnot(allequal(rownames(geno1.srt),snp.diff$rownames))
  stopifnot(allequal(colnames(geno1.srt),colnames(geno2.srt)))
  stopifnot(allequal(colnames(geno1.srt),samp.diff$subjectID))

  # adjust for SNPs with flipped A/B
  swapAB <- snp.diff$alleleA1!=snp.diff$alleleA2
  geno2.srt[swapAB,] <- 2-geno2.srt[swapAB,]

  # loop over SNPs
  for (i in 1:nrow(snp.diff)){
    r <- cor(geno1.srt[i,], geno2.srt[i,],use="pairwise.complete.obs")
    r2 <- r*r
    checkEquals(r2, snp.diff$dosageR2[i])
  }

  # loop over samples
  for (i in 1:nrow(samp.diff)){
    r <- cor(geno1.srt[,i], geno2.srt[,i],use="pairwise.complete.obs")
    r2 <- r*r
    checkEquals(r2, samp.diff$dosageR2[i])
  }

  # unlink tempfile
  close(gds)
  unlink(gds.sub)

}
