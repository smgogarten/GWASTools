
test_dosCorSelectGenotype <- function() {
  # make GenotypeData object
  # using test Illumina dataset
  gds.file <- system.file("extdata", "illumina_geno.gds", package="GWASdata")
  gds <- GdsGenotypeReader(gds.file)
  # select random 20 SNPs - make sure 2 are Y chr
  snp.selA <- sample(getSnpID(gds),18)
  snp.selY <- sample(getSnpID(gds)[which(getChromosome(gds, char=TRUE) %in% "Y")],2)
  snp.sel <- c(snp.selA, snp.selY)
  table(getChromosome(gds, char=TRUE)[getSnpID(gds) %in% snp.sel])
  # select random 10 scans
  samp.sel <- sample(getScanID(gds),10)
  close(gds)

  # get scan annotation to keep track of M vs F
  scan.file <- system.file("data", "illuminaScanADF.RData", package="GWASdata")
  scanIlm <- pData(getobj(scan.file))
  scanIlm.sel <- scanIlm[scanIlm$scanID %in% samp.sel,]

  # set one M to F to check Y SNP behavior (females may already be set to NA)
  male.sel <- sample(scanIlm.sel$scanID[scanIlm.sel$sex %in% "M"],1)
  scanIlm.sel$sex[scanIlm.sel$scanID %in% male.sel] <- "F"

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
  scanAnnot <- ScanAnnotationDataFrame(scanIlm.sel)

  genoData <- GenotypeData(gds, snpAnnot, scanAnnot)

  # shuffle scan order 
  scan.order <- getScanID(gds)[c(5:10,1:4)]
  # exclude Y chr SNPs
  snp.order <- getSnpID(gds)[!is.element(getChromosome(gds,char=TRUE),"Y")]

  geno.sel <- GWASTools:::.dosCorSelectGenotype(genoData,
                                    scanIDs=scan.order,
                                    snpIDs=snp.order)

  # now check selection
  rownames(geno.sel) <- NULL # remove snpIDs as row names  
  # loop over scans (columns)
  for (i in getScanID(gds)){
    geno.allsnps <- getGenotype(gds,snp=c(1,-1),scan=c(which(getScanID(gds)==i),1))
    geno.chk <- geno.allsnps[!is.element(getChromosome(gds,char=TRUE),"Y")]
    checkIdentical(geno.sel[,as.character(i)],geno.chk)
  }

  # check female genotypes on the Y - select all females on the Y chr SNPs
  scan.F <- getScanID(genoData)[getSex(genoData) %in% "F"]
  snp.Y <- getSnpID(genoData)[getChromosome(genoData, char=TRUE) %in% "Y"]
  geno.femY.sel <- GWASTools:::.dosCorSelectGenotype(genoData,
                                    scanIDs=scan.F,
                                    snpIDs=snp.Y)

  checkEquals(sum(!is.na(geno.femY.sel)), 0)
  
  # checkException - check error conditions. question - should these only happen for checks that are explicit in the function itself?

  # non integer scanIDs
  checkException(GWASTools:::.dosCorSelectGenotype(genoData, scanIDs=letters[1:3], snpIDs=snp.order), silent=TRUE)  

  # non integer snpIDs
  checkException(GWASTools:::.dosCorSelectGenotype(genoData, scanIDs=scan.order, snpIDs=letters[1:3]), silent=TRUE)
  
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

  # get scan annotation to keep track of M vs F
  scan.file <- system.file("data", "illuminaScanADF.RData", package="GWASdata")
  scanIlm <- pData(getobj(scan.file))
  scanIlm.sel <- scanIlm[scanIlm$scanID %in% samp.sel,]
  # table(scanIlm.sel$sex)

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
  # use sex information from Illumina scan annot
  scanAnnot <- ScanAnnotationDataFrame(data.frame(scanID=getScanID(gds),
                                                  subjectID=letters[1:nscan(gds)],
                                                  sex=scanIlm.sel$sex,
                                                  stringsAsFactors=FALSE))          
  
                                       
  genoData1 <- GenotypeData(gds, snpAnnot, scanAnnot)

  # make second GenotypeData as variation on genoData1
  geno1.matx <- getGenotype(gds,snp=c(1,-1), scan=c(1,-1))
  geno2.matx <- abs(geno1.matx-0.1)
  # remove one of the overlapping SNPs by changing the position
  pos.init <- getPosition(gds)
  pos.init[5] <- pos.init[5]+10
  mgr <- MatrixGenotypeReader(genotype=geno2.matx,snpID=getSnpID(gds),
                                         chromosome=getChromosome(gds),
                                         position=as.integer(pos.init),
                                         scanID=1:nscan(gds))
  # change scanIDs in second dataset
  scan1 <- pData(scanAnnot)
  scan2 <- scan1
  scan2$scanID <- 1:nrow(scan2)
  # remove one of the overlapping samples
  scan2$subjectID[4] <- "z"

  snp1 <- pData(snpAnnot)
  snp2 <- snp1
  # change allele A and alleleB mapping for one SNP
  snp2$position[5] <- as.integer(snp2$position[5] + 10)
  snp2$alleleA[3] <- snp1$alleleB[3]
  snp2$alleleB[3] <- snp1$alleleA[3]
  
  genoData2 <- GenotypeData(mgr, SnpAnnotationDataFrame(snp2),
                            scanAnnot=ScanAnnotationDataFrame(scan2))

  ## test genoData1 vs genoData1 --
  message("\nComparing identical datasets\n")
  out.equal <- dupDosageCorAcrossDatasets(genoData1=genoData1,
                                          genoData2=genoData1,
                                          snp.block.size=18,
                                          scan.block.size=9)
  snp.equal <- out.equal$snp
  samp.equal <- out.equal$samps

  # dosage correlation (r) should be 1 or NA
  checkEquals(snp.equal$dosageR[!is.na(snp.equal$dosageR)],
              rep(1,sum(!is.na(snp.equal$dosageR))))
  checkTrue(max(snp.equal$nsamp.dosageR)<=nrow(samp.equal))

  checkEquals(samp.equal$dosageR[!is.na(samp.equal$dosageR)],
              rep(1,sum(!is.na(samp.equal$dosageR))))
  checkTrue(max(samp.equal$nsnp.dosageR)<=nrow(snp.equal))  

  # check exception - no matching samples
  message("\nComparing no matching samples\n")  
  scan3 <- scan1
  scan3$subjectID <- paste0(scan3$subjectID,"x")
  genoData3 <-  GenotypeData(gds, snpAnnot, ScanAnnotationDataFrame(scan3))
  checkIdentical(NULL, dupDosageCorAcrossDatasets(genoData1=genoData1, genoData2=genoData3))
  # or check warnings - convert warnings to errors
  # options(warn=2)
  # checkException(dupDosageCorAcrossDatasets(genoData1=genoData1, genoData2=genoData3))
  # options(warn=0)  

  # check > 1 dup pair per sample - should still work, but select 1 dup pair per sample
  message("\nComparing >1 dup pair\n")
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
  message("\nComparing different datasets\n")
  out.diff <- dupDosageCorAcrossDatasets(genoData1=genoData1,
                                         genoData2=genoData2)
  snp.diff <- out.diff$snp
  samp.diff <- out.diff$samp

  # check there are 19 matching snps
  checkEquals(nrow(snp.diff), 19)

  # check there are 9 sample dup pairs
  checkEquals(nrow(samp.diff), 9)  
  
  # check output r values to manual calculations
  # use .dosCorSelectGenotype function to get genoData1 and genoData2 in the same snp and sample order

  message("\tselecting genotypes from genoData1")
  geno1.srt <- GWASTools:::.dosCorSelectGenotype(genoData1,
                                                 scanIDs=samp.diff$scanID1,
                                                 snpIDs=snp.diff$snpID1)
  message("\tselecting genotypes from genoData2")
  geno2.srt <- GWASTools:::.dosCorSelectGenotype(genoData2,
                                                scanIDs=samp.diff$scanID2,
                                                snpIDs=snp.diff$snpID2)

  ## row and col names are snpID and scanID, respectively
  # replace rownames with chrom + position
  # [reminder: match returns a vector of the positions of (first) matches of its first argument in its second]
  row1 <- with(snp1[match(rownames(geno1.srt),snp1$snpID),],
                    paste(chromosome, position))
  row2 <- with(snp2[match(rownames(geno2.srt),snp2$snpID),],
                    paste(chromosome, position))
  rownames(geno1.srt) <- row1
  rownames(geno2.srt) <- row2
  
  # replace colnames with subjectID
  col1 <- scan1$subjectID[match(colnames(geno1.srt),scan1$scanID)]
  col2 <- scan2$subjectID[match(colnames(geno2.srt),scan2$scanID)]
  colnames(geno1.srt) <- col1
  colnames(geno2.srt) <- col2  
  
  # check that the two matrices are in the same snp and scan order
  snp.diff$rownames <- paste(snp.diff$chromosome, snp.diff$position)
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
    # r2 <- r*r
    checkEquals(r, snp.diff$dosageR[i])
  }

  # loop over samples
  for (i in 1:nrow(samp.diff)){
    r <- cor(geno1.srt[,i], geno2.srt[,i],use="pairwise.complete.obs")
    # r2 <- r*r
    checkEquals(r, samp.diff$dosageR[i])
  }

  # unlink tempfile
  close(gds)
  unlink(gds.sub)

}
