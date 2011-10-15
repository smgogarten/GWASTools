alleleFrequency <- function(genoData, # object of type GenotypeData
         scan.exclude = NULL, # vector of scanIDs to exclude
         verbose = TRUE) {

  # This function gets the frequency of the A allele over specified scans

  # check that sex in included in genoData
  if (!hasSex(genoData)) stop("sex not found in genoData")

  # get vectors for included scans
  scanID <- getScanID(genoData)
  sex <- getSex(genoData)
  males <- which(!is.na(sex) & sex == "M" & !is.element(scanID, scan.exclude))
  females <- which(!is.na(sex) & sex == "F" & !is.element(scanID, scan.exclude))

  # get chromosome types
  chrom <- getChromosome(genoData, char=TRUE)
  xchr <- chrom == "X"
  ychr <- chrom == "Y"
  
  # Object to hold allelic frequency results
  snpID <- getSnpID(genoData)
  afreq <- matrix(0, length(snpID), 3, dimnames=list(snpID, c("M","F","all")))
  nonmiss <- matrix(0, length(snpID), 3, dimnames=list(snpID, c("M","F","all")))
	
  N <- length(scanID)
  for(i in 1:N){    
    if (verbose & (i %% 100 == 0)) {
      message(paste("reading scan", i, "of", N))
    }
    if (i %in% males) {
      geno <- getGenotype(genoData, snp=c(1,-1), scan=c(i,1))

      # all but sex chroms
      noxy <- !(xchr | ychr | is.na(geno))
      afreq[noxy,"M"] <- afreq[noxy,"M"] + geno[noxy]
      nonmiss[noxy,"M"] <- nonmiss[noxy,"M"] + 2 # 2 alleles per genotype

      # X and Y - don't count heterozygotes
      xy <- (xchr | ychr) & !is.na(geno) & (geno != 1)
      afreq[xy,"M"] <- afreq[xy,"M"] + (geno[xy]/2) # 1 allele coded as 2
      nonmiss[xy,"M"] <- nonmiss[xy,"M"] + 1

    } else if (i %in% females) {
      geno <- getGenotype(genoData, snp=c(1,-1), scan=c(i,1))

      # all but Y
      noy <- !(ychr | is.na(geno))
      afreq[noy,"F"] <- afreq[noy,"F"] + geno[noy]
      nonmiss[noy,"F"] <- nonmiss[noy,"F"] + 2 # 2 alleles per genotype
    }
  }
  afreq[,"all"] <- rowSums(afreq[,c("M","F")])
  nonmiss[,"all"] <- rowSums(nonmiss[,c("M","F")])
  afreq.frac <- afreq/nonmiss
  return(afreq.frac)
}
