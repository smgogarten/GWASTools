
hetBySnpSex <- function(genoData, # object of type GenotypeData
         scan.exclude = NULL, 
         verbose = TRUE) {
  
  # check that sex in included in genoData
  if (!hasSex(genoData)) stop("sex not found in genoData")

  # select males and females which are not excluded
  scanID <- getScanID(genoData)
  sex <- getSex(genoData)	
  males <- which(!is.na(sex) & sex == "M" & !is.element(scanID, scan.exclude))
  females <- which(!is.na(sex) & sex == "F" & !is.element(scanID, scan.exclude))

  # matrix to hold heterozygote frequency
  snpID <- getSnpID(genoData)
  het <- matrix(0, length(snpID), 2, dimnames=list(snpID, c("M","F")))
  nonmiss <- matrix(0, length(snpID), 2, dimnames=list(snpID, c("M","F")))

  N <- length(scanID)
  for(i in 1:N){
    if (verbose & (i %% 100 == 0)) {
      message(paste("reading scan", i, "of", N))
    }
    if (i %in% males) {
      geno <- getGenotype(genoData, snp=c(1,-1), scan=c(i,1))
      het[,"M"] <- het[,"M"] + (!is.na(geno) & geno == 1)
      nonmiss[,"M"] <- nonmiss[,"M"] + !is.na(geno)
    } else if (i %in% females) {
      geno <- getGenotype(genoData, snp=c(1,-1), scan=c(i,1))
      het[,"F"] <- het[,"F"] + (!is.na(geno) & geno == 1)
      nonmiss[,"F"] <- nonmiss[,"F"] + !is.na(geno)
    }
  }	
  return(het/nonmiss)			
}
