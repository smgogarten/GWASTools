ncdfSetMissingGenotypes <- function(
	parent.file,	# name of the parent netCDF file
	new.file,		# name of the new netCDF file to create
        regions, # data frame with regions
	sample.include=NULL,	# vector of sampleIDs for samples to include in new.file
        verbose=TRUE
) {

.Deprecated("setMissingGenotypes")
    
  stopifnot(all(c("scanID", "chromosome", "left.base", "right.base", "whole.chrom") %in% names(regions)))

  nc.old <- NcdfGenotypeReader(parent.file)
  snpID <- getSnpID(nc.old)
  chrom <- getChromosome(nc.old)
  pos <- getPosition(nc.old)
  scanID <- getScanID(nc.old)
  if (is.null(sample.include)) {
    sample.include <- scanID
  } else {
    stopifnot(all(sample.include %in% scanID))
  }
  
  snpdim <- dim.def.ncdf("snp", "count", snpID)
  sampledim <- dim.def.ncdf("sample", "count", 1:length(sample.include), unlim=TRUE)
  varID <- var.def.ncdf("sampleID", "id", dim=sampledim, missval=0, prec="integer")
  varpos <- var.def.ncdf("position", "bases", dim=snpdim, missval=-1, prec="integer")
  varchr <- var.def.ncdf("chromosome", "id", dim=snpdim, missval=-1, prec="integer")
  vargeno <- var.def.ncdf("genotype", "num_A_alleles", dim=list(snpdim,sampledim), missval=-1, prec="byte")
  nc.new <- create.ncdf(new.file, list(varID, varpos, varchr, vargeno))
  put.var.ncdf(nc.new, varID, sample.include)
  put.var.ncdf(nc.new, varpos, pos)
  put.var.ncdf(nc.new, varchr, chrom)
  att.put.ncdf(nc.new, 0, "array_name", getAttribute(nc.old, "array_name")) 
  att.put.ncdf(nc.new, 0, "genome_build", getAttribute(nc.old, "genome_build"))

  ind.old <- which(scanID %in% sample.include)
  for (ind.new in 1:length(sample.include)) {
    if (verbose & ind.new%%10==0) message(paste("sample", ind.new, "of", length(sample.include)))
    geno <- getGenotype(nc.old, snp=c(1,-1), scan=c(ind.old[ind.new],1))
    reg <- regions[regions$scanID %in% scanID[ind.old[ind.new]],]
    if (nrow(reg) > 0) {
      for (a in 1:nrow(reg)) {
        if (reg$whole.chrom[a]) {
          geno[chrom == reg$chromosome[a]] <- NA
        } else {
          geno[chrom == reg$chromosome[a] & reg$left.base[a] <= pos & pos <= reg$right.base[a]] <- NA
        }
      }
    }
    put.var.ncdf(nc.new, vargeno, geno, start=c(1,ind.new), count=c(-1,1))
  }
  
  close.ncdf(nc.new)
  close(nc.old)
  return(invisible(NULL))
}
