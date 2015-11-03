setMissingGenotypes <- function(
	parent.file,	# name of the parent GDS file
	new.file,		# name of the new GDS file to create
        regions, # data frame with regions
        file.type = c("gds", "ncdf"),
	sample.include=NULL,	# vector of sampleIDs for samples to include in new.file
        compress = "ZIP_RA",
        verbose=TRUE) {

  stopifnot(all(c("scanID", "chromosome", "left.base", "right.base", "whole.chrom") %in% names(regions)))

  ## get file type
  file.type <- match.arg(file.type)

  if (file.type == "gds") {
      old <- GdsGenotypeReader(parent.file)
  } else if (file.type == "ncdf") {
      old <- NcdfGenotypeReader(parent.file)
  }

  snpID <- getSnpID(old)
  chrom <- getChromosome(old)
  pos <- getPosition(old)
  scanID <- getScanID(old)
  snp.annotation <- data.frame(snpID, chromosome=chrom, position=pos)
  
  if (is.null(sample.include)) {
    sample.include <- scanID
  } else {
    stopifnot(all(sample.include %in% scanID))
  }

  ## create data file
  if (file.type == "gds") {
      if (hasVariable(old, "snp.rs.id")) {
          snp.annotation$snpName <- getVariable(old, "snp.rs.id")
      }
      if (hasVariable(old, "snp.allele")) {
          snp.annotation$alleleA <- getAlleleA(old)
          snp.annotation$alleleB <- getAlleleB(old)
      }
      gfile <- .createGds(snp.annotation, new.file, "genotype", compress=compress)
      put.attr.gdsn(index.gdsn(gfile, "snp.chromosome"), "autosome.start", min(autosomeCode(old)))
      put.attr.gdsn(index.gdsn(gfile, "snp.chromosome"), "autosome.end", max(autosomeCode(old)))
      put.attr.gdsn(index.gdsn(gfile, "snp.chromosome"), "X", XchromCode(old))
      put.attr.gdsn(index.gdsn(gfile, "snp.chromosome"), "XY", XYchromCode(old))
      put.attr.gdsn(index.gdsn(gfile, "snp.chromosome"), "Y", YchromCode(old))
      put.attr.gdsn(index.gdsn(gfile, "snp.chromosome"), "M", MchromCode(old))
      put.attr.gdsn(index.gdsn(gfile, "snp.chromosome"), "MT", MchromCode(old))
  } else if (file.type == "ncdf") {
      gfile <- .createNcdf(snp.annotation, new.file, "genotype", length(sample.include))
      att.put.ncdf(gfile, 0, "array_name", getAttribute(old, "array_name")) 
      att.put.ncdf(gfile, 0, "genome_build", getAttribute(old, "genome_build"))
  }

  ind.old <- which(scanID %in% sample.include)
  for (ind.new in 1:length(sample.include)) {
    if (verbose & ind.new%%10==0) message(paste("sample", ind.new, "of", length(sample.include)))
    geno <- getGenotype(old, snp=c(1,-1), scan=c(ind.old[ind.new],1))
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
    .addData(gfile, "genotype", list("genotype"=geno), sample.include[ind.new], ind.new)
  }

  .close(gfile, verbose=verbose)
  close(old)
  return(invisible(NULL))
}
