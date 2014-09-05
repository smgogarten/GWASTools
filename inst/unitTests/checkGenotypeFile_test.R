test_ncdfCheckGenotype <- function() {
  # snp annotation
  data(affy_snp_annot)
  snpAnnot <- affy_snp_annot[,c("snpID", "probeID")]
  names(snpAnnot) <- c("snpID", "snpName")

  # scan annotation
  data(affy_scan_annot)
  scanAnnot <- affy_scan_annot[,c("scanID", "genoRunID", "chpFile")]
  names(scanAnnot) <- c("scanID", "scanName", "file")

  # netCDF
  path <- system.file("extdata", "affy_raw_data", package="GWASdata")
  ncfile <- system.file("extdata", "affy_geno.nc", package="GWASdata")
  col.nums <- as.integer(c(2,3)); names(col.nums) <- c("snp", "geno")
  diagfile <- tempfile()
  res <- ncdfCheckGenotype(path, ncfile, snpAnnot, scanAnnot, sep.type="\t",
                         skip.num=1, col.total=6, col.nums=col.nums,
                         scan.name.in.file=-1, check.scan.index=1:3,
                         n.scans.loaded=3, diagnostics.filename=diagfile)
  checkTrue(all(res$chk == 1))
  checkTrue(all(res$snp.order == 1))
  checkTrue(all(res$geno.chk == 1))
  
  # now change a file and see if it detects the difference
  origfile <- system.file("extdata", "affy_raw_data",
                          scanAnnot$file[1], package="GWASdata")
  orig <- read.table(origfile, as.is=TRUE, header=TRUE, sep="\t")
  orig[1,3] <- 0
  tmpfile <- "tmp.txt"
  tmppath <- tempdir()
  tmp <- paste(tmppath, tmpfile, sep="/")
  write.table(orig, tmp, sep="\t", row.names=FALSE, quote=FALSE)
  scanAnnot$file[1] <- tmpfile
  res <- ncdfCheckGenotype(tmppath, ncfile, snpAnnot, scanAnnot, sep.type="\t",
                         skip.num=1, col.total=6, col.nums=col.nums,
                         scan.name.in.file=-1, check.scan.index=1,
                         n.scans.loaded=3, diagnostics.filename=diagfile)
  checkTrue(res$geno.chk[1] == 0)

  write.table(orig[1:100,], tmp, sep="\t", row.names=FALSE, quote=FALSE)
  res <- ncdfCheckGenotype(tmppath, ncfile, snpAnnot, scanAnnot, sep.type="\t",
                         skip.num=1, col.total=6, col.nums=col.nums,
                         scan.name.in.file=-1, check.scan.index=1,
                         n.scans.loaded=3, diagnostics.filename=diagfile)
  checkTrue(res$row.num[1] == 100)

  file.remove(tmp,diagfile)
}


test_checkGenotypeFile_ncdf <- function() {
  # snp annotation
  data(affy_snp_annot)
  snpAnnot <- affy_snp_annot[,c("snpID", "probeID")]
  names(snpAnnot) <- c("snpID", "snpName")

  # scan annotation
  data(affy_scan_annot)
  scanAnnot <- affy_scan_annot[,c("scanID", "genoRunID", "chpFile")]
  names(scanAnnot) <- c("scanID", "scanName", "file")

  # netCDF
  path <- system.file("extdata", "affy_raw_data", package="GWASdata")
  ncfile <- system.file("extdata", "affy_geno.nc", package="GWASdata")
  col.nums <- as.integer(c(2,3)); names(col.nums) <- c("snp", "geno")
  diagfile <- tempfile()
  res <- checkGenotypeFile(path, ncfile, file.type="ncdf", snpAnnot, scanAnnot, sep.type="\t",
                         skip.num=1, col.total=6, col.nums=col.nums,
                         scan.name.in.file=-1, check.scan.index=1:3,
                         n.scans.loaded=3, diagnostics.filename=diagfile)
  checkTrue(all(res$chk == 1))
  checkTrue(all(res$snp.order == 1))
  checkTrue(all(res$geno.chk == 1))
  
  # now change a file and see if it detects the difference
  origfile <- system.file("extdata", "affy_raw_data",
                          scanAnnot$file[1], package="GWASdata")
  orig <- read.table(origfile, as.is=TRUE, header=TRUE, sep="\t")
  orig[1,3] <- 0
  tmpfile <- "tmp.txt"
  tmppath <- tempdir()
  tmp <- paste(tmppath, tmpfile, sep="/")
  write.table(orig, tmp, sep="\t", row.names=FALSE, quote=FALSE)
  scanAnnot$file[1] <- tmpfile
  res <- checkGenotypeFile(tmppath, ncfile, file.type="ncdf", snpAnnot, scanAnnot, sep.type="\t",
                         skip.num=1, col.total=6, col.nums=col.nums,
                         scan.name.in.file=-1, check.scan.index=1,
                         n.scans.loaded=3, diagnostics.filename=diagfile)
  checkTrue(res$geno.chk[1] == 0)

  write.table(orig[1:100,], tmp, sep="\t", row.names=FALSE, quote=FALSE)
  res <- checkGenotypeFile(tmppath, ncfile, file.type="ncdf", snpAnnot, scanAnnot, sep.type="\t",
                         skip.num=1, col.total=6, col.nums=col.nums,
                         scan.name.in.file=-1, check.scan.index=1,
                         n.scans.loaded=3, diagnostics.filename=diagfile)
  checkTrue(res$row.num[1] == 100)

  file.remove(tmp,diagfile)
}

test_checkGenotypeFile_gds <- function() {
  # snp annotation
  data(affy_snp_annot)
  snpAnnot <- affy_snp_annot[,c("snpID", "probeID")]
  names(snpAnnot) <- c("snpID", "snpName")

  # scan annotation
  data(affy_scan_annot)
  scanAnnot <- affy_scan_annot[,c("scanID", "genoRunID", "chpFile")]
  names(scanAnnot) <- c("scanID", "scanName", "file")

  # netCDF
  path <- system.file("extdata", "affy_raw_data", package="GWASdata")
  ncfile <- system.file("extdata", "affy_geno.nc", package="GWASdata")
  gdsfile <- tempfile()
  convertNcdfGds(ncfile, gdsfile)
  col.nums <- as.integer(c(2,3)); names(col.nums) <- c("snp", "geno")
  diagfile <- tempfile()
  res <- checkGenotypeFile(path, gdsfile, file.type="gds", snpAnnot, scanAnnot, sep.type="\t",
                         skip.num=1, col.total=6, col.nums=col.nums,
                         scan.name.in.file=-1, check.scan.index=1:3,
                         n.scans.loaded=3, diagnostics.filename=diagfile)
  checkTrue(all(res$chk == 1))
  checkTrue(all(res$snp.order == 1))
  checkTrue(all(res$geno.chk == 1))
  
  # now change a file and see if it detects the difference
  origfile <- system.file("extdata", "affy_raw_data",
                          scanAnnot$file[1], package="GWASdata")
  orig <- read.table(origfile, as.is=TRUE, header=TRUE, sep="\t")
  orig[1,3] <- 0
  tmpfile <- "tmp.txt"
  tmppath <- tempdir()
  tmp <- paste(tmppath, tmpfile, sep="/")
  write.table(orig, tmp, sep="\t", row.names=FALSE, quote=FALSE)
  scanAnnot$file[1] <- tmpfile
  res <- checkGenotypeFile(tmppath, ncfile, file.type="ncdf", snpAnnot, scanAnnot, sep.type="\t",
                         skip.num=1, col.total=6, col.nums=col.nums,
                         scan.name.in.file=-1, check.scan.index=1,
                         n.scans.loaded=3, diagnostics.filename=diagfile)
  checkTrue(res$geno.chk[1] == 0)

  write.table(orig[1:100,], tmp, sep="\t", row.names=FALSE, quote=FALSE)
  res <- checkGenotypeFile(tmppath, ncfile, file.type="ncdf", snpAnnot, scanAnnot, sep.type="\t",
                         skip.num=1, col.total=6, col.nums=col.nums,
                         scan.name.in.file=-1, check.scan.index=1,
                         n.scans.loaded=3, diagnostics.filename=diagfile)
  checkTrue(res$row.num[1] == 100)

  file.remove(tmp,diagfile)
}
