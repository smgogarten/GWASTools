test_ncdfCheckIntensity <- function() {
  # snp annotation
  data(affy_snp_annot)
  snpAnnot <- affy_snp_annot[,c("snpID", "probeID")]
  names(snpAnnot) <- c("snpID", "snpName")

  # scan annotation
  data(affy_scan_annot)
  scanAnnot <- affy_scan_annot[,c("scanID", "genoRunID", "chpFile", "alleleFile")]
  names(scanAnnot) <- c("scanID", "scanName", "file", "inten.file")

  # netCDF
  path <- system.file("extdata", "affy_raw_data", package="GWASdata")
  intenpath <- system.file("extdata", "affy_raw_data", package="GWASdata")
  ncfile <- system.file("extdata", "affy_qxy.nc", package="GWASdata")
  col.nums <- as.integer(c(2,4)); names(col.nums) <- c("snp", "qs")
  diagfile <- tempfile()
  res <- ncdfCheckIntensity(path, intenpath, ncfile, snpAnnot, scanAnnot, sep.type="\t",
                          skip.num=1, col.total=6, col.nums=col.nums,
                          scan.name.in.file=-1, check.scan.index=1:3,
                          n.scans.loaded=3, affy.inten=TRUE,
                          diagnostics.filename=diagfile)
  checkTrue(all(res$chk == 1))
  checkTrue(all(res$snp.order == 1))
  checkTrue(all(res$qs.chk == 1))
  for (element in res$inten.chk) checkTrue(all(element == 1))

  # now change a file and see if it detects the difference
  # quality
  origfile <- system.file("extdata", "affy_raw_data",
                          scanAnnot$file[1], package="GWASdata")
  orig <- read.table(origfile, as.is=TRUE, header=TRUE, sep="\t")
  orig[1,4] <- 0
  tmpfile <- "tmp.txt"
  tmppath <- tempdir()
  tmp <- paste(tmppath, tmpfile, sep="/")
  write.table(orig, tmp, sep="\t", row.names=FALSE, quote=FALSE)
  scanAnnot$file[1] <- tmpfile
  res <- ncdfCheckIntensity(tmppath, intenpath, ncfile, snpAnnot, scanAnnot, sep.type="\t",
                          skip.num=1, col.total=6, col.nums=col.nums,
                          scan.name.in.file=-1, check.scan.index=1,
                          n.scans.loaded=3, affy.inten=TRUE,
                          diagnostics.filename=diagfile)
  checkTrue(res$qs.chk[1] == 0)

  # intensity
  scanAnnot <- affy_scan_annot[,c("scanID", "genoRunID", "chpFile", "alleleFile")]
  names(scanAnnot) <- c("scanID", "scanName", "file", "inten.file")
  origfile <- system.file("extdata", "affy_raw_data",
                          scanAnnot$inten.file[1], package="GWASdata")
  orig <- read.table(origfile, as.is=TRUE, header=TRUE, sep="\t")
  orig[1,2] <- 0
  tmpfile2 <- "tmp2.txt"
  tmp2 <- paste(tmppath, tmpfile2, sep="/")
  write.table(orig, tmp2, sep="\t", row.names=FALSE, quote=FALSE)
  scanAnnot$inten.file[1] <- tmpfile2
  checkException({
    res <- ncdfCheckIntensity(path, tmppath, ncfile, snpAnnot, scanAnnot, sep.type="\t",
                          skip.num=1, col.total=6, col.nums=col.nums,
                          scan.name.in.file=-1, check.scan.index=1,
                          n.scans.loaded=3, affy.inten=TRUE,
                          diagnostics.filename=diagfile)
  })

  file.remove(tmp, tmp2, diagfile)
}


test_checkIntensityFile <- function() {
  # snp annotation
  data(affy_snp_annot)
  snpAnnot <- affy_snp_annot[,c("snpID", "probeID")]
  names(snpAnnot) <- c("snpID", "snpName")

  # scan annotation
  data(affy_scan_annot)
  scanAnnot <- affy_scan_annot[,c("scanID", "genoRunID", "chpFile", "alleleFile")]
  names(scanAnnot) <- c("scanID", "scanName", "file", "inten.file")

  # netCDF
  path <- system.file("extdata", "affy_raw_data", package="GWASdata")
  intenpath <- system.file("extdata", "affy_raw_data", package="GWASdata")
  ncfile <- system.file("extdata", "affy_qxy.nc", package="GWASdata")
  col.nums <- as.integer(c(2,4)); names(col.nums) <- c("snp", "quality")
  diagfile <- tempfile()
  res <- checkIntensityFile(path, intenpath, ncfile, file.type="ncdf", snpAnnot, scanAnnot, sep.type="\t",
                          skip.num=1, col.total=6, col.nums=col.nums,
                          scan.name.in.file=-1, check.scan.index=1:3,
                          n.scans.loaded=3, affy.inten=TRUE,
                          diagnostics.filename=diagfile)
  checkTrue(all(res$chk == 1))
  checkTrue(all(res$snp.order == 1))
  checkTrue(all(res$qs.chk == 1))
  for (element in res$inten.chk) checkTrue(all(element == 1))

  # now change a file and see if it detects the difference
  # quality
  origfile <- system.file("extdata", "affy_raw_data",
                          scanAnnot$file[1], package="GWASdata")
  orig <- read.table(origfile, as.is=TRUE, header=TRUE, sep="\t")
  orig[1,4] <- 0
  tmpfile <- "tmp.txt"
  tmppath <- tempdir()
  tmp <- paste(tmppath, tmpfile, sep="/")
  write.table(orig, tmp, sep="\t", row.names=FALSE, quote=FALSE)
  scanAnnot$file[1] <- tmpfile
  res <- checkIntensityFile(tmppath, intenpath, ncfile, file.type="ncdf", snpAnnot, scanAnnot, sep.type="\t",
                          skip.num=1, col.total=6, col.nums=col.nums,
                          scan.name.in.file=-1, check.scan.index=1,
                          n.scans.loaded=3, affy.inten=TRUE,
                          diagnostics.filename=diagfile)
  checkTrue(res$qs.chk[1] == 0)

  # intensity
  scanAnnot <- affy_scan_annot[,c("scanID", "genoRunID", "chpFile", "alleleFile")]
  names(scanAnnot) <- c("scanID", "scanName", "file", "inten.file")
  origfile <- system.file("extdata", "affy_raw_data",
                          scanAnnot$inten.file[1], package="GWASdata")
  orig <- read.table(origfile, as.is=TRUE, header=TRUE, sep="\t")
  orig[1,2] <- 0
  tmpfile2 <- "tmp2.txt"
  tmp2 <- paste(tmppath, tmpfile2, sep="/")
  write.table(orig, tmp2, sep="\t", row.names=FALSE, quote=FALSE)
  scanAnnot$inten.file[1] <- tmpfile2
  checkException({
    res <- checkIntensityFile(path, tmppath, ncfile, file.type="ncdf", snpAnnot, scanAnnot, sep.type="\t",
                          skip.num=1, col.total=6, col.nums=col.nums,
                          scan.name.in.file=-1, check.scan.index=1,
                          n.scans.loaded=3, affy.inten=TRUE,
                          diagnostics.filename=diagfile)
  })

  file.remove(tmp, tmp2, diagfile)
}
