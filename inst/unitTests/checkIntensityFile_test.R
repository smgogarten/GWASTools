
test_checkIntensityFile_ncdf <- function() {
  # snp annotation
  data(affy_snp_annot)
  snpAnnot <- affy_snp_annot[,c("snpID", "probeID")]
  names(snpAnnot) <- c("snpID", "snpName")

  # scan annotation
  data(affy_scan_annot)
  scanAnnot <- affy_scan_annot[,c("scanID", "genoRunID", "chpFile")]
  names(scanAnnot) <- c("scanID", "scanName", "file")

  # quality
  path <- system.file("extdata", "affy_raw_data", package="GWASdata")
  ncfile <- system.file("extdata", "affy_qxy.nc", package="GWASdata")
  col.nums <- as.integer(c(2,4)); names(col.nums) <- c("snp", "quality")
  diagfile <- tempfile()
  res <- checkIntensityFile(path, ncfile, file.type="ncdf", snpAnnot, scanAnnot, sep.type="\t",
                          skip.num=1, col.total=6, col.nums=col.nums,
                          scan.name.in.file=-1, check.scan.index=1:3,
                          n.scans.loaded=3, affy.inten=FALSE,
                          diagnostics.filename=diagfile)
  checkTrue(all(res$chk == 1))
  checkTrue(all(res$snp.chk == 1))
  checkTrue(all(res$qs.chk == 1))

  # now change a file and see if it detects the difference
  origfile <- system.file("extdata", "affy_raw_data",
                          scanAnnot$file[1], package="GWASdata")
  orig <- read.table(origfile, as.is=TRUE, header=TRUE, sep="\t")
  orig[1,4] <- 0
  tmpfile <- "tmp.txt"
  tmppath <- tempdir()
  tmp <- paste(tmppath, tmpfile, sep="/")
  write.table(orig, tmp, sep="\t", row.names=FALSE, quote=FALSE)
  scanAnnot$file[1] <- tmpfile
  res <- checkIntensityFile(tmppath, ncfile, file.type="ncdf", snpAnnot, scanAnnot, sep.type="\t",
                          skip.num=1, col.total=6, col.nums=col.nums,
                          scan.name.in.file=-1, check.scan.index=1,
                          n.scans.loaded=3, affy.inten=FALSE,
                          diagnostics.filename=diagfile)
  checkTrue(res$qs.chk[1] == 0)

  # intensity
  scanAnnot <- affy_scan_annot[,c("scanID", "genoRunID", "alleleFile")]
  names(scanAnnot) <- c("scanID", "scanName", "file")
  res <- checkIntensityFile(path, ncfile, file.type="ncdf", snpAnnot, scanAnnot, sep.type="\t",
                          skip.num=1, col.total=6, col.nums=col.nums,
                          scan.name.in.file=-1, check.scan.index=1:3,
                          n.scans.loaded=3, affy.inten=TRUE,
                          diagnostics.filename=diagfile)
  checkTrue(all(res$chk == 1))
  checkTrue(all(res$snp.chk == 1))
  for (element in res$inten.chk) checkTrue(all(element == 1))
  
  # now change a file and see if it detects the difference
  origfile <- system.file("extdata", "affy_raw_data",
                          scanAnnot$file[1], package="GWASdata")
  orig <- read.table(origfile, as.is=TRUE, header=TRUE, sep="\t")
  orig[1,2] <- 0
  tmpfile2 <- "tmp2.txt"
  tmp2 <- paste(tmppath, tmpfile2, sep="/")
  write.table(orig, tmp2, sep="\t", row.names=FALSE, quote=FALSE)
  scanAnnot$file[1] <- tmpfile2
  res <- checkIntensityFile(tmppath, ncfile, file.type="ncdf", snpAnnot, scanAnnot, sep.type="\t",
                            skip.num=1, col.total=6, col.nums=col.nums,
                            scan.name.in.file=-1, check.scan.index=1,
                            n.scans.loaded=3, affy.inten=TRUE,
                            diagnostics.filename=diagfile)
  checkTrue(res$inten.chk$X[1] == 0)

  file.remove(tmp, tmp2, diagfile)
}

test_checkIntensityFile_gds <- function() {
  path <- system.file("extdata", "illumina_raw_data", package="GWASdata")
  data(illumina_snp_annot, illumina_scan_annot)
  snpAnnot <- illumina_snp_annot
  scanAnnot <- illumina_scan_annot
  scanAnnot <- scanAnnot[scanAnnot$file %in% list.files(path),]
  gdsfile <- tempfile()

  snpAnnot <- snpAnnot[,c("snpID", "rsID", "chromosome", "position", "alleleA", "alleleB")]
  names(snpAnnot)[1:2] <- c("snpID", "snpName")
  scanAnnot <- scanAnnot[,c("scanID", "genoRunID", "file")]
  names(scanAnnot) <- c("scanID", "scanName", "file")
  col.nums <- as.integer(c(1,2,5,16,17))
  names(col.nums) <- c("snp", "sample", "quality", "X", "Y")
  diagfile <- tempfile()
  res <- createDataFile(path, gdsfile, file.type="gds", variables=c("quality","X","Y"),
                        snpAnnot, scanAnnot, sep.type=",",
                        skip.num=11, col.total=21, col.nums=col.nums,
                        scan.name.in.file=1, diagnostics.filename=diagfile)
      
  res <- checkIntensityFile(path, gdsfile, file.type="gds",
                            snpAnnot, scanAnnot, sep.type=",",
                            skip.num=11, col.total=21, col.nums=col.nums,
                            scan.name.in.file=1, check.scan.index=1:3,
                            n.scans.loaded=3, affy.inten=FALSE,
                            diagnostics.filename=diagfile)
  checkTrue(all(res$chk == 1))
  checkTrue(all(res$snp.chk == 1))
  checkTrue(all(res$qs.chk == 1))
  for (element in res$inten.chk) checkTrue(all(element == 1))

  # now change a file and see if it detects the difference
  # quality
  origfile <- system.file("extdata", "illumina_raw_data",
                          scanAnnot$file[1], package="GWASdata")
  orig <- read.csv(origfile, as.is=TRUE, header=TRUE, skip=10)
  orig[1,5] <- 0
  tmpfile <- "tmp.txt"
  tmppath <- tempdir()
  tmp <- paste(tmppath, tmpfile, sep="/")
  write.csv(orig, tmp, row.names=FALSE, quote=FALSE)
  scanAnnot$file[1] <- tmpfile
  res <- checkIntensityFile(tmppath, gdsfile, file.type="gds",
                            snpAnnot, scanAnnot, sep.type=",",
                            skip.num=1, col.total=21, col.nums=col.nums,
                            scan.name.in.file=1, check.scan.index=1,
                            n.scans.loaded=3, affy.inten=FALSE,
                            diagnostics.filename=diagfile)
  checkTrue(res$qs.chk[1] == 0)

  # intensity
  orig <- read.csv(origfile, as.is=TRUE, header=TRUE, skip=10)
  orig[1,16] <- orig[1,16]*2
  write.csv(orig, tmp, row.names=FALSE, quote=FALSE)
  res <- checkIntensityFile(tmppath, gdsfile, file.type="gds",
                            snpAnnot, scanAnnot, sep.type=",",
                            skip.num=1, col.total=21, col.nums=col.nums,
                            scan.name.in.file=1, check.scan.index=1,
                            n.scans.loaded=3, affy.inten=FALSE,
                            diagnostics.filename=diagfile)
  checkTrue(res$inten.chk$X[1] == 0)

  file.remove(tmp, diagfile, gdsfile)
}
