test_apartSnpSelection <- function() {
  chrom <- c(rep(1, 100), rep(27, 10))
  pos <- 1:length(chrom)
  res <- apartSnpSelection(chrom, pos, min.dist=10)
  checkEquals(res[101:110], rep(FALSE, 10))
  
  chrom <- c(rep(1, 100), rep(NA, 10))
  pos <- 1:length(chrom)
  res <- apartSnpSelection(chrom, pos, min.dist=10)
  checkEquals(res[101:110], rep(FALSE, 10))
}
