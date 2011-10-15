test_saveas <- function() { 
  x <- 1:10
  path <- tempdir()
  saveas(x, "myx", path)
  newfile <- file.path(path, "myx.RData")
  checkTrue(file.exists(newfile))
  load(newfile)
  checkTrue("myx" %in% objects())
  unlink(newfile)
  
  saveas(x, "myx.Rdata", path)
  newfile <- file.path(path, "myx.Rdata")
  checkTrue(file.exists(newfile))
  unlink(newfile)
  
  saveas(x, "myx.rda", path)
  newfile <- file.path(path, "myx.rda")
  checkTrue(file.exists(newfile))
  unlink(newfile)
}
