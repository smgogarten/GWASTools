test_getobj <- function() {
  x <- 1:10
  file <- tempfile()
  save(x, file=file)
  y <- getobj(file)
  checkIdentical(x,y)
  unlink(file)
}
