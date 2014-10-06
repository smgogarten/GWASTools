test_GdsReader <- function() {
  file <- tempfile()
  gds <- createfn.gds(file)
  add.gdsn(gds, "var1", 1:100)
  add.gdsn(gds, "var2", letters)
  dat <- matrix(0, nrow=100, ncol=26)
  add.gdsn(gds, "var3", dat)
  closefn.gds(gds)

  obj <- GdsReader(file)
  checkIdentical(c("var1", "var2", "var3"), getVariableNames(obj))
  checkTrue(hasVariable(obj, "var1"))
  checkTrue(!hasVariable(obj, "foo"))
  checkIdentical(1:100, getVariable(obj, "var1"))
  checkIdentical(letters, getVariable(obj, "var2"))
  var3 <- getVariable(obj, "var3", start=c(2,1), count=c(5,10))
  checkIdentical(c(5L,10L), dim(var3))
  checkIdentical(dat[2:6,1:10], var3)
  checkIdentical(getVariable(obj, "var3"),
                 getVariable(obj, "var3", start=c(1,1), count=c(-1,-1)))
  checkIdentical(dim(dat), getDimension(obj, "var3"))

  checkTrue(!hasVariable(obj, "foo"))
  checkIdentical(NULL, getVariable(obj, "foo"))
  
  # check data types
  checkTrue(is(getVariable(obj, "var1"), "vector"))
  checkTrue(is(var3, "matrix"))
  
  
  # file errors
  checkException(GdsReader())
  checkException(GdsReader("foo")) 
  
  close(obj)
  
  # using an existing gds object
  gds <- openfn.gds(file)
  obj <- GdsReader(gds)
  checkIdentical(c("var1", "var2", "var3"), getVariableNames(obj))
  closefn.gds(gds)
  
  unlink(file)
}



test_GdsReader_attr <- function() {
  file <- tempfile()
  gds <- createfn.gds(file)
  node <- add.gdsn(gds, "var1", -1:10)
  put.attr.gdsn(node, "missing.value", -1)
  closefn.gds(gds)

  obj <- GdsReader(file)
  checkIdentical(-1, getAttribute(obj, "missing.value", "var1"))
  checkIdentical(c(NA,0:10), getVariable(obj, "var1"))

  close(obj)
  unlink(file)
}


test_GdsReader_folders <- function() {
  file <- tempfile()
  gds <- createfn.gds(file)
  add.gdsn(gds, "var1", 1:10)
  add.gdsn(gds, "var2", list(a=1:2, b=1:2))
  add.gdsn(gds, "var3", list(x=1:5, y=1:5, z=1:5))
  closefn.gds(gds)

  obj <- GdsReader(file)
  checkIdentical(c("var1", "var2/a", "var2/b", "var3/x", "var3/y", "var3/z"),
                 getVariableNames(obj))
  checkIdentical(1:2, getVariable(obj, "var2/b"))
  close(obj)
  unlink(file)
}
