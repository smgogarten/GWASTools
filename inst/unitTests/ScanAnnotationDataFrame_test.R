test_ScanAnnotationDataFrame <- function() {
  # test object
  scanID <- 1:10
  sex <- sample(c("M","F"), 10, replace=TRUE)
  plate <- letters[1:10]
  x <- data.frame(scanID=scanID, sex=sex, plate=plate,
                  stringsAsFactors=FALSE)
  obj <- ScanAnnotationDataFrame(x)
  checkIdentical(x, as(obj, "data.frame"))

  # retrieve named columns
  checkIdentical(scanID, getScanID(obj))
  checkTrue(hasSex(obj))
  checkIdentical(sex, getSex(obj))

  # other columns
  checkTrue(hasVariable(obj, "plate"))
  checkIdentical(plate, getVariable(obj,"plate"))
  checkTrue(!hasVariable(obj, "foo"))
  checkIdentical(NULL, getVariable(obj, "foo"))
  vars <- c("scanID", "plate")
  checkIdentical(x[,vars], getVariable(obj, vars))

  # check indexing
  sel <- 1:5 # numeric
  checkIdentical(scanID[sel], getScanID(obj, sel))
  sel <- scanID > 5 # logical
  checkIdentical(scanID[sel], getScanID(obj, sel))

  # missing required columns
  x <- data.frame(sex=sex, plate=plate)
  checkException(ScanAnnotationDataFrame(x))

  # bad column format
  x <- data.frame(scanID=as.character(scanID))
  checkException(ScanAnnotationDataFrame(x))
  x <- data.frame(scanID=scanID, sex=rep("A",10))
  checkException(ScanAnnotationDataFrame(x))

  # scanID not unique
  x <- data.frame(scanID=rep(1L,10))
  checkException(ScanAnnotationDataFrame(x))

  # get and set methods
  x <- data.frame(scanID=scanID, sex=sex)
  pData(obj) <- x
  checkIdentical(x, getAnnotation(obj))
  meta <- getMetadata(obj)
  meta["scanID", "labelDescription"] <- "id"
  varMetadata(obj) <- meta
  checkIdentical(meta, getMetadata(obj))
}
