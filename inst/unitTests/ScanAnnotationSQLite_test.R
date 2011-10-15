test_ScanAnnotationSQLite <- function() {
  # create new database
  dbpath <- tempfile()  
  obj <- ScanAnnotationSQLite(dbpath)
  
  # test object
  scanID <- 1:10
  sex <- sample(c("M","F"), 10, replace=TRUE)
  plate <- letters[1:10]
  x <- data.frame(scanID=scanID, sex=sex, plate=plate,
                  stringsAsFactors=FALSE)
  writeAnnotation(obj, x)
  checkIdentical(x, getAnnotation(obj))
  checkEquals(10, nscan(obj))
  metadata <- data.frame(varname=c("scanID", "sex", "plate"),
                         description=c("id", "sex", "plate"),
                         stringsAsFactors=FALSE)
  writeMetadata(obj, metadata)
  checkIdentical(metadata, getMetadata(obj))

  # test errors
  xx <- data.frame(sex=sex, plate=plate)
  checkException(writeAnnotation(obj, xx))
  xmeta <- data.frame(var="scanID", descr="id")
  checkException(writeMetadata(obj, xmeta))

  # add a column
  subjID <- paste("subj", 1:10, sep="")
  x <- cbind(x, subjID, stringsAsFactors=FALSE)
  writeAnnotation(obj, x)
  checkIdentical(x, getAnnotation(obj))
  checkException(validObject(obj))
  newmeta <- data.frame(varname="subjID", description="subj id",
                         stringsAsFactors=FALSE)
  writeMetadata(obj, newmeta, append=TRUE)
  checkIdentical(rbind(metadata, newmeta), getMetadata(obj))

  # retrieve named columns
  checkIdentical(scanID, getScanID(obj))
  checkIdentical(sex, getSex(obj))
  checkTrue(hasSex(obj))

  # other columns
  checkTrue(hasVariable(obj, "subjID"))
  checkIdentical(subjID, getVariable(obj,"subjID"))
  checkTrue(!hasVariable(obj, "foo"))
  checkIdentical(NULL, getVariable(obj, "foo"))
  vars <- c("scanID", "subjID")
  checkIdentical(x[,vars], getVariable(obj, vars))

  # check indexing
  sel <- 1:5 # numeric
  checkIdentical(scanID[sel], getScanID(obj, index=sel))
  sel <- scanID > 5 # logical
  checkIdentical(scanID[sel], getScanID(obj, index=sel))

  # check condition
  cond <- "LIMIT 10"
  checkIdentical(scanID[1:10], getScanID(obj, condition=cond))
  cond <- "WHERE sex='M'"
  checkIdentical(scanID[sex == "M"], getScanID(obj, condition=cond))

  close(obj)
  file.remove(dbpath)
}
