#load an Rdata file and return the object
#(using load directly returns only the object name)
getobj <- function(Rdata) {
  objname <- load(Rdata)
  if (length(objname) > 1) {
    warning(paste("Multiple objects stored in file", Rdata,
                  "\nReturning only the first object"))
  }
  return(get(objname))
}

#save an R object obj as name in an Rdata file called path/name.RData
saveas <- function(obj, name, path=".") {
  assign(name, obj)
  if (grepl("[.]RData$", name, ignore.case=TRUE) |
      grepl("[.]rda$", name, ignore.case=TRUE)) {
    file <- file.path(path, name)
  }
  else {
    file <- file.path(path, paste(name, "RData", sep="."))
  }
  save(list=name, file=file)
}

#returns TRUE if x and y exist and all elements are equal
#if there are NA values, returns TRUE if is.na(x) == is.na(y) and
#all other elements are equal, otherwise returns FALSE
#retuns FALSE if x or y (but not both) is NULL
allequal <- function(x, y) {
  if ((is.null(x) & !is.null(y)) | (!is.null(x) & is.null(y))) return(FALSE)

  if (!all(is.na(x) == is.na(y))) return(FALSE)

  return(all(x[!is.na(x)] == y[!is.na(y)]))
}

#####
# Read and write the first n lines of a file
#####
readWriteFirst <- function (filein, fileout, n)
{
	# read first n lines of filein and write them to fileout, where filein and fileout are file names
	incon <- file(filein,"r")  # open connection for reading
	outcon <- file(fileout,"w")  # open for writing
	for(i in 1:n) 
		write(readLines(incon, n=1), file=outcon, append=TRUE)
	close(incon)
	close(outcon)
}

