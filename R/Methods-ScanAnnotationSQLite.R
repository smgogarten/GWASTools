# Methods for ScanAnnotationSQLite

# constructor
ScanAnnotationSQLite <- function(dbpath, ...) {
  dbcon <- dbConnect(SQLite(), dbname=dbpath)
  
  # if there are no tables yet, add them  
  scans.tables <- dbListTables(dbcon)
  if (!("Annotation" %in% scans.tables)) {
    sql <- "CREATE table Annotation (
              scanID integer primary key)"
    dbSendQuery(dbcon, sql)
  }
  if (!("Metadata" %in% scans.tables)) {
    metadata <- data.frame(varname=c("scanID"),
                           description=NA, stringsAsFactors=FALSE)
    dbWriteTable(dbcon, "Metadata", metadata, row.names=FALSE)
  }
  
  new("ScanAnnotationSQLite", dbpath=dbpath, dbcon=dbcon, ...)
}

# validity
setValidity("ScanAnnotationSQLite",
            function(object) {
              # check for required tables
              REQUIRED.TABLES <- c("Annotation", "Metadata")
              scans.tables <- dbListTables(object@dbcon)
              missing.tables <- setdiff(REQUIRED.TABLES, scans.tables)
              if (length(missing.tables) != 0) {
                msg <- c("tables not found in SQL database: ",
                         paste(missing.tables, collapse=", "))
                return(paste(msg, collapse=""))
              }
              
              # check for required columns
              REQUIRED.FIELDS <- c(object@idCol)
              scans.fields <- getVariableNames(object)
              missing.fields <- setdiff(REQUIRED.FIELDS, scans.fields)
              if (length(missing.fields) != 0) {
                msg <- c("fields not found in SQL table 'Annotation':",
                         paste(missing.fields, collapse=", "))
                return(paste(msg, collapse=""))
              }
              
              REQUIRED.META <- c("varname", "description")
              meta.fields <- dbListFields(object@dbcon, "Metadata")
              missing.fields <- setdiff(REQUIRED.META, meta.fields)
              if (length(missing.fields) != 0) {
                msg <- c("fields not found in SQL table 'Metadata': ",
                         paste(missing.fields, collapse=", "))
                return(paste(msg, collapse=""))
              }

              # check format of required columns
              # scanID should be a unique integer vector
              scanID <- getScanID(object)
              if (length(scanID) != length(unique(scanID)) ||
                  !is.integer(scanID)) {
                return(paste(object@idCol, "must be a unique integer vector"))
              }
              # sex should be M/F
              if (hasSex(object)) {
                sex <- getSex(object)
                if (!all(sex %in% c("M","F") | is.na(sex))) {
                  return(paste(object@sexCol, "should have values M/F"))
                }
              }

              # check that metadata matches annotation
              varnames <- dbGetQuery(object@dbcon, "SELECT varname from Metadata")[[1]]
              missing.names <- setdiff(scans.fields, varnames)
              if (length(missing.names) != 0) {
                msg <- c("variables not found in both annotation and metadata tables: ",
                         paste(missing.names, collapse=", "))
                return(paste(msg, collapse=""))
              }
              
              TRUE
            })

setMethod("show",
          signature(object = "ScanAnnotationSQLite"),
          function(object) {
            cat("An object of class", class(object), "\n")
            cat(paste("with", nscan(object), "scans and",
                      length(getVariableNames(object)), "variables\n"))
          })

setMethod("open",
          signature(con = "ScanAnnotationSQLite"),
          function(con, ...) {
            con@dbcon <- dbConnect(SQLite(), dbname=con@dbpath)
          })

setMethod("close",
          signature(con = "ScanAnnotationSQLite"),
          function(con, ...) {
            dbDisconnect(con@dbcon)
          })

setMethod("nscan", "ScanAnnotationSQLite",
          function(object) {
            sql <- "SELECT count(*) FROM Annotation"
            getQuery(object, sql)[[1]]
          })

setMethod("getVariableNames",
          signature(object = "ScanAnnotationSQLite"),
          function(object) {
            dbListFields(object@dbcon, "Annotation")
          })

setMethod("hasVariable",
          signature(object = "ScanAnnotationSQLite"),
          function(object, varname) {
            varname %in% getVariableNames(object)
          })

setMethod("getVariable",
          signature(object = "ScanAnnotationSQLite"),
          function(object, varname, index, condition) {
            # check that variable exists
            varexist <- hasVariable(object, varname)
            if (!all(varexist)) {
              warning(paste(paste(varname[!varexist], collapse=","), "not found"))
              return(NULL)
            }
            
            if (missing(condition)) {
              condition <- ""
            }
            
            sql <- paste("SELECT", paste(varname, collapse=","),
                         "FROM Annotation", condition)
            res <- getQuery(object, sql)
            # res is a data.frame
            # return only the first column for 1D data
            if (!missing(index) & ncol(res) == 1) {
              res[index,1]
            } else if (!missing(index) & ncol(res) > 1) {
              res[index,]
            } else if (missing(index) & ncol(res) == 1) {
              res[,1]
            } else {
              res
            }
          })

setMethod("getScanID",
          signature(object = "ScanAnnotationSQLite"),
          function(object, ...) {
            getVariable(object, object@idCol, ...)
          })
  
setMethod("hasSex",
          signature(object = "ScanAnnotationSQLite"),
          function(object) {
            hasVariable(object, object@sexCol)
          })
        
setMethod("getSex",
          signature(object = "ScanAnnotationSQLite"),
          function(object, ...) {
            getVariable(object, object@sexCol, ...)
          })

setMethod("getAnnotation",
          signature(object = "ScanAnnotationSQLite"),
          function(object, ...) {
            dbReadTable(object@dbcon, "Annotation",
                        row.names=NULL, ...)
          })

setMethod("getMetadata",
          signature(object = "ScanAnnotationSQLite"),
          function(object, ...) {
            dbReadTable(object@dbcon, "Metadata",
                        row.names=NULL, ...)
          })

setMethod("getQuery",
          signature(object = "ScanAnnotationSQLite"),
          function(object, statement, ...) {
            dbGetQuery(object@dbcon, statement, ...)
          })

# can specify overwrite=TRUE to add columns (must rewrite entire table)
setMethod("writeAnnotation",
          signature(object = "ScanAnnotationSQLite"),
          function(object, value, append=FALSE, overwrite=TRUE, ...) { 
            stopifnot(object@idCol %in% names(value))
            if (append) overwrite <- FALSE
            dbWriteTable(object@dbcon, "Annotation", value,
                         row.names=FALSE, append=append, overwrite=overwrite, ...)
          })

# can specify append=TRUE to add rows
setMethod("writeMetadata",
          signature(object = "ScanAnnotationSQLite"),
          function(object, value, append=FALSE, overwrite=TRUE, ...) {
            stopifnot(c("varname", "description") %in% names(value))
            if (append) overwrite <- FALSE
            dbWriteTable(object@dbcon, "Metadata", value,
                         row.names=FALSE, append=append, overwrite=overwrite, ...)
          })
