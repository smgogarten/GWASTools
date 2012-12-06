# Methods for SnpAnnotationSQLite

# constructor
SnpAnnotationSQLite <- function(dbpath, ...) {
  dbcon <- dbConnect(SQLite(), dbname=dbpath)
  
  # if there are no tables yet, add them  
  snps.tables <- dbListTables(dbcon)
  if (!("Annotation" %in% snps.tables)) {
    sql <- "CREATE table Annotation (
              snpID integer primary key,
              chromosome integer,
              position integer)"
    dbSendQuery(dbcon, sql)
  }
  if (!("Metadata" %in% snps.tables)) {
    metadata <- data.frame(varname=c("snpID", "chromosome", "position"),
                           description=rep(NA, 3), stringsAsFactors=FALSE)
    dbWriteTable(dbcon, "Metadata", metadata, row.names=FALSE)
  }
  
  new("SnpAnnotationSQLite", dbpath=dbpath, dbcon=dbcon, ...)
}

# validity
setValidity("SnpAnnotationSQLite",
            function(object) {
              # check for required tables
              REQUIRED.TABLES <- c("Annotation", "Metadata")
              snps.tables <- dbListTables(object@dbcon)
              missing.tables <- setdiff(REQUIRED.TABLES, snps.tables)
              if (length(missing.tables) != 0) {
                msg <- c("tables not found in SQL database: ",
                         paste(missing.tables, collapse=", "))
                return(paste(msg, collapse=""))
              }
              
              # check for required columns
              REQUIRED.FIELDS <- c(object@idCol, object@chromosomeCol,
                                   object@positionCol)
              snps.fields <- getVariableNames(object)
              missing.fields <- setdiff(REQUIRED.FIELDS, snps.fields)
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
              # snpID should be a unique integer vector
              snpID <- getSnpID(object)
              if (length(snpID) != length(unique(snpID)) ||
                  !is.integer(snpID)) {
                return(paste(object@idCol, "must be a unique integer vector"))
              }
              # chromosome should be an integer vector
              # get only first 10 rows to save time
              if (!is.integer(getChromosome(object, condition="LIMIT 10"))) {
                return(paste(object@chromosomeCol, "must be an integer vector"))
              }
              # position should be an integer vector
              if (!is.integer(getPosition(object, condition="LIMIT 10"))) {
                return(paste(object@positionCol, "must be an integer vector"))
              }

              # check that metadata matches annotation
              varnames <- dbGetQuery(object@dbcon, "SELECT varname from Metadata")[[1]]
              missing.names <- setdiff(snps.fields, varnames)
              if (length(missing.names) != 0) {
                msg <- c("variables not found in both annotation and metadata tables: ",
                         paste(missing.names, collapse=", "))
                return(paste(msg, collapse=""))
              }
              
              TRUE
            })

setMethod("show",
          signature(object = "SnpAnnotationSQLite"),
          function(object) {
            cat("An object of class", class(object), "\n")
            cat(paste("with", nsnp(object), "snps and",
                      length(getVariableNames(object)), "variables\n"))
          })

setMethod("open",
          signature(con = "SnpAnnotationSQLite"),
          function(con, ...) {
            con@dbcon <- dbConnect(SQLite(), dbname=con@dbpath)
          })

setMethod("close",
          signature(con = "SnpAnnotationSQLite"),
          function(con, ...) {
            dbDisconnect(con@dbcon)
          })

setMethod("nsnp", "SnpAnnotationSQLite",
          function(object) {
            sql <- "SELECT count(*) FROM Annotation"
            getQuery(object, sql)[[1]]
          })

setMethod("getVariableNames",
          signature(object = "SnpAnnotationSQLite"),
          function(object) {
            dbListFields(object@dbcon, "Annotation")
          })

setMethod("hasVariable",
          signature(object = "SnpAnnotationSQLite"),
          function(object, varname) {
            varname %in% getVariableNames(object)
          })

setMethod("getVariable",
          signature(object = "SnpAnnotationSQLite"),
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

setMethod("getSnpID",
          signature(object = "SnpAnnotationSQLite"),
          function(object, ...) {
            getVariable(object, object@idCol, ...)
          })

setMethod("getChromosome",
          signature(object = "SnpAnnotationSQLite"),
          function(object, char=FALSE, ...) {
            var <- getVariable(object, object@chromosomeCol, ...)
            
            # convert to characters
            if (char) {
              # default is unknown code
              chromChar <- rep("U", length(var))
              autosome <- var %in% object@autosomeCode
              chromChar[autosome] <- as.character(var[autosome])
              xchrom <- var == object@XchromCode & !is.na(var)
              chromChar[xchrom] <- "X"
              ychrom <- var == object@YchromCode & !is.na(var)
              chromChar[ychrom] <- "Y"
              xychrom <- var == object@XYchromCode & !is.na(var)
              chromChar[xychrom] <- "XY"
              mchrom <- var == object@MchromCode & !is.na(var)
              chromChar[mchrom] <- "M"
              var <- chromChar
            }
            var
          })

setMethod("getPosition",
          signature(object = "SnpAnnotationSQLite"),
          function(object, ...) {
            getVariable(object, object@positionCol, ...)
          })

setMethod("getAlleleA",
          signature(object = "SnpAnnotationSQLite"),
          function(object, ...) {
            getVariable(object, object@alleleACol, ...)
          })

setMethod("getAlleleB",
          signature(object = "SnpAnnotationSQLite"),
          function(object, ...) {
            getVariable(object, object@alleleBCol, ...)
          })

setMethod("getAnnotation",
          signature(object = "SnpAnnotationSQLite"),
          function(object, ...) {
            dbReadTable(object@dbcon, "Annotation",
                        row.names=NULL, ...)
          })

setMethod("getMetadata",
          signature(object = "SnpAnnotationSQLite"),
          function(object, ...) {
            dbReadTable(object@dbcon, "Metadata",
                        row.names=NULL, ...)
          })

setMethod("getQuery",
          signature(object = "SnpAnnotationSQLite"),
          function(object, statement, ...) {
            dbGetQuery(object@dbcon, statement, ...)
          })

# can specify overwrite=TRUE to add columns (must rewrite entire table)
setMethod("writeAnnotation",
          signature(object = "SnpAnnotationSQLite"),
          function(object, value, append=FALSE, overwrite=TRUE, ...) { 
            stopifnot(c(object@idCol, object@chromosomeCol, object@positionCol) %in%
                      names(value))
            if (append) overwrite <- FALSE
            dbWriteTable(object@dbcon, "Annotation", value,
                         row.names=FALSE, append=append, overwrite=overwrite, ...)
          })

# can specify append=TRUE to add rows
setMethod("writeMetadata",
          signature(object = "SnpAnnotationSQLite"),
          function(object, value, append=FALSE, overwrite=TRUE, ...) {
            stopifnot(c("varname", "description") %in% names(value))
            if (append) overwrite <- FALSE
            dbWriteTable(object@dbcon, "Metadata", value,
                         row.names=FALSE, append=append, overwrite=overwrite, ...)
          })

setMethod("autosomeCode", "SnpAnnotationSQLite",
          function(object) {
            object@autosomeCode
          })
      
setMethod("XchromCode", "SnpAnnotationSQLite",
          function(object) {
            object@XchromCode
          })
          
setMethod("YchromCode", "SnpAnnotationSQLite",
          function(object) {
            object@YchromCode
          })
              
setMethod("XYchromCode", "SnpAnnotationSQLite",
          function(object) {
            object@XYchromCode
          })
              
setMethod("MchromCode", "SnpAnnotationSQLite",
          function(object) {
            object@MchromCode
          })
          
