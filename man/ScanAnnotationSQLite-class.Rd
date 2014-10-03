\name{ScanAnnotationSQLite}
\docType{class}

\alias{ScanAnnotationSQLite-class}
\alias{ScanAnnotationSQLite}
\alias{hasVariable,ScanAnnotationSQLite-method}
\alias{getVariable,ScanAnnotationSQLite-method}
\alias{getScanID,ScanAnnotationSQLite-method}
\alias{hasSex,ScanAnnotationSQLite-method}
\alias{getSex,ScanAnnotationSQLite-method}
\alias{getVariableNames,ScanAnnotationSQLite-method}
\alias{getAnnotation,ScanAnnotationSQLite-method}
\alias{getMetadata,ScanAnnotationSQLite-method}
\alias{writeAnnotation,ScanAnnotationSQLite-method}
\alias{writeMetadata,ScanAnnotationSQLite-method}
\alias{getQuery,ScanAnnotationSQLite-method}
\alias{open,ScanAnnotationSQLite-method}
\alias{close,ScanAnnotationSQLite-method}
\alias{nscan,ScanAnnotationSQLite-method}
\alias{show,ScanAnnotationSQLite-method}

\title{Class ScanAnotationSQLite}

\description{
  The ScanAnnotationSQLite class stores annotation data associated with
  scans, as well as metadata describing each column, in an SQLite database.
}

\section{Constructor}{
  \describe{
    \item{}{
      \code{ScanAnnotationSQLite(dbpath)}:

      \code{dbpath} is the path to a SQLite database with tables
      "Annotation" and "Metadata."  "Annotation" must contain at least the
      following column:
      \itemize{
	\item "scanID": vector containing unique scan ids.
      }
      If a column representing sex is present, it must have the
      following format:
      \itemize{
	\item "sex": character vector with values 'M' or 'F'.
      }
      "Metadata" must contain at least the following columns:
      \itemize{
	\item "varname": name of variable in annotation
	\item "description": description of column in annotation
      }
      
      If the database does not yet exist, a database is created
      with tables "Annotation" and "Metadata."

      The \code{ScanAnnotationSQLite} constructor creates and returns
      a ScanAnnotationSQLite instance.
    }
  }
}

\section{Accessors}{
  In the code snippets below, \code{object} is a ScanAnnotationSQLite
  object.

  \describe{
    \item{}{
      \code{open(object)}: Opens a connection to the database.
    }
    \item{}{
      \code{close(object)}: Closes the database connection.
    }
    \item{}{
      \code{nscan(object)}: The number of scans in the database.
    }
    \item{}{
      \code{getScanID(object, index, condition)}: A unique vector of scan
      IDs.  The optional \code{index} is a logical or
      integer vector specifying elements to extract.    The optional
      \code{condition} is a character string with an SQL clause used to
      select data (e.g., "LIMIT 10", "WHERE sex='M'").
    }
    \item{}{
      \code{getSex(object, index, condition)}: A character vector of sex, with values 'M'
      or 'F'.    The optional \code{index} is a logical or
      integer vector specifying elements to extract.  The optional
      \code{condition} is a character string with an SQL clause used to
      select data.
    }
    \item{}{
      \code{hasSex(object)}: Returns \code{TRUE} if the column 'sex' is present in
      \code{object}. 
    }
    \item{}{
      \code{getVariable(object, varname, index, condition)}: A vector of the
      column \code{varname}.  The optional \code{index} is a logical or
      integer vector specifying elements to extract.  The optional
      \code{condition} is a character string with an SQL clause used to
      select data (e.g., "LIMIT 10", "WHERE sex='M'").
      Returns \code{NULL} if
      \code{varname} is not found in \code{object}.
    }
    \item{}{
      \code{hasVariable(object, varname)}: Returns \code{TRUE} if
      \code{varname} is a column in \code{object}, \code{FALSE} if not.
    }
    \item{}{
      \code{getVariableNames(object)}: Returns a character vector with
      the names of all columns in \code{object}.
    }
    \item{}{
      \code{getAnnotation(object)}: Returns all annotation variables
      as a data frame.
    }
    \item{}{
      \code{getMetadata(object)}: Returns metadata describing the
      annotation variables as a data frame.
    }
    \item{}{
      \code{getQuery(object, statement)}: Returns result of the SQL
      query \code{statement}.
    }
    \item{}{
      \code{writeAnnotation(object, value, append=FALSE,
	overwrite=TRUE)}: Writes \code{value} to the scan annotation
      table.  \code{value} must be a data.frame containing a column "scanID".
    }
    \item{}{
      \code{writeMetadata(object, value, append=FALSE,
	overwrite=TRUE)}: Writes \code{value} to the metadata table.
      \code{value} should be a data.frame containing 
       columns "varname" and "description".
    }
  }
}

\author{Stephanie Gogarten}

\seealso{
  \code{\link{SnpAnnotationSQLite}},
  \code{\link{ScanAnnotationDataFrame}},
  \code{\link{GenotypeData}}, \code{\link{IntensityData}}
}

\examples{
library(GWASdata)
dbpath <- tempfile()
scanAnnot <- ScanAnnotationSQLite(dbpath)

data(illumina_scan_annot)
writeAnnotation(scanAnnot, illumina_scan_annot)

# list columns
vars <- getVariableNames(scanAnnot)

# add metadata
metadf <- data.frame(varname=vars, description=rep(NA, length(vars)),
  row.names=vars, stringsAsFactors=FALSE)
metadf["scanID", "description"] <- "unique id"
writeMetadata(scanAnnot, metadf)

scanID <- getScanID(scanAnnot)
sex <- getSex(scanAnnot)
if (hasVariable(scanAnnot, "plate")) plate <- getVariable(scanAnnot, "plate")
subjectID <- getVariable(scanAnnot, "subjectID", condition="WHERE sex='M'")

# display data
head(getAnnotation(scanAnnot))
getMetadata(scanAnnot)

close(scanAnnot)
file.remove(dbpath)
}

\keyword{methods}
\keyword{classes}
    