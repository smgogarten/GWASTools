\name{SnpAnnotationSQLite}
\docType{class}

\alias{SnpAnnotationSQLite-class}
\alias{SnpAnnotationSQLite}
\alias{hasVariable,SnpAnnotationSQLite-method}
\alias{getVariable,SnpAnnotationSQLite-method}
\alias{getSnpID,SnpAnnotationSQLite-method}
\alias{getChromosome,SnpAnnotationSQLite-method}
\alias{getPosition,SnpAnnotationSQLite-method}
\alias{getVariableNames,SnpAnnotationSQLite-method}
\alias{getAnnotation,SnpAnnotationSQLite-method}
\alias{getMetadata,SnpAnnotationSQLite-method}
\alias{writeAnnotation,SnpAnnotationSQLite-method}
\alias{writeMetadata,SnpAnnotationSQLite-method}
\alias{getQuery,SnpAnnotationSQLite-method}
\alias{open,SnpAnnotationSQLite-method}
\alias{close,SnpAnnotationSQLite-method}
\alias{nsnp,SnpAnnotationSQLite-method}
\alias{show,SnpAnnotationSQLite-method}
\alias{XchromCode,SnpAnnotationSQLite-method}
\alias{XYchromCode,SnpAnnotationSQLite-method}
\alias{YchromCode,SnpAnnotationSQLite-method}
\alias{MchromCode,SnpAnnotationSQLite-method}

\title{Class SnpAnotationSQLite}

\description{
  The SnpAnnotationSQLite class stores annotation data associated with
  SNPs, as well as metadata describing each column, in an SQLite database.
}

\section{Constructor}{
  \describe{
    \item{}{
      \code{SnpAnnotationSQLite(dbpath)}:

      \code{dbpath} is the path to a SQLite database with tables
      "Annotation" and "Metadata."  "Annotation" must contain at least the
      following columns:
      \itemize{
	\item "snpID": integer vector containing unique SNP ids.
	\item "chromosome": integer vector containing chromosome codes.
	\item "position": integer vector containing position (in base
	pairs) on the chromosome.
      }
      Default values for chromosome codes are 1-22, 23=X, 24=XY, 25=Y,
      26=M.  The defaults may be changed with the arguments
      \code{XchromCode}, \code{XYchromCode}, \code{YchromCode}, and
      \code{MchromCode}.
      
      "Metadata" must contain at least the following columns:
      \itemize{
	\item "varname": name of variable in annotation
	\item "description": description of column in annotation
      }
      
      If the database does not yet exist, a database is created
      with tables "Annotation" and "Metadata."

      The \code{SnpAnnotationSQLite} constructor creates and returns
      a SnpAnnotationSQLite instance.
    }
  }
}

\section{Accessors}{
  In the code snippets below, \code{object} is a SnpAnnotationSQLite
  object.

  \describe{
    \item{}{
      \code{open(object)}: Opens a connection to the database.
    }
    \item{}{
      \code{close(object)}: Closes the database connection.
    }
    \item{}{
      \code{nsnp(object)}: The number of SNPs in the database.
    }
    \item{}{
      \code{getSnpID(object, index, condition)}: A unique integer vector of snp
      IDs.  The optional \code{index} is a logical or
      integer vector specifying elements to extract.    The optional
      \code{condition} is a character string with an SQL clause used to
      select data (e.g., "LIMIT 10", "WHERE chromosome=1").
    }
    \item{}{
      \code{getChromosome(object, index, condition, char=FALSE)}: A vector of
      chromosomes.  The optional \code{index} is a logical or
      integer vector specifying elements to extract.  The optional
      \code{condition} is a character string with an SQL clause used to
      select data (e.g., "LIMIT 10", "WHERE chromosome=1").
      If \code{char=FALSE} (default), returns an integer vector.
      If \code{char=TRUE}, returns a character vector with elements in
      (1:22,X,XY,Y,M,U).  "U" stands for "Unknown" and is the value
      given to any chromosome code not falling in the other categories.
    }
    \item{}{
      \code{getPosition(object, index, condition)}: An integer vector of base pair
      positions.  The optional \code{index} is a logical or
      integer vector specifying elements to extract.    The optional
      \code{condition} is a character string with an SQL clause used to
      select data (e.g., "LIMIT 10", "WHERE chromosome=1").
    }
    \item{}{
      \code{getVariable(object, varname, index, condition)}: A vector of the
      column \code{varname}.  The optional \code{index} is a logical or
      integer vector specifying elements to extract.  The optional
      \code{condition} is a character string with an SQL clause used to
      select data (e.g., "LIMIT 10", "WHERE chromosome=1").
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
	overwrite=TRUE)}: Writes \code{value} to the SNP annotation
      table.  \code{value} must be a data.frame containing columns
      "snpID", "chromosome", and "position".
    }
    \item{}{
      \code{writeMetadata(object, value, append=FALSE,
	overwrite=TRUE)}: Writes \code{value} to the metadata table.
      \code{value} should be a data.frame containing
       columns "varname" and "description".
     }
     \item{}{
       \code{XchromCode(object)}: Returns the integer code for the X
       chromosome.
     }
     \item{}{
       \code{XYchromCode(object)}: Returns the integer code for the
       pseudoautosomal region.
     }
     \item{}{
       \code{YchromCode(object)}: Returns the integer code for the Y
       chromosome.
     }
     \item{}{
       \code{MchromCode(object)}: Returns the integer code for 
       mitochondrial SNPs.
     }
   }
 }

\author{Stephanie Gogarten}

\seealso{
  \code{\link{ScanAnnotationSQLite}},
  \code{\link{SnpAnnotationDataFrame}},
  \code{\link{GenotypeData}}, \code{\link{IntensityData}}
}

\examples{
library(GWASdata)
dbpath <- tempfile()
snpAnnot <- SnpAnnotationSQLite(dbpath)

data(affy_snp_annot)
writeAnnotation(snpAnnot, affy_snp_annot)

# list columns
vars <- getVariableNames(snpAnnot)

# add metadata
metadf <- data.frame(varname=vars, description=rep(NA, length(vars)),
  row.names=vars, stringsAsFactors=FALSE)
metadf["snpID", "description"] <- "integer id"
writeMetadata(snpAnnot, metadf)

# get snpID and chromosome
snpID <- getSnpID(snpAnnot)
chrom <- getChromosome(snpAnnot)

# get positions only for chromosome 22
pos22 <- getPosition(snpAnnot, condition="WHERE chromosome = 22")

# get rsID
if (hasVariable(snpAnnot, "rsID")) rsID <- getVariable(snpAnnot, "rsID")

# display data
head(getAnnotation(snpAnnot))
getMetadata(snpAnnot)

close(snpAnnot)
file.remove(dbpath)
}

\keyword{methods}
\keyword{classes}
    