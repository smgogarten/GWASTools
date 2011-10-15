\name{mendelList}
\alias{mendelList}
\alias{mendelListAsDataFrame}
\title{Mendelian Error Checking}
\description{
  \code{mendelList} creates a "mendelList" object (a list of trios).
  \code{mendelListAsDataFrame} converts a "mendelList" object to a data frame.
}

\usage{
mendelList(familyid, offspring, father, mother, sex, scanID)

mendelListAsDataFrame(mendel.list)
}

\arguments{
  \item{familyid}{A vector of family ID numbers.}
  \item{offspring}{A vector of offspring ID numbers.}
  \item{father}{A vector of father ID numbers.}
  \item{mother}{A vector of mother ID numbers.}
  \item{sex}{A vector to specify whether each scan is male "M" or female "F".}
  \item{scanID}{A vector of unique scan identification numbers.  These will be used to identify scans in output.}
  
  \item{mendel.list}{An object of class "mendelList".}
}

\details{
	The lengths of \code{familyid}, \code{offspring}, \code{father},
	\code{mother}, \code{sex}, and \code{scanID} must all be identical.
	The "mendelList" object is required as input for the \code{\link{mendelErr}} function.
}

\value{
  \code{mendelList} returns a "mendelList" object.  A "mendelList" object is a list of lists.
  The first level list is all the families.
  The second level list is offspring within families who have one or both parents genotyped.
  Within the second level are data.frame(s) with columns "offspring", "father", and "mother" which each contain
  the scanID for each member of the trio (a missing parent is denoted by -1).
  When replicates of the same offsping ID occur (duplicate scans for the same subject),
  this data.frame has multiple rows representing all combinations
  of scanIDs for that trio.

  \code{mendelListAsDataFrame}  returns a data.frame with variables "offspring", "father", and "mother" which each contain
the scanID for each member of the trio (a missing parent is denoted by -1).  This takes every data.frame
from the "mendelList" object and puts them all into one large data frame.  This can be easier to work with
for certain analyses.
}

\author{Xiuwen Zheng, Matt Conomos}

\seealso{\code{\link{mendelErr}}}

\examples{
# get sample annotation
library(GWASdata)
data(affy_scan_annot)

# generate trio list
men.list <- with(affy_scan_annot, mendelList(family, subjectID, father, mother, sex, scanID))

class(men.list)
# [1] "mendelList"

# convert into a data.frame
men.df <- mendelListAsDataFrame(men.list)

class(men.df)
# [1] "data.frame"
}

\keyword{Mendelian}

