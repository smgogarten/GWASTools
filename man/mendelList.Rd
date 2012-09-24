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
  \item{familyid}{A vector of family identifiers.}
  \item{offspring}{A vector of offspring subject identifiers.}
  \item{father}{A vector of father identifiers.}
  \item{mother}{A vector of mother identifiers.}
  \item{sex}{A vector to specify whether each subject is male "M" or female "F".}
  \item{scanID}{A vector of scanIDs indicating unique genotyping
  instances for the \code{offspring} vector.  In the case of duplicate samples, the same
  \code{offspring} identifier may correspond to multiple \code{scanID} values.}
  
  \item{mendel.list}{An object of class "mendelList".}
}

\details{
  The lengths of \code{familyid}, \code{offspring}, \code{father},
  \code{mother}, \code{sex}, and \code{scanID} must all be
  identical.  These vectors should include all genotyped samples, i.e.,
  samples present in the \code{father} and \code{mother} vectors should also
  appear in the \code{offspring} vector if there are genotypes for these
  samples, and their unique scan IDs should be given in the
  \code{scanID} vector.
  
  Identifiers may be character strings or integers, but not factors.
	
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

\author{Xiuwen Zheng, Matthew P. Conomos}

\seealso{\code{\link{mendelErr}}}

\examples{
# data frame of sample information.  No factors!
dat <- data.frame(family=c(1,1,1,1,2,2,2), offspring=c("a","a","b","c","d","e","f"),
  father=c("b","b",0,0,"e",0,0), mother=c("c","c",0,0,"f",0,0),
  sex=c("M","M","M","F","F","M","F"), scanID=1:7,
  stringsAsFactors=FALSE)
dat

men.list <- mendelList(dat$family, dat$offspring, dat$father, dat$mother,
                       dat$sex, dat$scanID)
men.list

# If fathers and mothers do not have separate entries in each vector,
# mendelList returns a "NULL":
dat <- dat[c(1,5),]
dat
mendelList(dat$family, dat$offspring, dat$father, dat$mother,
           dat$sex, dat$scanID)

men.df <- mendelListAsDataFrame(men.list)
men.df
}

\keyword{Mendelian}

