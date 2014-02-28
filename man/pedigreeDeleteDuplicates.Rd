\name{pedigreeDeleteDuplicates}
\alias{pedigreeDeleteDuplicates}

\title{
Remove duplicates from a pedigree
}
\description{
\code{pedigreeDeleteDuplicates} removes duplicates
from a pedigree.
}
\usage{
pedigreeDeleteDuplicates(pedigree, duplicates)
}


  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\arguments{
  \item{pedigree}{A dataframe containing the pedigree information
    for the samples to be examined with columns labeled "family",
    "individ", "mother", "father" and "sex" containing the identifiers
    of the family, individual, individual's mother, individual's
    father and individual's sex (coded as "M" or "F") .}

  \item{duplicates}{dataframe with columns "family" (family id) and "individ" (individual id).}

  \item{verbose}{Logical value specifying whether or not to display messages.}	

}



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\details{
  The output of \code{\link{pedigreeCheck}} can be provided
  to \code{pedigreeDeleteDuplicates} in order to generate a new
  pedigree with duplicates removed.
}
\value{
The output of \code{pedigreeFindDuplicates} is a list containing two dataframes:
\item{dups.mismatch}{A dataframe containing the family id, individual id and number of copies for any duplicates with mismatching pedigree data}
\item{dups.match}{A dataframe containing the family id, individual id and number of copies for any duplicates with matching pedigree data.}

The output of \code{pedigreeDeleteDuplicates} is a pedigree
identical to \code{pedigree}, but with duplicates removed.
}



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\author{Cecelia Laurie
}


\seealso{\code{\link{pedigreeCheck}}, \code{\link{pedigreePairwiseRelatedness}}
}

\examples{
family <- c(1,1,1,1,2,2,2,2)
individ <- c(1,2,3,3,4,5,6,6)
mother <- c(0,0,1,1,0,0,4,4)
father <- c(0,0,2,2,0,0,5,5)
sex <- c("F","M","F","F","F","F","M","M")
pedigree <- data.frame(family, individ, mother, father, sex, stringsAsFactors=FALSE)
duplicates <- pedigreeCheck(pedigree)$duplicates
pedigree.no.dups <- pedigreeDeleteDuplicates(pedigree, duplicates) 
}

\keyword{manip}
