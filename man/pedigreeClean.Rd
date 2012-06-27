\name{pedigreeClean}
\alias{pedigreeClean}

\title{Basic pedigree data checking
}

\description{This function checks a pedigree for completeness and gross errors. 
}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\usage{
pedigreeClean(pedigree)
}

\arguments{

\item{pedigree}{A dataframe containing the pedigree information for the samples to be examined with columns labeled "family", "individ", "mother", "father" and "sex" containing the identifiers of the family, individual, individual's mother, individual's father and individual's sex (coded as "M" or "F"). 
  Identifiers can be integer, numeric or character but identifiers for mother and father for founders are assumed to be 0. Identifier for unknown parent is assumed to be 0.}

}


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\details{

The function performs a basic check on pedigree data for gross errors, checking for missing (or blank) id's, mis-coded sex, individual id's of 0 and for individuals that appear as both mothers and fathers. 

}


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\value{
A list with the following components:
\item{family.missing}{A vector of integers containing the row positions of entries where family id's are missing (\code{NA}) or blank}
\item{individ.missing}{A vector of integers containing the row positions of entries where individual id's are missing (\code{NA}) or blank}
\item{father.missing}{A vector of integers containing the row positions of entries where father id's are missing (\code{NA}) or blank}
\item{mother.missing}{A vector of integers containing the row positions of entries where mother id's are missing (\code{NA}) or blank}

\item{sexcode.error}{A vector of integers containing the row positions of entries with mis-specified sex}
\item{zero.individ}{A vector of integers containing the row positions with individ (individual id) equal to 0}
\item{both.mother.father}{A list where the names of the list are family id's and the entries of the list are vectors containing the id's appearing as both mothers and fathers.}

Returns \code{NULL} if no errors were found.
}


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\author{Cecelia Laurie
}

\seealso{\code{\link{pedigreeCheck}}, \code{\link{pedigreeFindDuplicates}}, \code{\link{pedigreePairwiseRelatedness}}
}


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\examples{
family <- c("a","a","a","b","b","c","")
individ <- c("A","B","C","A","B",0,"")
mother <- c("B","C",0,0,0,NA,0)
father <- c("C","D",0,0,"",0,"D")
sex <- c("F","2","M","F","F","M","F")
samp <- data.frame(family, individ, mother, father, sex, stringsAsFactors=FALSE)
pedigreeClean(samp)

# $family.missing
# [1] 7
#
# $individ.missing
# [1] 7
#
# $father.missing
# [1] 5
#
# $mother.missing
# [1] 6
#
# $sexcode.error
# [1] 2
#
# $zero.individ
# [1] 6
#
# $both.mother.father
# $both.mother.father$a
# [1] "C"
}


\keyword{manip}

