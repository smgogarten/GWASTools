\name{pedigreeClean}
\alias{pedigreeClean}

\title{Basic pedigree data checking
}

\description{This function checks a pedigree for completeness and gross errors
}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\usage{
pedigreeClean(pedigree, verbose = TRUE)
}

\arguments{

\item{pedigree}{A dataframe containing the pedigree information for the samples to be examined with columns labeled "family", "individ", "mother", "father" and "sex" containing the identity numbers of the family, individual, individual's mother, individual's father and individual's gender (coded as "M" or "F") .}

\item{verbose}{Logical value specifying whether or not to show progress information.}	

}


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\details{

The function performs a basic check on pedigree data for gross errors, checking for missing id's, non-integer id's, mis-coded gender, id's of 0 and for individuals that appear as both mothers and fathers. 


%Output: checks for NA entries and outputs vector of row positions where NAs exist,
%distinquishing NA entries in family id from NA entries where family id ok
%, checks all columns other than sex are numeric - outputs names of columns with problems, 
%check sex is coded "M" and "F" and outputs row numbers with problems, and
%checks individual id's are not 0 and outputs row numbers with individual id = 0
%checks that no individual appears in both mother and father column 


}


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\value{
A list with the following components:
\item{fam.na}{A vector of integers containing the row positions of entries with missing family id's}
\item{other.na}{A vector of integers containing the row positions of entries with missing individual, mother or father id's}
\item{cols.not.numeric}{A vector of integers containing the row positions of non-numeric id's}
\item{rows.sexcode.error}{A vector of integers containing the row positions of entries with mis-specified gender}
\item{zero.individ}{A vector of integers containing the row positions with an id equal to 0.}
\item{mofa}{A vector containing the id's of individuals appearing as both mothers and fathers}

Returns \code{NULL} if no errors were found.
}


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\author{Cecelia Laurie
}

\seealso{\code{\link{pedigreeCheck}}, \code{\link{pedigreeFindDuplicates}}, \code{\link{pedigreePairwiseRelatedness}}
}


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\examples{
family <- c(1,1,1,NA,2,2,2,2)
individ <- c(1,2,3,0,4,5,6,NA)
mother <- c(0,0,1,1,0,0,4,4)
father <- c(0,0,2,2,0,0,5,4)
sex <- c("F","M","F","F","F","F","M",1)
pedigree <- data.frame(family, individ, mother, father, sex)
pedigreeClean(pedigree)

#$fam.na
#[1] 4

#$other.na
#[1] 8

#$cols.not.numeric
#NULL

#$rows.sexcode.error
#[1] 8

#$zero.individ
#[1] 4

#$mofa
#[1] 2
}


\keyword{manip}

