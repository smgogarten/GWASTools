\name{pedigreeCheck}
\alias{pedigreeCheck}

\title{Testing for internal consistency of pedigrees}

\description{
Find inconsistencies within pedigrees.
}


\usage{
pedigreeCheck(pedigree)
}



\arguments{

   \item{pedigree}{A dataframe containing the pedigree information for the samples to be examined with columns labeled "family", "individ", "mother", "father" and "sex" containing the identity numbers of the family, individual, individual's mother, individual's father and individual's gender (coded as "M" or "F") .}
   
  
}






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




\details{

The function \code{pedigreeCheck} finds any of a number of possible inconsitencies within pedigree data: single individual families, gender mismatches with mother or father, impossible relationships (either where a child is a parent of self or an individual is both a child and a parent of the same person) and families that consist of two unrelated persons. It also looks for multiple subfamilies within a family.



}
\value{
The output for \code{pedigreeCheck} is a list with the following elements:

\item{one.person}{A vector of family ids for one-person families}
\item{mismatch.sex }{A vector of of family ids where sex of mother and/or father is incorrect }
\item{impossible.related }{A vector of of family ids where either child is mother of self or an individual is both child and mother of same person }
\item{duos }{A vector of of family ids where 'family' consists of only 2 unrelated persons }
\item{subfamilies.ident }{A matrix with columns for the family id (of
  'families' with multiple subfamilies), subfamily identifier and
  individual of each person in the identified subfamilies. Note that
  subfamilies are not identified for any families already identified
  with problems, and that the individual id's include individuals
  identified as a mother or father who may not be listed as individuals
  in the pedigree.
  
 If no inconsistencies are found, the output is \code{NULL}.
}



}
\author{Cecelia Laurie
}
\note{Subfamilies are not identified for any families already identified with problems, and individual id's in subfamilies may include individual ids listed for a mother or father who may not be listed as individuals in the pedigree.
}


\seealso{\code{\link{pedigreeClean}}, \code{\link{pedigreeFindDuplicates}}, \code{\link{pedigreePairwiseRelatedness}}
}

\examples{
family <- c(1,1,1,2,2,2,3)
individ <- c(1,2,3,4,5,6,7)
mother <- c(0,0,1,0,0,4,0)
father <- c(0,0,2,0,0,5,0)
sex <- c("F","M","F","F","F","M","M")
pedigree <- data.frame(family, individ, mother, father, sex)
pedigreeCheck(pedigree)
#$one.person
#[1] 3

#$mismatch.sex
#[1] 2

#$impossible.related
#NULL

#$subfamilies.ident
#data frame with 0 columns and 0 rows
}

\keyword{manip}

