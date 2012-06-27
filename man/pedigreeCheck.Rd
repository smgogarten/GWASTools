\name{pedigreeCheck}
\alias{pedigreeCheck}

\title{Testing for internal consistency of pedigrees}

\description{
Find inconsistencies within pedigrees.
}


\usage{
pedigreeCheck(pedigree, unknown.parents=FALSE)
}



\arguments{

   \item{pedigree}{A dataframe containing the pedigree information for the samples to be examined with columns labeled "family", "individ", "mother", "father" and "sex" containing the identifiers of the family, individual, individual's mother, individual's father and individual's sex (coded as "M" or "F") .
      Identifiers can be integer, numeric or character but identifiers for mother and father for founders are assumed to be 0. Identifier for unknown parent is assumed to be 0.
    }
    \item{unknown.parents}{A logical indicating whether to check for
      unknown parents.}
   
  
}






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




\details{

The function \code{pedigreeCheck} finds any of a number of possible inconsistencies within pedigree data: single individual families, sex mismatches with mother or father, and impossible relationships (either where a child is a parent of self or an individual is both a child and a parent of the same person). 
   It also looks for multiple subfamilies within a family. If the option \code{unknown.parents=TRUE}, it will identify (and output) individuals with one known parent and one unknown parent.
  The function first checks the pedigree for gross errors.
}
\value{
The output for \code{pedigreeCheck} is a list with the following elements if the option \code{unknown.parents=FALSE} (the default):

\item{one.person}{A vector of family ids for one-person families}
\item{mismatch.sex }{A vector of family ids where sex of mother and/or father is incorrect }
\item{impossible.related }{A vector of family ids where either child is mother of self or an individual is both child and mother of same person }
\item{subfamilies.ident }{A matrix with columns for the family id (of
  'families' with multiple subfamilies), subfamily identifier and
  individual of each person in the identified subfamilies. Note that
  subfamilies are not identified for any families already identified
  with problems, and that the individual id's include individuals
  identified as a mother or father who may not be listed as individuals
  in the pedigree.}
  
 If no inconsistencies are found, the output is \code{NULL}.

 If \code{unknown.parents = TRUE} the output contains as well:
\item{unknown.parents}{ a sub-data.frame of pedigree identifying individuals with one parent known and one parent unknown.}
}


\author{Cecelia Laurie
}
\note{Subfamilies are not identified for any families already identified with problems, and individual id's in subfamilies may include individual ids listed for a mother or father who may not be listed as individuals in the pedigree.
}


\seealso{\code{\link{pedigreeClean}}, \code{\link{pedigreeFindDuplicates}}, \code{\link{pedigreePairwiseRelatedness}}
}

\examples{
family <- c(1,1,1,2,2,2,3,4,4,4,4,4)
individ <- c(1,2,3,1,2,3,1,1,2,3,4,5)
mother <- c(2,0,1,2,1,0,1,2,0,0,5,0)
father <- c(3,0,3,0,3,0,2,3,0,0,0,0)
sex <- c("F","F","M","F","M","M","F","M","F","M","F","F")
samp <- data.frame(family, individ, mother, father, sex, stringsAsFactors=FALSE)
pedigreeCheck(samp, unknown.parents=TRUE)

# $one.person
# [1] 3
#
# $mismatch.sex
# [1] 2
#
# $impossible.related
# [1] 1 2 3
#
# $subfamilies.ident
#   family subfamily individ
# 1      4         1       1
# 2      4         1       2
# 3      4         1       3
# 4      4         2       4
# 5      4         2       5
#
# $unknown.parents
#    family individ mother father sex
# 4       2       1      2      0   F
# 11      4       4      5      0   F
}

\keyword{manip}

