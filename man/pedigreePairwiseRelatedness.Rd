\name{pedigreePairwiseRelatedness}
\alias{pedigreePairwiseRelatedness}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\title{
Calculate theoretical pairwise relatedness values from pedigrees
}
\description{
This function calculates the pairwise relatedness values from pedigree data.
}



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\usage{
pedigreePairwiseRelatedness(pedigree, use.any.ids = FALSE)
}



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\arguments{
  \item{pedigree}{A dataframe containing the pedigree information for the samples to be examined with columns labeled "family", "individ", "mother", "father" and "sex" containing the identity numbers of the family, individual, individual's mother, individual's father and individual's gender (coded as "M" or "F") .}

  \item{use.any.ids}{A logical value specifying whether pairs of individuals should be created using only id's listed in the "individ" column (if 
  \code{FALSE} or if pairs should be created using any id's contained in  "individ", "mother" or "father" columns.}
}


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\details{  The function assumes (and checks) that there are no one person families, no mismatched mother/father sexes and no impossible relationships. Relatedness is not calculated for inbred families.

}



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\value{

A list with the following components:

\item{inbred.fam}{A vector of id's of families with inbreeding (to be handled by hand) }
\item{relativeprs}{A dataframe with columns "Individ1", "Individ2", "relation",  "kinship coefficient" and "family" containing the id's of the pair of individuals, the relationship between the individuals if closely related (possible values are "UN" = unrelated, "PO" = parent/offspring, "FS" = full siblings, "HS" = half siblings, and "FC" = first cousins), kinship coefficient and family id. }

}


\author{Cecelia Laurie

}




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\seealso{\code{\link{pedigreeClean}}, \code{\link{pedigreeCheck}}, \code{\link{pedigreeFindDuplicates}}

}


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\examples{
family <- c(1,1,1,1,2,2,2,2)
individ <- c(1,2,3,4,5,6,7,8)
mother <- c(0,0,1,1,0,0,5,5)
father <- c(0,0,2,2,0,0,6,0)
sex <- c("F","M","F","F","F","M","M","M")
pedigree <- data.frame(family, individ, mother, father, sex)
pedigreePairwiseRelatedness(pedigree)

#$inbred.fam
#NULL

#$relativeprs
#   Individ1 Individ2 relation kinship family
#1         1        2        U   0.000      1
#2         1        3       PO   0.250      1
#3         1        4       PO   0.250      1
#4         2        3       PO   0.250      1
#5         2        4       PO   0.250      1
#6         3        4       FS   0.250      1
#11        5        6        U   0.000      2
#21        5        7       PO   0.250      2
#31        5        8       PO   0.250      2
#51        6        7       PO   0.250      2
#61        6        8        U   0.000      2
#8         7        8       HS   0.125      2
}

\keyword{manip }


