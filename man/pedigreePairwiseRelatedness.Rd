\name{pedigreePairwiseRelatedness}
\alias{pedigreePairwiseRelatedness}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\title{
Assign relatedness from pedigree data
}
\description{
This function assigns relationships from pedigree data. Output includes the theoretical pairwise kinship coefficients.
}



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\usage{
pedigreePairwiseRelatedness(pedigree, use.any.ids=FALSE)
}



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\arguments{
  \item{pedigree}{A dataframe containing the pedigree information for the samples to be examined with columns labeled "family", "individ", "mother", "father" and "sex" containing the identifiers for family, individual, individual's mother, individual's father and individual's sex (coded as "M" or "F") . 
     Identifiers can be integer, numeric or character but identifiers for mother and father for founders are assumed to be 0. Identifier for unknown parent is assumed to be 0. Error messages are returned for \code{NA} or blank identifiers.}

  \item{use.any.ids}{A logical value specifying whether pairs of individuals should be created using only id's listed in the "individ" column (if 
  \code{FALSE}) or if pairs should be created using any id's contained in  "individ", "mother" or "father" columns.}
}


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\details{  Assigns relationships between individuals in a pedigree, including "U" = unrelated, "PO" = parent/offspring, "FS" = full siblings, "HS" = half siblings, "Av" = avuncular, "GpGc" = grandparent-grandchild, and "FC" = first cousins, among others).

Relatedness is not calculated for inbred families but kinship coefficients are.

}



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\value{

A list with the following components:

\item{inbred.fam}{A vector of id's of families with inbreeding
  (relationships are not assigned). }
\item{inbred.KC}{A dataframe for inbred families with columns
  "Individ1","Individ2", "kinship" and "family" containing the id's of
  the pair of individuals, kinship coefficient and family id. } 
\item{relativeprs}{A dataframe with columns "Individ1", "Individ2", "relation",  "kinship" and "family" containing the id's of the pair of individuals, the relationship between the individuals if closely related (possible values are "U" = unrelated, "PO" = parent/offspring, "FS" = full siblings, "HS" = half siblings, "Av" = avuncular, "GpGc" = grandparent-grandchild, and "FC" = first cousins, among others), kinship coefficient and family id. }

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
pedigree <- data.frame(family, individ, mother, father, sex, stringsAsFactors=FALSE)
pedigreePairwiseRelatedness(pedigree)

#$inbred.fam
#NULL
#
#$inbred.KC
#NULL
#
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

family <- rep(1,3)
individ <- c("I1","I3","I4")
mother <- c(0,"I1","I1")
father <- c(0,"I2","I2")
sex <- c("F","M","F")
samp1 <- data.frame(family, individ, mother, father, sex, stringsAsFactors=FALSE)
# note father 2 is known but 'not genotyped'

out <- pedigreePairwiseRelatedness(samp1)
out$relativeprs
#  Individ1 Individ2 relation kinship family
#1       I1       I3       PO    0.25      1
#2       I1       I4       PO    0.25      1
#4       I3       I4       FS    0.25      1

out2 <- pedigreePairwiseRelatedness(samp1, use.any.ids=TRUE)
out2$relativeprs
#  Individ1 Individ2 relation kinship family
#1       I1       I3       PO    0.25      1
#2       I1       I4       PO    0.25      1
#3       I1       I2        U    0.00      1
#4       I3       I4       FS    0.25      1
#5       I3       I2       PO    0.25      1
#6       I4       I2       PO    0.25      1

family <- rep(2,7)
individ <- paste("I",c(1,2,3,4,5,6,7),sep="")
mother <- c(0,0,0,"I1","I1","I3","I5")
father <- c(0,0,0,"I2","I2","I4","I4")
sex <- c("F","M","F","M","F","F","F")
samp2 <- data.frame(family, individ, mother, father, sex, stringsAsFactors=FALSE)
pedigreePairwiseRelatedness(samp2)
# inbred family
}

\keyword{manip }


