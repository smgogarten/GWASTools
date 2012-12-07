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

   \item{pedigree}{A dataframe containing the pedigree information for the samples to be examined with columns labeled "family", "individ", "mother", "father" and "sex" containing the identifiers of the family, individual, individual's mother, individual's father and individual's sex (coded as "M" or "F") .
      Identifiers can be integer, numeric or character but identifiers for mother and father for founders are assumed to be 0.  }    
}




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




\details{

The function \code{pedigreeCheck} finds any of a number of possible errors and inconsistencies within pedigree data.
 If no problems are encountered, the output is \code{NULL}. If problems are encountered, output contains information
  for the errors encountered (a sub-list of the output values described below) and the following message is printed:
"All row numbers refer to rows in the full pedigree (not just within a family). Correct current problems and rerun pedigreeCheck.
There may be additional problems not investigated because of the current problems."

}
\value{
The output for \code{pedigreeCheck} is \code{NULL} or a sub-list of the following: 
\item{family.missing.rows}{A vector of integers containing the row positions of entries in the full pedigree where family id's are missing (NA) or blank}
\item{individ.missing_or_0.rows}{A vector of integers containing the row positions of entries in the full pedigree where individual id's are missing (NA), blank, or 0}
\item{father.missing.rows}{A vector of integers containing the row positions of entries in the full pedigree where father id's are missing (NA) or blank}
\item{mother.missing.rows}{A vector of integers containing the row positions of entries in the full pedigree where mother id's are missing (NA) or blank}
\item{sexcode.error.rows}{A vector of integers containing the row positions of entries in the full pedigree where the 'sex' variable is mis-coded}

\item{both.mother.father}{A data.frame with the variables 'family','parentID','mother.row',and 'father.row' where
    'family' = family identifier, 'parentID' = identifier of parent that appears as both mother and father,
     'father.row' = row positions(s) in full pedigree in which parent appears as father, and 'mother.row' = row position(s) in full pedigree in which parent appears as mother
     (if mutliple rows, row numbers are concatenated with separator = ';')}
\item{parent.no.individ.entry}{A data.frame with the variables 'row.num', 'family', 'no_individ_entry', and 'parentID', 
  where 'row.num' = row position of entry in the full pedigree where mother and/or father IDs are not included in the pedigree,
'family' = family identifier,  'no_individ_entry' has values 'father', 'mother' or 'both' indicating which parent is not in the pedigree, and
    'parentID' = the identifier(s) for individuals not in the pedigree (if more than one, identifiers are concatenated with separator =';')}
\item{unknown.parent.rows}{A data.frame with variables 'row.num' = row position in full pedigree where one parent is known and one parent is unknown and 'family' = family identifier.}
\item{duplicates}{A data.frame with variables 'family' = family identifier, 'individ' = individual identifier, 'copies' = number of copies of individual and  'match'= T/F depending upon whether all copies have identical pedigree information} 
\item{one.person.fams}{A data.frame identifying singeltons (one person families) with variables 'family' = family identifier and 'founder' = T/F depending up whether the singleton is a founder or not }

\item{mismatch.sex }{A data.frame with variables 'family' = family identifier and 'individ' = individual identifier for individuals that occur as mothers but sex is "M" or occur as fathers but sex is "F"  }
\item{impossible.related.rows }{A list where each entry in the list contains a set of row positions in the full pedigree which together indicate impossible relationships: where either a child is mother of self or an individual is both child and mother of the same person. Names of list entries are associated family identifiers. }

\item{subfamilies.ident }{A data.frame with variables 'family' = family identifier, "subfamily" = sub-family identifier within family, and
   'individ' = individual identifier of members of identified sub-family. }
  
 If no inconsistencies are found, the output is \code{NULL}.

}


\author{Cecelia Laurie
}
\note{All row numbers in output refer to row positions in the full pedigree (not just within family). User should correct current problems and rerun pedigreeCheck.
There may be additional problems not investigated because of the current problems.
}


\seealso{\code{\link{pedigreePairwiseRelatedness}}
}

\examples{

#basic errors
family <- c("a","a","a","b","b","c","")
individ <- c("A","B","C","A","B",0,"")
mother <- c("B","C",0,0,0,NA,0)
father <- c("C","D",0,0,"",0,"D")
sex <- c("F","2","M","F","F","M","F")
samp <- data.frame(family, individ, mother,father,sex,stringsAsFactors=FALSE)
pedigreeCheck(samp)
# there are other problems not investigated since 
#    the above are basic problems to be cleared up first

## 'duplicates', 'both.mother.father', 'parent.no.individ.entry'
family <- c("b","b","b","b","c","c",rep("d",5))
individ <- c("A","B","C","A","B","B",1:5)
mother <- c("B",0,0,"D",0,0,0,0,1,2,1)
father <- c("C",0,0,"C",0,0,0,0,2,1,2)
sex <- c("F","F","M","M","F","F","F","M","F","F","M")
samp <- data.frame(family, individ, mother,father,sex,stringsAsFactors=FALSE)
pedigreeCheck(samp)
# there are other problems (such as mismatch.sex) but not investigated 
#     directly because already had both.mother.father inconsistency

# 'parent.no.individ.entry', 'one.person.fams', 'unknown.parent.rows',
#    'mismatch.sex','impossible.related.rows'
family <- c(1,1,1,2,2,2,3,4,4,4,5,5,5,5,6,6,6)
individ <- c(1,2,3,1,2,3,1,1,3,2,1,2,3,4,1,2,3)
mother <- c(2,0,1,2,1,0,1,2,0,2,2,4,0,0,2,1,0)
father <- c(3,0,3,0,3,0,2,3,1,0,3,1,0,0,3,3,0)
sex <- c("F","F","M","F","F","M","F","F","F","F","M","F","M","F","F","M","F")
samp <- data.frame(family, individ,mother,father,sex,stringsAsFactors=FALSE)
pedigreeCheck(samp)
# 'mismatch.sex' and 'impossible.related.rows' are only investigated 
#      for families where there are no other inconsistencies

## 'subfamilies.ident'
family <- rep(1,12)
individ <- 1:12
mother <- c(0,0,2,2,0,0,5,0,7,0,0,10)
father <- c(0,0,1,1,0,0,6,0,8,0,0,11)
sex <- c("M",rep("F",4),"M","F","M","M","F","M","M")
samp <- data.frame(family,individ,mother,father,sex,stringsAsFactors=FALSE)
pedigreeCheck(samp)
# 'subfamilies.ident' is only investigated for families 
#     where there are no other inconsistencies
}

\keyword{manip}

