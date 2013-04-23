### INITIAL CLEANING: Input is dataframe with columns labeled "family", "individ",
#"mother", "father", "sex"

# ASSUME "mother", "father" are 0 for founders
#Output: checks for NA or blank entries and outputs vector of row positions where NAs exist
#check sex is coded "M" and "F" and outputs row numbers with problems, and
#checks individual id's are not 0 and outputs row numbers with individual id = 0
#checks that no individual appears in both mother and father column 
# all output is rows of pedigree where error occurs except 'mother.father.same' gives family id

pedigreeClean<-function(pedigree)
{
  .Defunct("pedigreeCheck")
}


