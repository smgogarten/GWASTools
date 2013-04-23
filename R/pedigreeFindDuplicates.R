##This function identifies duplicates of individual within family
# and checks for matching pedigree information
#Input: dataframe with columns labeled "family", "individ","mother", "father", "sex"
#Output: two dataframes containing family id, individual id, number of copies
#dups$mismatch contains the info for duplicates that have mismatches in pedigree info
#dups$match contains the info for duplicates that match in pedigree info
# after correcting problems, the final matching duplicate dataframe will be needed later


pedigreeFindDuplicates<-function(pedigree, verbose=TRUE) 
{
  .Defunct("pedigreeCheck")
}
