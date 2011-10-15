



##Deleting duplicates in pedigrees
#  input pedigree dataframe pedigree, input dataframe identifying duplicates
#  where duplicates is dataframe with duplicates$family (family id) and duplicates$individ (individual id)
#  (there could be other columns)
#  duplicates could be the output from ident.dups
  
  
  
pedigreeDeleteDuplicates <- function(pedigree, duplicates) 
{
	dr<-NULL
	for (i in 1:dim(duplicates)[1]) 
	{
		sel<- is.element(pedigree$family,   duplicates$family[i]) & 
		      is.element(pedigree$individ,  duplicates$individ[i])
		      
		rst<- row.names( pedigree[sel , ])
		dr<- c(dr,rst[2:length(rst)])  
	}
	
	out<-pedigree[ !is.element(row.names(pedigree),dr), ]
	
}
