### INITIAL CLEANING: Input is dataframe with columns labeled "family", "individ",
#"mother", "father", "sex"

#Output: checks for NA entries and outputs vector of row positions where NAs exist,
#distinquishing NA entries in family id from NA entries where family id ok
#, checks all columns other than sex are numeric - outputs names of columns with problems, 
#check sex is coded "M" and "F" and outputs row numbers with problems, and
#checks individual id's are not 0 and outputs row numbers with individual id = 0
#checks that no individual appears in both mother and father column 

pedigreeClean<-function(pedigree, verbose = TRUE)
{
	## Finding entries that are NA
	ent<-which(is.na(pedigree),arr.ind=TRUE)
	entries.na<-unique(ent[ ,1])
	fam.na<-which(is.na(pedigree$family))
	other.na<-setdiff(entries.na,fam.na)
	
	
	##Finding columns that are not numeric - providing the names
	cols.not.numeric<-NULL
	tsamp<-pedigree[,c("family","individ","mother","father")]
	n<-dim(tsamp)[2]
	for (i in 1:n) {
	 if (!is.numeric(tsamp[,i])) cols.not.numeric<-c(cols.not.numeric,names(tsamp)[i])} 
	
	##Checking sex code
	rows.sexcode.error<-which(!is.element(pedigree$sex,c("M","F")))
	
	##Checking individual id's of 0
	zero.individ<-which(pedigree$individ==0)
	
	##Checking that same person is not both mother and father
	mofa<-NULL
	u<-unique(pedigree$family)
	for (i in u){tsamp<-pedigree[is.element(pedigree$family,i),]
	mo<-unique(tsamp$mother[!is.element(tsamp$mother,0)&!is.na(tsamp$mother)])
	fa<-unique(tsamp$father[!is.element(tsamp$father,0)&!is.na(tsamp$father)])
	prob<-intersect(mo,fa)
	if(length(prob)!=0)mofa<-c(mofa,i) }
	
	if(length(entries.na)==0 && length(cols.not.numeric)==0 && length(rows.sexcode.error)==0 && length(zero.individ)==0 && length(mofa)==0) 
	{   
		out.list <- NULL
	} 
	else 
	{
		out.list<-list(fam.na,other.na,cols.not.numeric,rows.sexcode.error,zero.individ,mofa)
		names(out.list)<-c("fam.na","other.na","cols.not.numeric","rows.sexcode.error","zero.individ","mofa");
	}
	return(out.list)
}


