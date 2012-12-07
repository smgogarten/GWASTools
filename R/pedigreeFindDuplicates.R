##This function identifies duplicates of individual within family
# and checks for matching pedigree information
#Input: dataframe with columns labeled "family", "individ","mother", "father", "sex"
#Output: two dataframes containing family id, individual id, number of copies
#dups$mismatch contains the info for duplicates that have mismatches in pedigree info
#dups$match contains the info for duplicates that match in pedigree info
# after correcting problems, the final matching duplicate dataframe will be needed later


pedigreeFindDuplicates<-function(pedigree, verbose=TRUE) 
{
  .Deprecated("pedigreeCheck")
  
	er<-pedigreeClean(pedigree)
      if(!is.null(er)) stop("basic errors - run pedigreeClean to identify")

	###Check duplicates within family
	dup<-data.frame()
	u<-unique(pedigree$family)
	un<-length(u)
	for(i in 1:un) {
	  fsamp<-pedigree[is.element(pedigree$family,u[i]),c("individ","mother","father","sex")] #get info for family i
	  ui<-unique(fsamp$individ)
	  uin<-length(ui)
	  if(uin <1) stop("error")
	   for( j in 1:uin) {
		 x<-fsamp[is.element(fsamp$individ,ui[j]), c("mother","father","sex")] #get rows for individual j in family i
		 if(dim(x)[1]<=1) next
		 chk<-all(length(unique(x$mother))==1 && length(unique(x$father))==1 && length(unique(x$sex))==1)
		 if(chk==FALSE) dup<-rbind(dup,c(u[i],ui[j],dim(x)[1],1)) else dup<-rbind(dup,c(u[i],ui[j],dim(x)[1],0))
		  
		}
	 }
	if(dim(dup)[1]==0) {
	  if (verbose)
		  message("There are no duplicates")
	  out.list<- NULL
	} 
	else 
	{
	names(dup)<-c("family","individ","copies","match")
	dup$match[is.element(dup$match,0)]<-"Yes"
	dup$match[is.element(dup$match,1)]<-"No"
	
	dups.match<-dup[is.element(dup$match,"Yes"),c("family","individ","copies")]
	dups.mismatch<-dup[is.element(dup$match,"No"),c("family","individ","copies")]
	
	out.list<-list(dups.mismatch,dups.match); names(out.list)<-c("dups.mismatch", "dups.match")
	}
	return(out.list)
}
