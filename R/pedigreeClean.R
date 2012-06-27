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
	cl<-unlist(lapply(pedigree,class))
      if(is.element("factor",cl)) stop("some variables are factors")	

	## Finding entries that are missing
      wfam<-which(is.na(pedigree$family) | is.element(pedigree$family,""))
      wind<-which(is.na(pedigree$individ) | is.element(pedigree$individ,""))
      wfa<-which(is.na(pedigree$father) | is.element(pedigree$father,""))
      wmo<-which(is.na(pedigree$mother) | is.element(pedigree$mother,""))
	
	
	##Checking sex code
	rows.sexcode.error<-which(!is.element(pedigree$sex,c("M","F")))
	
	##Checking individual id's of 0
	zero.individ<-which(pedigree$individ %in% 0)
	
	##Checking that same person is not both mother and father
	mofa<-NULL
	u<-unique(pedigree$family[!is.na(pedigree$family) & !is.element(pedigree$family,"")])
      mf<-list()
      cnt<-1
	for (i in u){tsamp<-pedigree[is.element(pedigree$family,i),]
	  mo<-unique(tsamp$mother[!is.element(tsamp$mother,0)&!is.na(tsamp$mother) & !is.element(tsamp$mother, "")])
	  fa<-unique(tsamp$father[!is.element(tsamp$father,0)&!is.na(tsamp$father) & !is.element(tsamp$father,"")])
	  prob<-intersect(mo,fa)
	  if(length(prob)!=0){mofa<-c(mofa,i); mf[[cnt]]<-prob;cnt<-cnt+1}
      }
      if(length(mofa)!=0) names(mf)<-mofa

 
      out<-list(wfam,wind,wfa,wmo,rows.sexcode.error,zero.individ)
      lens<-unlist(lapply(out,length))
      mfln<-length(mf)
      if(all(lens==0) & mfln==0){
		out.list <- NULL;return(out.list)} else {
		out.list<-list(wfam,wind,wfa,wmo,rows.sexcode.error,zero.individ,mf)
		names(out.list)<-c("family.missing","individ.missing","father.missing","mother.missing","sexcode.error","zero.individ","both.mother.father");
		return(out.list)
      }
}


