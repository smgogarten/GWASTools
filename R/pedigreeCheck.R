###
##Testing for internal consistency of families
###

#Requires access to the function "relate.r"

#Input: dataframe with family,individ,mother,father,sex
#Output is a list: 
# $one.person is vector of family ids for one-person families
# $mismatch.sex consists of family ids where sex of mother and/or father is incorrect
# $impossible.related consists of family ids where either child is mother of self or an individual is both child and mother of same person
# $duos consists of family ids where 'family' consists of only 2 unrelated persons
# $subfamilies.ident is a matrix with family id (of 'families' with multiple subfamilies), subfamily identifier, individual ids of persons in the subfamily
# (Note: subfamilies are not identified for any families already identified with problems)
#(NOTE: In subfamilies.ident the individual id's include individuals identified as mother or father who may not be in individ)

pedigreeCheck<-function(pedigree)
{
		
	## local function for testing for relatedness
	relate<-function(ind,x)
	{
		clust<-ind
		if(!is.element(ind, x$individ))
			stop("Error - there is no individual with id",ind,"in family",unique(x$family),"\n");
		
		mo<-x$mother[is.element(x$individ, ind)]
		fa<-x$father[is.element(x$individ, ind)]
		if(mo!=0) clust<-c(clust,mo)
		if(fa!=0) clust<-c(clust,fa)
		clust<-c(clust,x$individ[is.element(x$mother,ind)|is.element(x$father,ind)])
		return(clust)
	}
		
		
		
		
		
	#Initial checks to verify initial cleaning 
	#(no NAs, sex code ok, no 0 individual ids, all numeric (except sex), no duplicates) 
	
	## Checking entries not NA
	entries.na<-which(is.na(pedigree),arr.ind=TRUE)
	if(length(entries.na)!=0) {stop("Error - There are NAs. Run pedigreeClean to locate"); }
	
	##Checking columns that not numeric 
	
	cols.not.numeric<-NULL
	tsamp<-pedigree[,c("family","individ","mother","father")]
	n<-dim(tsamp)[2]
	for (i in 1:n) 
	{
		if (!is.numeric(tsamp[,i])) {stop("Error. Some columns are not numeric. Run pedigreeClean to locate");}
	}
	##Checking sex code
	rows.sexcode.error<-which(!is.element(pedigree$sex,c("M","F")))
	if(length(rows.sexcode.error)!=0) {stop("Error in sex code. Run pedigreeClean to locate");}
	
	##Checking individual id's of 0
	zero.individ<-which(pedigree$individ==0)
	if(length(zero.individ)!=0) {stop("Error - There are individual ids which are 0. Run pedigreeClean to locate");}
	
	##Checking duplicates
	u<-unique(pedigree$family)
	un<-length(u)
	for(i in 1:un) 
	{
		 fsamp<-pedigree[is.element(pedigree$family,u[i]),c("individ","mother","father","sex")] #get info for family i
		 ui<-unique(fsamp$individ)
		 uin<-length(ui)
		 if(uin <1) message("error")
		 for( j in 1:uin) 
		 {
			 x<-fsamp[is.element(fsamp$individ,ui[j]), c("mother","father","sex")] #get rows for individual j in family i
			 if(dim(x)[1]<=1) next else { stop("Error - there are duplicates. Run pedigreeFindDuplicates to locate and pedigreeDeleteDuplicates to remove")}
		 } 
	}
	
	
	## Finding one person families and doing parent sex check (one,bsex)
	u<-unique(pedigree$family)
	un<-length(u)
	one<-NULL
	bsex<-NULL
	badfam<-NULL
	for (i in 1:un) 
	{
		fsamp<-pedigree[is.element(pedigree$family,u[i]),c("individ","mother","father","sex")] #get info for family i
		if(dim(fsamp)[1]==1) {one<-c(one,u[i]);next}
		moth<-fsamp$mother[fsamp$mother !=0]
		ind<-fsamp[is.element(fsamp$individ,moth),]
		if(!all(ind$sex=="F")) {bsex<-c(bsex,u[i]);next}
		fath<-fsamp$father[fsamp$father !=0]
		ind<-fsamp[is.element(fsamp$individ,fath),]
		if(!all(ind$sex=="M")) {bsex<-c(bsex,u[i]);next}
		
		##Testing for impossible relationships
		#child is parent of self or individual is child and parent of same person
		
		ui<-fsamp$individ
		uin<-length(ui)
	
		moin<-NULL
		fain<-NULL
		for( j in 1:uin) 
		{
			moin<-c(moin,paste(fsamp$mother[j],fsamp$individ[j]))
			fain<-c(fain,paste(fsamp$father[j],fsamp$individ[j]))
		}
		for (j in 1:uin)
		{
			inmo<-paste(fsamp$individ[j],fsamp$mother[j])
			infa<-paste(fsamp$individ[j],fsamp$father[j])
			if(is.element(inmo,union(moin,fain)) | is.element(infa,union(moin,fain))) {badfam<-c(badfam,u[i]);break}
		}
	}#end of family loop
	
	os<-length(one)==0 && length(bsex)==0 && length(badfam)==0
	
	bf<-c(one,bsex,badfam)
	ug<-setdiff(u,bf)
	n<-length(ug)
	
	if(n==0)
	{
		out.list<-list(one,bsex,badfam)
		names(out.list)<-c("one.person","mismatch.sex","impossible.related")
	} 
	else 
	{
	
		subfams<-data.frame()
		sfam<-NULL  #record of families that have subfamilies
	
		for (i in 1:n) 
		{
			fam<-pedigree[is.element(pedigree$family,ug[i]),c("individ","mother","father")]
		
		#add 'extra' mothers and fathers to individ list
		ui<-unique(fam$individ)
		um<-fam$mother
		uf<-fam$father
		um<-setdiff(um,0); uf<-setdiff(uf,0)
		em<-setdiff(um,intersect(ui,um))
		ef<-setdiff(uf,intersect(ui,uf))
		new<-c(em,ef);ln<-length(new)
		if(ln !=0) {
		individ<-c(fam$individ,new)
		mother<-c(fam$mother,rep(0,ln))
		father<-c(fam$father,rep(0,ln))
		fam<-data.frame(individ,mother,father)} 
		
		#check relatedness
		  flag1<-0 ; counter<-0; Klus<-NULL
		  while(flag1==0) 
		  {
				ui<-setdiff(unique(fam$individ),Klus)
				clus<-NULL
				uclus<-NULL
				choice<-ui
				flag2<-0; counter<-counter+1
				  while(flag2==0) {
				   ind<-choice[1]
				   nclus<-relate(ind,fam)  #find everyone related to ind, includes ind
				   clus<-union(clus,nclus)  #all related so far
				   uclus<-union(uclus,ind)   #individuals already 'used' in relatedness check
				   choice<-setdiff(clus,uclus)
				   if(length(choice)==0) flag2<-1
				   if(all(is.element(ui,clus))) flag2<-1
				   } #end: finds a subfamily
				 Klus<-union(Klus,clus)
				 if(all(is.element(fam$individ,Klus))) 
					 	flag1<-1  #everyone accounted for
				 if(flag1==1 && counter==1) 
					 	break
				 sfam<-union(sfam,ug[i])
				 f<-rep(ug[i],length(clus))
				 sf<-rep(counter,length(clus))
				 newsf<-cbind(f,sf,clus)
				 subfams<-rbind(subfams,newsf)
		  } 
		} 
		if(length(sfam)!=0) { names(subfams)<-c("family","subfamily","individ") }
	
	
		len<-c(length(one),length(bsex),length(badfam))
		if(all(len==0)&& length(sfam)==0)
			out.list<- NULL
		else
		{
	
			out.list<-list(one,bsex,badfam,subfams)
			names(out.list)<-c("one.person","mismatch.sex","impossible.related","subfamilies.ident")
		}
	}
	return(out.list) 
}
