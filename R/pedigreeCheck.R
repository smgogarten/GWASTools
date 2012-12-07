###Input is dataframe with columns labeled "family", "individ", "mother", "father", "sex"
# ASSUME "mother", "father" are 0 for founders

# any row numbers in output refer to row numbers from full pedigree

#Output: checks for NA or blank entries and outputs vector of row positions where NAs exist
#  for individ also checks entries of 0 - "family.missing","individ.missing_or_0","father.missing","mother.missing"
#check sex is coded "M" and "F" and outputs row numbers with problems - "sexcode.error"

# any family with any of the above errors are deleted from further analysis

#checks that no individual appears in both mother and father column 
# "both.mother.father" is data.frame of family and ID along with row(s) where occurs as mother and row(s)where occurs as father
#     if multiple rows, row numbers are concatenated using ";"

# "unknown.parent": checks for one known, one unknown parent - outputs row(s) of pedigree along with family id

#  "parent.no.individ.entry": checks for mother/father IDs not appearing in individ : output is dataframe with
#     family, row number in pedigree, indicator if mother or father or both are not in individ ("mother","father","both"),
#      and parent ID (if both, parent IDs are concatenated with ":")

# "one.person" is data.frame of one-person families: "family" = family id, "founder" = T/F is founder or other
# "dup" is data.frame of duplicates: "family","individ", "copies","match" where "copies" is number of copies
#       "match" is T/F if all pedigree information matches over all duplicates   

# any family with any of the above error/output are deleted from further analysis

###
##Testing for internal consistency of families
###

# "mismatch.sex": data.frame "family","individ","sex" "as.parent" where
 #     "family","individ","sex" is as in pedigree and "as.parent" is either 'mother' or 'father' depending upon whether individ appeared as mother or father 
# "impossible.related": list where each element of list is vector of row numbers (referring to full pedigree) involved in impossible relationships: 
 #       name of list element = family name
 #     either child is mother of self or an individual is both child and mother of same person

# any family with any of the above error/output are deleted from further analysis

#
# $subfamilies.ident is a matrix with family id (of 'families' with multiple subfamilies), subfamily identifier, individual ids of persons in the subfamily
############################

	## local function for testing for relatedness
	relate<-function(ind,x)
	{
		clust<-ind
		if(!is.element(ind, x$individ))
			stop("Error - there is no individual with id ",ind," in family ",unique(x$family),"\n");
		
		mo<-x$mother[is.element(x$individ, ind)]
		fa<-x$father[is.element(x$individ, ind)]
		if(mo!=0) clust<-c(clust,mo)
		if(fa!=0) clust<-c(clust,fa)
		clust<-c(clust,x$individ[is.element(x$mother,ind)|is.element(x$father,ind)])
		return(clust)
	}


pedigreeCheck<-function(pedigree)
{
      if(!is.element(class(pedigree),"data.frame")) stop("input must be data.frame")
      if(nrow(pedigree) ==0) stop("input data.frame has no data")
	cl<-unlist(lapply(pedigree,class))
      if(is.element("factor",cl)) stop("some variables are factors")
      req<-c("family","individ","mother","father","sex")
      ck<-all(is.element(req,names(pedigree)))
      if(!ck) stop(paste("data frame does not include all required variables: ",req))	

      mssg<-"All row numbers refer to row in the full pedigree (not just within family). \n Correct current problems and rerun pedigreeCheck.
There may be additional problems not investigated because of the current problems.\n"

	## Finding entries that are missing or otherwise unacceptable
      wfam<-which(is.na(pedigree$family) | is.element(pedigree$family,""))
      wind<-which(is.na(pedigree$individ) | is.element(pedigree$individ,"")| pedigree$individ %in% 0)
      wfa<-which(is.na(pedigree$father) | is.element(pedigree$father,""))
      wmo<-which(is.na(pedigree$mother) | is.element(pedigree$mother,""))
	
	
	##Checking sex code
	rows.sexcode.error<-which(!is.element(pedigree$sex,c("M","F")))

      ## bad families so far	
      bfam.ind<-union(wind,union(wfa,union(wmo,rows.sexcode.error)))
      if(length(bfam.ind) !=0){
            bfam.ind<-setdiff(bfam.ind,wfam)  
            bfam<-unique(pedigree$family[bfam.ind])
      } else { bfam<-NULL}

     out.listA<-list(wfam,wind,wfa,wmo,rows.sexcode.error)
     nms<-c("family.missing.rows","individ.missing_or_0.rows","father.missing.rows","mother.missing.rows","sexcode.error.rows")
     lensA<-unlist(lapply(out.listA,function(x) length(x))) 
     selt<-lensA>0
     if(sum(selt)!=0){outA<-out.listA[selt];names(outA)<-nms[selt]  } else{ outA<-NULL}   
      

##Checking that same person is not both mother and father
	mofa.all<-NULL
      ss<-NULL
      unkp<-NULL
      onef<-NULL #founders
      oneo<-NULL # other
      dup<-NULL
      u<-unique(pedigree$family[!is.na(pedigree$family) & !is.element(pedigree$family,"") & !is.element(pedigree$family,bfam)])

       if(length(u) ==0)  {message(mssg); return(outA)}

          
      for(fm in u){
        sel<-which(is.element(pedigree$family,fm))
        ped<-pedigree[sel,]
        mo<-!is.na(ped$mother) &!is.element(ped$mother,"")& !is.element(ped$mother,0) 
        fa<-!is.na(ped$father) &!is.element(ped$father,"") & !is.element(ped$father,0) 
        moid<-ped$mother[mo]
        faid<-ped$father[fa]
        prob<-intersect(moid,faid)
        if(length(prob)!=0){
          mofa<-data.frame("family"=rep(fm,length(prob)),"parentID"=prob,"mother.row"=rep(NA,length(prob)),"father.row"=rep(NA,length(prob)),stringsAsFactors=FALSE)
          for(k in 1:length(prob)){
            id<-prob[k]
            wm<-which(is.element(ped$mother,id))
            if(length(wm>1)) wmm<-paste(sel[wm],collapse=";")
            wf<-which(is.element(ped$father,id))
            if(length(wf>1)) wff<-paste(sel[wf],collapse=";")           
            mofa[k,]<-c(fm,id,wmm,wff)
          }
         mofa.all<-rbind(mofa.all,mofa)
        }

     # one unknown parent NOTE fa and mo indicate known parent IDs
       unk<-(!mo & fa) | (!fa & mo)
       if(sum(unk)!=0){
         unkr<-sel[unk]
         unkf<-rep(fm,length(unkr))
         unkk<-data.frame("row.num"=unkr, "family"=unkf,stringsAsFactors=FALSE)
         unkp<-rbind(unkp,unkk)
       }

       ## checking that there are individ entries for any mother or father ids
         momo<-mo & !is.element(ped$mother,ped$individ)
         fafa<- fa & !is.element(ped$father,ped$individ)
         bo<-momo & fafa
         mo2<-momo & !bo
         fa2<-fafa & !bo
         rn<-momo|fafa; ln<-sum(rn)
         if(ln!=0){
            s<-data.frame("row.num"=sel[rn],"family"=rep(fm,ln),"no_individ_entry" =rep(NA,ln),"parentID"=rep(NA,ln),stringsAsFactors=FALSE) 
            s$no_individ_entry[is.element(s$row.num,sel[mo2])]<-"mother"
            s$no_individ_entry[is.element(s$row.num,sel[fa2])]<-"father"  
            s$no_individ_entry[is.element(s$row.num,sel[bo])]<-"both"
            for(k in 1:nrow(s)){
              rw<-s$row.num[k]
              if(s$no_individ_entry[k] %in% "both"){
                 ids<-paste(pedigree$mother[rw],pedigree$father[rw],sep=";")
                 s$parentID[k]<-ids  } else {
                 s$parentID[k]<-pedigree[rw,s$no_individ_entry[k]]
              }
           }              
            ss<-rbind(ss,s)
         }
   
        # one person families
           if(nrow(ped)==1) {
               if(!(is.element(ped$mother,0) & is.element(ped$father,0))){
                  oneo<-c(oneo,fm)} else {onef<-c(onef,fm)}
           }
 
        # duplicates

          ui<-unique(ped$individ)
          uin<-length(ui)
          for(j in 1:uin){
             tmp<-ped[is.element(ped$individ,ui[j]),]
             if(nrow(tmp)<=1) next  
             chk<-all(length(unique(tmp$mother))==1 & length(unique(tmp$father))==1 & length(unique(tmp$sex))==1)
             dp<-data.frame("family"=fm,"individ"=ui[j],"copies"=nrow(tmp),stringsAsFactors=FALSE)
             if(chk==TRUE) dp$match<-TRUE else dp$match<-FALSE
             dup<-rbind(dup,dp)
          }
          
     } # end of family loop
 
     bfam2<-bfam
     if(!is.null(mofa.all))bfam2<-union(bfam2,unique(mofa.all$family))
     if(!is.null(dup))bfam2<-union(bfam2,unique(dup$family))
     if(!is.null(ss))bfam2<-union(bfam2,unique(ss$family))
     if(!is.null(onef)) bfam2<-union(bfam2,onef)
     if(!is.null(oneo)) bfam2<-union(bfam2,oneo)
     if(!is.null(unkp)) bfam2<-union(bfam2,unkp$family)

     oneall<-c(onef,oneo)
     one<-NULL
     if(!is.null(oneall)){
         oneL<-c(rep(TRUE,length(onef)),rep(FALSE,length(oneo)))
         one<-data.frame("family"=oneall,"founder"=oneL,stringsAsFactors=FALSE)
     }

     outlistB<-list(mofa.all,dup,ss,one,unkp)
     lensB<-unlist(lapply(outlistB,is.null))
     nms<-c("both.mother.father","duplicates","parent.no.individ.entry","one.person.fams","unknown.parent.rows")
     if(any(!lensB)) {outB<-outlistB[!lensB]; names(outB)<-nms[!lensB]} else {outB<-NULL}

     
      badsex<-NULL  

      u<-unique(pedigree$family[!is.na(pedigree$family) & !is.element(pedigree$family,"") & !is.element(pedigree$family,bfam2)])

      if(length(u) ==0) {message(mssg); return(c(outA,outB))}
      imposs<-vector("list",length(u))   
      for(i in 1:length(u)){ 
         fm<-u[i]
         sel<-which(is.element(pedigree$family,fm))
         ped<-pedigree[sel,]

#     parent sex check
         mo<-!is.element(ped$mother,0)
         indsel<-is.element(ped$individ,ped$mother[mo])
         ind<-ped[indsel,]
         wF<-!is.element(ind$sex,"F")
         if(sum(wF)!=0) {
           minds<-ind[wF,]
           tmp<-data.frame("family"=rep(fm,nrow(minds)),"individ"=minds$individ,stringsAsFactors=FALSE)
           badsex<-rbind(badsex,tmp)
         }
         fa<-!is.element(ped$father,0)
         indsel<-is.element(ped$individ,ped$father[fa])
         ind<-ped[indsel,]
         wM<-!is.element(ind$sex,"M")
         if(sum(wM)!=0) {
           minds<-ind[wM,]
           tmp<-data.frame("family"=rep(fm,nrow(minds)),"individ"=minds$individ, stringsAsFactors=FALSE)
           badsex<-rbind(badsex,tmp)
         }
     
  ##Testing for impossible relationships
	#child is parent of self or individual is child and parent of same person

	   rws<-NULL	
	   ui<-ped$individ
	   uin<-length(ui)
	
         moin<-paste(ped$mother,ped$individ)
         fain<-paste(ped$father,ped$individ)
	   for (j in 1:uin)
		{
			inmo<-paste(ped$individ[j],ped$mother[j])
			infa<-paste(ped$individ[j],ped$father[j])
                  w1<-which(is.element(moin,c(inmo,infa)))
                  w2<-which(is.element(fain,c(inmo,infa)))
                  w<-union(w1,w2)
                  if(length(w)!=0){
                      rw<-c(j,w); imrw<-sel[rw]
                      rws<-union(rws,imrw)
                  }
             }
         if(length(rws)!=0) {imposs[[i]]<-rws} 
                                                
      } # end of family loop
      
      imp<-unlist(lapply(imposs,is.null))
      if(any(!imp)){
          unms<-u[!imp]
          imposs<-imposs[!imp]
          names(imposs)<-unms
      } else imposs<-NULL

      bfam3<-bfam2
      if(!is.null(badsex)){
          tf<-badsex$family; bfam3<-union(bfam3,tf)
      }
      if(any(!imp)) bfam3<-union(bfam3,unms)

      out.listC<-list(badsex,imposs)
      lensC<-c(is.null(badsex),all(imp))
      nms<-c("mismatch.sex","impossible.related.rows")
      if(any(!lensC)){outC<-out.listC[!lensC];names(outC)<-nms[!lensC]} else {outC<-NULL}

 ###  subfamilies
     
     u<-unique(pedigree$family[!is.na(pedigree$family) & !is.element(pedigree$family,"") & !is.element(pedigree$family,bfam3)])

       if(length(u) ==0) {message(mssg); return(c(outA,outB,outC))}

     sfam<-NULL
     subfams<-NULL
         
     for(fm in u){
         sel<-which(is.element(pedigree$family,fm))
         ped<-pedigree[sel,]
 		#check relatedness
		  flag1<-0 ; counter<-0; Klus<-NULL
		  while(flag1==0) 
		  {
				ui<-setdiff(unique(ped$individ),Klus)
				clus<-NULL
				uclus<-NULL
				choice<-ui
				flag2<-0; counter<-counter+1
				  while(flag2==0) {
				   ind<-choice[1]
				   nclus<-relate(ind,ped)  #find everyone related to ind, includes ind
				   clus<-union(clus,nclus)  #all related so far
				   uclus<-union(uclus,ind)   #individuals already 'used' in relatedness check
				   choice<-setdiff(clus,uclus)
				   if(length(choice)==0) flag2<-1
				   if(all(is.element(ui,clus))) flag2<-1
				   } #end: finds a subfamily
				 Klus<-union(Klus,clus)
				 if(all(is.element(ped$individ,Klus))) 
					 	flag1<-1  #everyone accounted for
				 if(flag1==1 && counter==1) 
					 	break
				 sfam<-union(sfam,fm)
				 f<-rep(fm,length(clus))
				 sf<-rep(counter,length(clus))
				 newsf<-data.frame(f,sf,clus,stringsAsFactors=FALSE)
                         names(newsf)<-c("family","subfamily","individ")
				 subfams<-rbind(subfams,newsf)
		  } 
		 
       } # end family loop       

       if(!is.null(subfams)){
           subfams<-subfams[order(subfams$family,subfams$subfamily,subfams$individ),]
           outD<-list(subfams);names(outD)<-"subfamilies.ident"} else outD<-NULL
   
       out<-c(outA,outB,outC,outD)
       if(!is.null(out)) message(mssg)
       return(out)
}
    
   
        



