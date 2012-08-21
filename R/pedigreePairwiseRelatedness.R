###Relationship pairs with deeper pedigrees
#Input - pedigree dataframe with individ, mother, father, sex
#  Assumes (but checks) initial cleaning done and no one person families, 
# no mismatched mother/father sex, no impossible relationships
#Output - vector of families with inbreeding (to be handled by hand)
#         dataframe of Individ1, Individ2, relation,  kinship coefficient, family

#Options: opt=FALSE (default) creates pairs using only Individ1 and Individ2 from list of id's in individ column
#         opt=TRUE  creates pairs using any id's contained in individ, mother, or father           

pedigreePairwiseRelatedness<-
			function( pedigree,  use.any.ids=FALSE) 
{
	
	
		
	unrelated<-function(Y,A)
	{
		
		
		##Find kinship coefficients for each pair of individuals in a family
		## and identify unrelated pairs
		#Input: Y = data frame individ, mother,father for a given family
		#       A = adjacency matrix
		#Assume no dups, consistency checked
		#Output: $kinship is matrix of kinship coefficients
		#        $unrelated is dataframe of unrelated pairs
		#Here we assume the id labels in the dataframe Y are consecutive integers starting at 1
		# and include any mother/father id's not included in original individ list
		
			
		#find depth of pedigree along with 'ordering' individuals
		#depth is no. of gens after founder gen (corresponds to power where A^depth==0)
		B<-A
		v<-NULL
		n<-dim(A)[1]
		vleft<-1:n
		flag<-0
		depth<-0
		while(flag==0){
		if(all(B==0)){flag<-1;break}
		depth<-depth+1
		for (i in vleft){if (all(B[i,]==0)) v<-c(v,i)}
		vleft<-setdiff(1:n,v)
		B<-B%*%A    } #end of the while
		v<-c(v,vleft)
		
		##order according to v, then recode to 1:n so mother/father ids correspond
		Yo<-Y[v,]
		individ<-c(1:dim(Yo)[1])
		mother<-match(Yo$mother,Yo$individ,nomatch=0)
		father<-match(Yo$father,Yo$individ,nomatch=0)
		YY<-data.frame(individ,mother,father,stringsAsFactors=FALSE)
		
		#Calculate kinship coefficient for each pair
		kc<-diag(1/2,n,n)
		n1<-n-1
		for (i in 1:n1) {i1<-i+1
		  for (j in i1:n) {
		  jm<-YY$mother[j]; jf<-YY$father[j]
		  if (jm==0) kcm<-0 else kcm<-kc[i,jm]
		  if(jf==0) kcf<-0 else kcf<-kc[i,jf]
		  kc[i,j]<-1/2*(kcm+kcf)
		  kc[j,i]<-kc[i,j]
		} }
		## recode to original input individual ids
		KC<-matrix(0,n,n)
		for (i in 1:n){
		 for (j in 1:n){ KC[v[i],v[j]]<-kc[i,j] } }
		
		##find unrelated pairs kinship coeff = 0
		tKC<-KC
		tKC[row(tKC)>col(tKC)]<-1
		wunr<-which(tKC==0,arr.ind=TRUE)
		unr<-data.frame(wunr,stringsAsFactors=FALSE)
		names(unr)<-c("Individ1","Individ2")
		out.list<-list(KC,unr)
		names(out.list)<-c("kinship","unrelated")
		return(out.list) 
	 }

	
	
	
	samp<- pedigree
	opt<- use.any.ids
	er<- pedigreeCheck(samp)
	if(!is.null(er)){stop("ERROR: There are consistency errors. Run pedigreeCheck to diagnose")}
	  
	u<-unique(samp$family)
	un<-length(u)
	relativeprs<-NULL
	inbreed<-NULL
      inbred.kc<-NULL
	for (i in 1:un){
	x<-samp[is.element(samp$family,u[i]),c("individ","mother","father")] #get family
	
	#add 'extra' mothers and fathers to individ list  (assumed to be "founders")
	ui<-x$individ  #recall no duplicates
	nind<-length(ui)
	um<-x$mother
	uf<-x$father
	individ<-x$individ
	mother<-x$mother
	father<-x$father


	#find non-zero mother/father (known) entries that are not in individ list
	um<-setdiff(um,0); uf<-setdiff(uf,0)
	em<-setdiff(um,intersect(ui,um))
	ef<-setdiff(uf,intersect(ui,uf))
	new<-c(em,ef);ln<-length(new)

	#identify children with one unknown parent, assign individ id and assume is 'founder' (no parents)
	wum<-which(is.element(x$mother,0)) 
	nwuf<-which(!is.element(x$father,0))
	wm<-intersect(wum,nwuf)
	wuf<-which(is.element(x$father,0))
	nwum<-which(!is.element(x$mother,0))
	wf<-intersect(wuf,nwum)
	lm<-length(wm); lf<-length(wf)

      if(lf!=0){ exf<-paste("unfa",individ[wf],sep=""); father[wf]<-exf} else exf<-NULL
	if(lm!=0){exm<-paste("unmo",individ[wm],sep="");mother[wm]<-exm} else exm<-NULL

	if(ln !=0) {
	individ<-c(individ,new)
	mother<-c(mother,rep(0,ln))
	father<-c(father,rep(0,ln))}
	if(lf!=0){
	individ<-c(individ,exf)
	mother<-c(mother,rep(0,lf))
	father<-c(father,rep(0,lf))}
	if(lm!=0){
	individ<-c(individ,exm)
	mother<-c(mother,rep(0,lm))
	father<-c(father,rep(0,lm))}
	
	XX<-data.frame(individ,mother,father,stringsAsFactors=FALSE)  
	
	#recode so that individ is 1:number of individuals, mother/father ids correspond
	individ<-c(1:dim(XX)[1])
	mother<-match(XX$mother,XX$individ,nomatch=0)
	father<-match(XX$father,XX$individ,nomatch=0)
	Y<-data.frame(individ,mother,father,stringsAsFactors=FALSE)
	
	#find offspring,parent pairs directly from pedigree
	p<-Y[!is.element(Y$mother,0),c("individ","mother")]
	q<-Y[!is.element(Y$father,0),c("individ","father")]
	names(p)<-c("offspring","parent")
	names(q)<-c("offspring","parent")
	po<-rbind(p,q)
	
	#find adjacency matrix
	n<-dim(Y)[1]
	A<-matrix(0,n,n)
	ipo<-as.matrix(po)
	A[ipo]<-1
	
	#Find indices of unrelated pairs
	dg<-unrelated(Y,A)
	unr<-dg$unrelated
	punr<-paste(unr$Individ1,unr$Individ2)
	
	##Flag families that have inbreeding
	w<-which(!is.element(Y$mother,0) & !is.element(Y$father,0))
	mofa<-paste(Y$mother[w],Y$father[w])
	famo<-paste(Y$father[w],Y$mother[w])
	if(any(!is.element(mofa,punr) & !is.element(famo,punr))) {
            inbreed<-c(inbreed,u[i])
           	pprs<-combn(Y$individ,2)
	      tp<-t(pprs)
	      inbprs<-data.frame(tp,stringsAsFactors=FALSE)
	      names(inbprs)<-c("Individ1","Individ2")
         	inbprs$kinship<-dg$kinship[tp]
            inbprs$family<-u[i]
	#decode back to original individ id's
	     if(!opt)uui<-ui else uui<-c(ui,new,exf,exm)
	     w1<-as.list(rep(NA,length(uui)))
	     w2<-as.list(rep(NA,length(uui)))
	     for (j in 1:length(uui)) {
	         w1[[j]]<-which(is.element(inbprs$Individ1,j))
	         w2[[j]]<-which(is.element(inbprs$Individ2,j))
           }
	     for (j in 1:length(uui)){
	         inbprs$Individ1[w1[[j]]]<-uui[j]
	         inbprs$Individ2[w2[[j]]]<-uui[j]  
           }
           inbred.kc<-rbind(inbred.kc,inbprs)
           next
      }
	
	###SIBS
	#Find which parents have more than one child - columns of adjacency matrix have more than 1 one
	FS<-NULL
	HS<-NULL
	
	HSr<-NULL
	
	s<-rep(NA,n)
	for (j in 1:n) s[j]<-sum(A[,j])
	wpar<-which(s>1)
	if(length(wpar)!=0) {
	## Identify sib pairs
	
	for(j in 1:length(wpar)){ 
	 cww<-which(A[,wpar[j]]==1) #identify children of a given parent
	 
	  pp<-combn(cww,2); pn<-dim(pp)[2]
	  for (k in 1:pn) {
	   z<-Y[is.element(Y$individ,pp[,k]),c("mother","father")]
	   mu<-length(unique(z$mother))
	   if(mu==1) zq<-z$father else zq<-z$mother
	   if(length(unique(zq))==1) {FS<-rbind(FS,pp[,k]);next} 
	   
	   p12<-paste(zq[1],zq[2]); p21<-paste(zq[2],zq[1])
	   if (is.element(p12,punr) | is.element(p21,punr)) HS<-rbind(HS,pp[,k]) else HSr<-rbind(HSr,pp[,k])  
	   } 
	  } } #end of building up FS,HS,HSr
	
	
	# There will be duplicates but that will not affect the final output matrix
	
	#end finding sibs
	
	##AVUNCULAR
	avF<-NULL
	
	#full Avuncular
	if(length(FS)!=0) fs<-dim(FS)[1] else fs<-0      
	if(fs!=0) {
	 for (ii in 1:fs) { au1<-FS[ii,1];au2<-FS[ii,2]  #potential aunts or uncles
	  c1<-which(A[,au1]==1) #children of au1
	  c2<-which(A[,au2]==1) #children of au2
	  au<-c(rep(au1,length(c2)),rep(au2,length(c1)))
	  nn<-c(c2,c1)
	  temp<-cbind(au,nn)
	  avF<-rbind(avF,temp) }
	 } #end of full avuncular
	
	#half Avuncular
	avH<-NULL
	avO<-NULL
	if(length(HS)!=0) hs<-dim(HS)[1] else hs<-0
	if(hs!=0) {
	 for (ii in 1:hs) { 
	cc1<-NULL; cc2<-NULL;co1<-NULL;co2<-NULL
	au1<-HS[ii,1];au2<-HS[ii,2]  #potential aunts or uncles
	  c1<-which(A[,au1]==1) #children of au1
	
	  if(length(c1)!=0){
		for (j in 1:length(c1)) {m<-Y$mother[c1[j]]; f<-Y$father[c1[j]]
		  if(m==au1) chk<-f else chk<-m
		  an<-paste(au2,chk);na<-paste(chk,au2)
		  if(is.element(an,punr)|is.element(na,punr)) cc1<-c(cc1,c1[j]) else co1<-c(co1,c1[j])
		 } }
	  c2<-which(A[,au2]==1) #children of au2
	  if(length(c2)!=0){
		for (j in 1:length(c2)) {m<-Y$mother[c2[j]]; f<-Y$father[c2[j]]
		  if(m==au2) chk<-f else chk<-m
		  an<-paste(au1,chk);na<-paste(chk,au1)
		  if(is.element(an,punr)|is.element(na,punr)) cc2<-c(cc2,c2[j]) else co2<-c(co2,c2[j])
		 } }
	   au<-c(rep(au1,length(cc2)),rep(au2,length(cc1)))
	   nn<-c(cc2,cc1)
	   temp<-cbind(au,nn)
	   avH<-rbind(avH,temp)
	   auo<-c(rep(au1,length(co2)),rep(au2,length(co1)))
	   nno<-c(co2,co1)
	   tempo<-cbind(auo,nno)
	   avO<-rbind(avO,tempo)
	 } }#end of half avuncular with some identification of 'other' avuncular
	
	#Other avuncular
	
	osib<-HSr
	if(length(osib)!=0) os<-dim(osib)[1] else os<-0      
	if(os!=0) {
	 for (ii in 1:os) { au1<-osib[ii,1];au2<-osib[ii,2]  #potential aunts or uncles
	  c1<-which(A[,au1]==1) #children of au1
	  c2<-which(A[,au2]==1) #children of au2
	  au<-c(rep(au1,length(c2)),rep(au2,length(c1)))
	  nn<-c(c2,c1)
	  temp<-cbind(au,nn) 
	  avO<-rbind(avO,temp)}
	 } #end of other avuncular
	
	###GRANDPARENT/GRANDCHILD
	A2<-A%*%A  #entry is 1 if path of length 2 from i to j; no entries>1 if no inbreeding
	s<-rep(NA,n)
	for (j in 1:n) s[j]<-sum(A2[,j])
	wpar<-which(s>=1)   #potential grandparents
	gpgc<-NULL
	if(length(wpar)!=0) {
	for(j in 1:length(wpar)){ 
	 cww<-which(A2[,wpar[j]]==1) #identify grandchildren of a given grandparent
	temp<-cbind(rep(wpar[j],length(cww)),cww)
	gpgc<-rbind(gpgc,temp) }
	 }
	
	###COUSINS
	#secgen is a vector of all second generation pairs
	secgen<-punr
	pun<-length(punr)  #punr is the pasted unrelateds
	fsn<-0;hsn<-0;hsnr<-0
	if(length(FS)!=0) {FSpr<-paste(FS[,1],FS[,2]); fsn<-length(FSpr);secgen<-c(secgen,FSpr)} else fsn<-0
	
	if(length(HS)!=0) {HSpr<-paste(HS[,1],HS[,2]); hsn<-length(HSpr);secgen<-c(secgen,HSpr)} else hsn<-0
	if(length(HSr)!=0) {HSprr<-paste(HSr[,1],HSr[,2]); hsnr<-length(HSprr); secgen<-c(secgen,HSprr)} else hsnr<-0
	
	#Identifies range of positions in secgen where each relation 'sits'
	Un<-pun; pu<-1:Un
	FSn<-pun+fsn
	if (fsn==0) fs<-NULL else fs<-(pun+1):FSn
	HSn<-FSn+hsn
	if(hsn==0) hs<-NULL else hs<-(FSn+1):HSn
	HSnr<-HSn+hsnr;
	if(hsnr==0) hsr<-NULL else hsr<-(HSn+1):HSnr
	
	fcous<-NULL
	hfcous<-NULL
	dfcous<-NULL
	ocous<-NULL
	
	#establish potential cousin pairs to examine
	s<-rep(NA,n)
	for (j in 1:n) s[j]<-sum(A2[,j])
	wpar<-which(s>1)   #potential grandparents with more than one grandchild
	if(length(wpar)!=0) {
	
	for(j in 1:length(wpar)){ 
	 cww<-which(A2[,wpar[j]]==1) #identify grandchildren of a given grandparent (at least 2)
	  pp<-combn(cww,2); pn<-dim(pp)[2]
	 for (K in 1:pn) {M<-Y$mother[pp[,K]];F<-Y$father[pp[,K]]
	   m1<-M[1];m2<-M[2];f1<-F[1];f2<-F[2]
	  
	   if(m1==m2 | f1==f2) next  #sibs - already accounted for
	   if(m1<m2) mm<-paste(m1,m2) else mm<-paste(m2,m1)
	   if(m1<f2) mf<-paste(m1,f2) else mf<-paste(f2,m1)
	   if(m2<f1) fm<-paste(m2,f1) else fm<-paste(f1,m2)
	   if(f1<f2) ff<-paste(f1,f2) else ff<-paste(f2,f1)
	
	#Decide relationship of these parent (of grandchildren) pairs  - Note secgen covers all possibilities
	   Mmm<-match(mm,secgen); Mmf<-match(mf,secgen) ; Mfm<-match(fm,secgen); Mff<-match(ff,secgen)
	#match will give position in secgen that is first match for first argument
	   rpar<-c(Mmm,Mmf,Mfm,Mff)
	
	#decide which (and how many) of the parent combinations are unrelated
	#at most three are unrelated since at least one of the parents in each parent pair is connected to common grandparent
	   wu<-which(is.element(rpar,pu)) 
	   woth<-setdiff(1:4,wu)  #which parent combinations are related somehow 
	   if(length(wu)==3){
		 if(is.element(rpar[woth],fs)) {fcous<-rbind(fcous,pp[,K]);next}  #first cousins since one combo is full sib
		 if(is.element(rpar[woth],hs)){hfcous<-rbind(hfcous,pp[,K]);next}}
	
	   if(length(wu)==2){
		 if(all(is.element(rpar[woth],fs))) {dfcous<-rbind(dfcous,pp[,K]);next}}
	
		 ocous<-rbind(ocous,pp[,K]) 
	 } #end of loop on K
	  }#end of loop on j 
	 }#end of if (length(wpar ..)
	
	 # There will be duplicates but this will not affect final output matrix
	
	#end finding cousins 
	
	#Half sib + first cousin
	hsfc<-NULL
	if(length(HSr)!=0 && length(FS)!=0){ hsr<-dim(HSr)[1]
	kklist<-NULL
	for (kk in 1:hsr) {
		M<-Y$mother[HSr[kk,]];F<-Y$father[HSr[kk,]]
	   m1<-M[1];m2<-M[2];f1<-F[1];f2<-F[2]
	   if(m1==m2) {chk1<-paste(f1,f2);chk2<-paste(f2,f1)} else {chk1<-paste(m1,m2);chk2<-paste(f1,f2)}
	   if(is.element(chk1,FSpr) | is.element(chk2,FSpr)) {kklist<-c(kklist,kk) ;hsfc<-rbind(hsfc,HSr[kk,]) } 
	}#end loop
	
	#Delete these specially identified half sibs from the HSr list
	
	w<-1:hsr
	ww<-setdiff(w,kklist)
	tHSr<-HSr[ww,]
	if(length(tHSr)==2) tHSr<-matrix(c(tHSr[1],tHSr[2]),1,2)
	HSr<-tHSr 
	  
	}#end if on HSr
	
	
	
	#FIND RELATED PAIR MATRIX
	#identify "Other" as the rest of the pairs
	#create dataframe with all pairs and column for type of relationship
	
	pprs<-combn(Y$individ,2)
	tp<-t(pprs)
	relprs<-data.frame(tp,stringsAsFactors=FALSE)
	names(relprs)<-c("Individ1","Individ2")
	relprs$relation<-rep("Other",dim(pprs)[2])
	relprs$kinship<-dg$kinship[tp]
	R<-paste(relprs$Individ1,relprs$Individ2)
	
	mU<-match(punr,R);relprs$relation[mU]<-"U"
	
	ppo<-paste(po[,1],po[,2])
	pop<-paste(po[,2],po[,1])
	mpo<-match(ppo,R) #gives indices relprs that match with PO
	mmpo<-match(pop,R)
	relprs$relation[mpo]<-"PO"
	relprs$relation[mmpo]<-"PO"
	
	if(length(FS)!=0){
	FSpr<-paste(FS[,1],FS[,2]) 
	mFS<-match(FSpr,R) #gives indices relprs that match with FS
	relprs$relation[mFS]<-"FS"}
	if(length(HS)!=0){
	HSpr<-paste(HS[,1],HS[,2]) 
	mHS<-match(HSpr,R) #gives indices relprs that match with HS
	relprs$relation[mHS]<-"HS"}
	if(length(HSr)!=0){
	HSprr<-paste(HSr[,1],HSr[,2])
	mHSr<-match(HSprr,R) #gives indices relprs that match with HSr
	relprs$relation[mHSr]<-"HSr"}
	
	if(length(avF)!=0) {
	 pavF<-paste(avF[,1],avF[,2])
	 PavF<-paste(avF[,2],avF[,1])
	 mavF<-match(pavF,R); MavF<-match(PavF,R)
	 relprs$relation[mavF]<-"Av" ;relprs$relation[MavF]<-"Av" }
	if(length(avH)!=0) {
	 pavH<-paste(avH[,1],avH[,2])
	 PavH<-paste(avH[,2],avH[,1])
	 mavH<-match(pavH,R); MavH<-match(PavH,R)
	 relprs$relation[mavH]<-"HAv"; relprs$relation[MavH]<-"HAv" }
	if(length(avO)!=0) {
	 pavO<-paste(avO[,1],avO[,2])
	 PavO<-paste(avO[,2],avO[,1])
	 mavO<-match(pavO,R);MavO<-match(pavO,R)
	 relprs$relation[mavO]<-"OAv"; relprs$relation[MavO]<-"OAv" }
	
	if(length(gpgc)!=0){
	  pgpgc<-paste(gpgc[,1],gpgc[,2])
	  Pgpgc<-paste(gpgc[,2],gpgc[,1])
	  mgpgc<-match(pgpgc,R);Mgpgc<-match(Pgpgc,R)
	  relprs$relation[mgpgc]<-"GpGc"; relprs$relation[Mgpgc]<-"GpGc" }
	
	if(length(fcous)!=0){
	  pfcous<-paste(fcous[,1],fcous[,2])
	  mfcous<-match(pfcous,R)
	  relprs$relation[mfcous]<-"FC" }
	if(length(hfcous)!=0){
	  phfcous<-paste(hfcous[,1],hfcous[,2])
	  mhfcous<-match(phfcous,R)
	  relprs$relation[mhfcous]<-"HFC" }
	if(length(dfcous)!=0){
	  pdfcous<-paste(dfcous[,1],dfcous[,2])
	  mdfcous<-match(pdfcous,R)
	  relprs$relation[mdfcous]<-"DFC" }
	if(length(ocous)!=0){
	  pocous<-paste(ocous[,1],ocous[,2])
	  mocous<-match(pocous,R)
	  relprs$relation[mocous]<-"OC" }
	
	if(length(hsfc)!=0){
	  phsfc<-paste(hsfc[,1],hsfc[,2])
	  mhsfc<-match(phsfc,R)
	  relprs$relation[mhsfc]<-"HSFC" }
	
	#delete any pairs that involve the 'added' individs
	#recall the id's in mother or father are contained in new 
	#exf and exm contain extra id's for missing identifiers
	#if opt=FALSE, delete all
	#if opt=TRUE, delete only the extras for missing identifiers
	
	tnew<-NULL
	ladd<-length(new)+length(exf)+length(exm)
	lex<-length(c(exf,exm))
	lui<-length(ui)+1
	lp<-length(ui)+length(new)
	lend<-length(ui)+ladd
	if(!opt && (ladd !=0)) tnew<-lui:lend
	if(opt && (lex!=0)) tnew<-(lp+1):lend 
	w1<-which(is.element(relprs[,1],tnew))
	w2<-which(is.element(relprs[,2],tnew))
	w<-union(w1,w2)
	W<-c(1:dim(relprs)[1])
	keep<-setdiff(W,w)
	relprs<-relprs[keep,]
	
	#decode back to original individ id's
	if(!opt)uui<-ui else uui<-c(ui,new,exf,exm)
	w1<-as.list(rep(NA,length(uui)))
	w2<-as.list(rep(NA,length(uui)))
	for (j in 1:length(uui)) {
	 w1[[j]]<-which(is.element(relprs$Individ1,j))
	 w2[[j]]<-which(is.element(relprs$Individ2,j)) }
	for (j in 1:length(uui)){
	 relprs$Individ1[w1[[j]]]<-uui[j]
	 relprs$Individ2[w2[[j]]]<-uui[j]  }
	
	#add family id column
	relprs$family<-rep(u[i],dim(relprs)[1])
	
	
	#add onto previous family
	relativeprs<-rbind(relativeprs,relprs)
	} #end of family loop
	out.list<-list(inbreed,inbred.kc,relativeprs)
	names(out.list)<-c("inbred.fam","inbred.KC","relativeprs")
	return(out.list)
}

