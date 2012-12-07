pedigreeMaxUnrelated<-function(pedigree,pref=NULL){

##INPUT: 'pedigree' is a pedigree dataframe with family, individ, mother, father, sex and selset columns. 
 # The selset column indicates whether the individual is in the specified subset of interest (1) or not (0)
 # optional 'pref' column is integer with lower values indicating higher preference
##OUTPUT: a dataframe such that $family is the family id and corresponding
#   entries in $Individ are the individual identifiers of a maximal set of individuals in specified subset, 
#   such that all pairs of individuals are unrelated
# Note that such maximal sets are not unique

# It is assumed that pedigree is the full known pedigree - partial pedigrees could falsely indicate unrelated relationship
# An error message is generated if there are pedigree inconsistencies. 
#   This includes (among other things) one person families, parent IDs not included as individual entries, and unknown parents.
#   The user will need to add one-person families to the output.

################ SUBFUNCTION KCpairs ######################
KCpairs<-function(Y,A){
#Input: Y = data frame individ, mother,father for a given family
#       A = adjacency matrix
#Assume no dups, consistency checked
#Output: kinship coefficient matrix
#Here we assume the id labels in the dataframe Y are consecutive integers starting at 1

  #find depth of pedigree along with 'ordering' individuals
  #depth is no. of gens after founder gen (corresponds to power where A^depth==0)
  B<-A
  v<-NULL
  if(!is.matrix(A))stop("Second input is not a matrix")
  n<-dim(A)[1]
  vleft<-1:n
  flag<-0
  depth<-0
  while(flag==0){
    if(all(B==0)){flag<-1;break}
    depth<-depth+1
    for (i in vleft){if (all(B[i,]==0)) v<-c(v,i)}
    vleft<-setdiff(1:n,v)
    B<-B%*%A    
  } #end of the while
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
    } 
  }
## recode to original input individual ids
  KC<-matrix(0,n,n)
  for (i in 1:n){
    for (j in 1:n){ KC[v[i],v[j]]<-kc[i,j] 
    } 
  }
  return(KC)
 
}
####################
################## SUBFUNCTION fam.maxset.unrelated ############
fam.maxset.unrelated<-function(KC,g,p){

## Input is KC kinship coefficient matrix 
##  g is vector of indices for selected individuals
##  p is vector of preference values corresponding to g
##Output maxset is a maximal set of indices such that 
#  all pairs are unrelated; 
#  if no pair of selected individuals are unrelated, 
#    output is one selected individual in the family

  pfval<-sort(unique(p))
  M<-KC[g,g]
## index j of M corresponds with index g[j] for KC (hence with individ g[j])
  m<-length(g)
  un<-list();ct<-vector()
##for each selected individ, count number of unrelated individs
#  order from largest count to smallest count
  for (j in 1:m) { 
    un[[j]]<-which(M[,j]==0)
    ct[j]<-length(un[[j]])
  }
    ind<-1:m
    mu<-cbind(ind,ct,p)
    uct<-sort(unique(ct),decreasing=TRUE)
    mu2<-NULL   # order by ct but retain ordering within a count group
    for(ii in 1:length(uct)){
       kt<-uct[ii]
       sel<-is.element(mu[,2],kt)
       mu2<-rbind(mu2,mu[sel,])
    }
    mud<-data.frame(mu2)
    X<-mud[!is.element(mud$ct,c(0,1)),]
    mn<-length(X$ind)

         ## if there are no individs in g with 2 or more unrelateds
    if(mn==0){
      sel<-is.element(mud$ct,1)
      if(sum(sel)==0){ # all selected are related
        J<-mud$ind[1]
        maxsetf<-g[J]
        return(maxsetf)
      } else {
      temp<-mud[is.element(mud$ct,1),]
      J<-temp$ind[1]
      maxsetf<-c(g[J],g[un[[J]]])
      return(maxsetf)
      }      
    }

##consider the individuals with 2 or more unrelateds
    nms<-1;jj<-1; Mmaxset<-NULL
    while(X$ct[jj]>nms && jj<=mn){
        flag<-0; id<-X$ind[jj]
        ss<-un[[id]];nss<-length(ss)
        for(k in nss:2) {
           cb<-combn(ss,k)
           for (L in 1:dim(cb)[2]) {
             gg<-cb[,L];ng<-length(gg)
             tM<-M[gg,gg]; D<-diag(.5,ng,ng)
             if(all(tM-D==0)) {
                maxset<-c(id,gg);nms<-length(gg) 
                if(length(maxset)>length(Mmaxset)) Mmaxset<-maxset
                if(length(pfval)>1){
                  if(length(maxset)==length(Mmaxset)) {
                     pvnew<-p[maxset]
                     pvold<-p[Mmaxset]
                     for(pp in pfval[1:(length(pfval)-1)]){
                        pn<-sum(pvnew==pp)
                        po<-sum(pvold==pp)
                        if(pn > po) {Mmaxset<-maxset; break}
                        if(pn<po) break
                     }
                  }
                }

                flag<-1
             }
             if(flag==1)break
           }
           if(flag==1) break
        }
        if(flag==0){ 
          maxset<-c(id,ss[1]);nms<-1
          if(length(maxset)>length(Mmaxset)) Mmaxset<-maxset
          if(length(pfval)>1){
             if(length(maxset)==length(Mmaxset)) {
                pvnew<-p[maxset]
                pvold<-p[Mmaxset]
                for(pp in pfval[1:(length(pfval)-1)]){
                  pn<-sum(pvnew==pp)
                  po<-sum(pvold==pp)
                  if(pn > po) {Mmaxset<-maxset; break}
                  if(pn<po) break
                }
             }
           }

        }
        jj<-jj+1
    }
    maxsetf<-g[Mmaxset]; return(maxsetf)
}
##################################################
###### MAIN PROGRAM ###################
if(!is.element(class(pedigree),"data.frame")) stop("Error: input must be a data.frame")
req<-c("family","individ","mother","father","sex","selset")
reqn<-paste(req,collapse = ", ")
if(!all(is.element(req,names(pedigree)))) stop(paste("Error: some of the required variable names (",reqn,") are not present.",sep=""))
if(!all(is.element(pedigree$selset,c(0,1)))) stop("Error: some values for 'selset' variable are not 0 or 1")

if(!is.null(pref)){
    if(!is.element(class(pref),"character")) stop("preference variable name must be character or NULL")
    if(!is.element(class(pedigree[,pref]),c("integer","numeric"))) stop("preference variable type needs to be 'integer' or 'numeric'")
    if(any(is.na(pedigree[,pref]))) stop("preference variable can not have NA values")
    if(!all(pedigree[,pref]>=1)) stop("preference variable values must be >=1")
}
chk<-pedigreeCheck(pedigree)
if(!is.null(chk)) stop("Error: there are pedigree inconsistencies. Run pedigreeCheck to diagnose")


u<-unique(pedigree$family)
un<-length(u)
max.ids<-NULL

for (i in 1:un){
   x<-pedigree[is.element(pedigree$family,u[i]),] #get family
   if(is.null(pref)){
       x$pref<-2
       sel<-x$mother==0 & x$father==0
       x$pref[sel]<-1
   }
      
   x<-x[order(x$pref),]

   ##If there is only one selected individual in family, record that person,go to next family
    # If there are no selected individuals, go to next family

   if(sum(x$selset==1)==1){
       fam.id<-u[i];Individ<-x$individ[x$selset==1];mx<-data.frame(fam.id,Individ)
       max.ids<-rbind(max.ids,mx);next
   }
   if(sum(x$selset==1)==0) next

   ### Work with pedigree to find kinship coefficients for every pair of individuals in the family

   XX<-x[,c("individ","mother","father")]  

    #recode so that individ is 1:number of individuals, mother/father ids correspond
   individ<-c(1:dim(XX)[1])
   mother<-match(XX$mother,XX$individ,nomatch=0)
   father<-match(XX$father,XX$individ,nomatch=0) 
    # 'no matches' should be only the founders so want 0 - pedigreeCheck should take care of this    
   Y<-data.frame(individ,mother,father,stringsAsFactors=FALSE) # but should be not strings

   #find offspring,parent pairs directly from pedigree
   p<-Y[!is.element(Y$mother,0),c("individ","mother")]
   q<-Y[!is.element(Y$father,0),c("individ","father")]
   names(p)<-c("offspring","parent")
   names(q)<-c("offspring","parent")
   po<-rbind(p,q)

  #find adjacency matrix
   n<-nrow(Y)
   A<-matrix(0,n,n)
   ipo<-as.matrix(po)
   A[ipo]<-1

  #Find kinship coefficient matrix and selected individs
   KC<-KCpairs(Y,A)
   g<-which(x$selset==1)
   p<-x$pref[g]   # preference levels
  
  ##Find maximal unrelated set of individuals in family i
   mset<-fam.maxset.unrelated(KC,g,p)

  #Decode back to original ids
   msetid<-x$individ[mset]
   ms<-length(msetid)

  #add to max.ids dataframe
  fam.id<-rep(u[i],ms)
  Individ<-msetid
  mx<-data.frame(fam.id,Individ,stringsAsFactors=FALSE)
  max.ids<-rbind(max.ids, mx)
}   # end fam loop
  names(max.ids)<-c("family","Individ")
  return(max.ids)
} #function end

