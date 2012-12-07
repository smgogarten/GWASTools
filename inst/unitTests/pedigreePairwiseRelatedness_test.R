test_pedigreePairwiseRelatedness <- function() {
## pedigreePairwiseRelatedness

family <- c(1,1,1,1,2,2,2,2,2)
individ <- c(1,2,3,4,5,6,7,8,9)
mother <- c(0,0,1,1,0,0,5,5,0)
father <- c(0,0,2,2,0,0,6,9,0)
sex <- c("F","M","F","F","F","M","M","M","M")
samp <- data.frame(family, individ, mother, father, sex, stringsAsFactors=FALSE)
pairs1<-t(combn(samp$individ[samp$family==1],2))
pairs2<-t(combn(samp$individ[samp$family==2],2))
pairs<-rbind(pairs1,pairs2)
prs<-paste(pairs[,1],pairs[,2])
relprs<-data.frame("Individ1"=pairs[,1],"Individ2"=pairs[,2],stringsAsFactors=FALSE)
relprs$relation<-"U"
relprs$kinship<-0
po<-c("1 3","1 4", "2 3","2 4", "5 7", "5 8","6 7","8 9")
FS<-"3 4"
HS<-"7 8"
posel<-is.element(prs,po)
relprs$kinship[posel]<-0.25
relprs$relation[posel]<-"PO"
fssel<-is.element(prs,FS)
relprs$kinship[fssel]<-0.25
relprs$relation[fssel]<-"FS"
hssel<-is.element(prs,HS)
relprs$kinship[hssel]<-0.125
relprs$relation[hssel]<-"HS"
relprs$family<-c(rep(1,choose(4,2)),rep(2,choose(5,2)))

expected.result<-list(inbred.fam=NULL,inbred.KC=NULL,relativeprs=relprs)
result<-pedigreePairwiseRelatedness(samp)
checkEquals(expected.result,result)
#######################################

# used for example
family <- c(1,1,1,1,2,2,2,2,2,2,2)
individ <- c(1,2,3,4,5,6,7,8,9,10,11)
mother <- c(0,0,1,1,0,0,5,5,0,0,10)
father <- c(0,0,2,2,0,0,6,9,0,0,7)
sex <- c("F","M","F","F","F","M","M","M","M","F","F")
samp <- data.frame(family, individ, mother, father, sex, stringsAsFactors=FALSE)
pairs1<-t(combn(samp$individ[samp$family==1],2))
pairs2<-t(combn(samp$individ[samp$family==2],2))
n1<-nrow(pairs1)
n2<-nrow(pairs2)
pairs<-rbind(pairs1,pairs2)
prs<-paste(pairs[,1],pairs[,2])
relprs<-data.frame("Individ1"=pairs[,1],"Individ2"=pairs[,2],stringsAsFactors=FALSE)
relprs$relation<-"U"
relprs$kinship<-0
po<-c("1 3","1 4", "2 3","2 4", "5 7", "5 8","6 7","8 9", "7 11","10 11")
FS<-"3 4"
HS<-"7 8"
Hav<-"8 11"
gpgc<-c("5 11","6 11")
posel<-is.element(prs,po)
relprs$kinship[posel]<-0.25
relprs$relation[posel]<-"PO"
fssel<-is.element(prs,FS)
relprs$kinship[fssel]<-0.25
relprs$relation[fssel]<-"FS"
hssel<-is.element(prs,HS)
relprs$kinship[hssel]<-0.125
relprs$relation[hssel]<-"HS"
havsel<-is.element(prs,Hav)
relprs$kinship[havsel]<-0.0625
relprs$relation[havsel]<-"HAv"
gpsel<-is.element(prs,gpgc)
relprs$kinship[gpsel]<-0.125
relprs$relation[gpsel]<-"GpGc"
relprs$family<-c(rep(1,n1),rep(2,n2))

expected.result<-list(inbred.fam=NULL,inbred.KC=NULL,relativeprs=relprs)
result<-pedigreePairwiseRelatedness(samp)
checkEquals(expected.result,result)

###############################

#inbreeding
family <- rep(2,7)
individ <- paste("I",c(1,2,3,4,5,6,7),sep="")
mother <- c(0,0,0,"I1","I1","I3","I5")
father <- c(0,0,0,"I2","I2","I4","I4")
sex <- c("F","M","F","M","F","F","F")
samp <- data.frame(family, individ, mother, father, sex, stringsAsFactors=FALSE)
pairs<-t(combn(samp$individ,2))
ks<-c(0,0,.25,.25,.125,.25,0,.25,.25,.125,.25,0,0,.25,0,.25,.25,.3750,.125,.3750,.1875)
fm<-rep(2,nrow(pairs))
kcpairs<-data.frame("Individ1"=pairs[,1],"Individ2"=pairs[,2],"kinship"=ks,"family"=fm,stringsAsFactors = FALSE)
expected.result<-list(inbred.fam=2,inbred.KC = kcpairs,relativeprs=NULL)
result<-pedigreePairwiseRelatedness(samp)
checkEquals(expected.result,result)
#################################

## somewhat more complicated
family<-rep(1,10)
individ<-paste("I",1:10,sep="")
mother<-c(0,0,"I1","I1",0,0,"I6","I6","I4","I7")
father<-c(0,0,"I2","I2",0,0,"I5","I5","I8","I3")
sex<-c("F","M","M","F","M","F","F","M","F","F")
samp<-data.frame(family,individ,mother,father,sex,stringsAsFactors=F)
pairs<-t(combn(samp$individ,2))
prs<-paste(pairs[,1],pairs[,2])
relprs<-data.frame("Individ1"=pairs[,1],"Individ2"=pairs[,2],stringsAsFactors=FALSE)
relprs$relation<-"U"
relprs$kinship<-0
PO<-c("I1 I3", "I1 I4", "I2 I3","I2 I4","I5 I7", "I5 I8","I6 I7","I6 I8", "I3 I10","I7 I10","I4 I9", "I8 I9")
FS<-c("I3 I4", "I7 I8")
AV<-c("I4 I10","I3 I9","I8 I10", "I7 I9")
dfc<-"I9 I10"
gpgc<-c("I1 I9","I1 I10","I2 I9","I2 I10","I5 I10","I6 I10", "I5 I9","I6 I9")
posel<-is.element(prs,PO)
fssel<-is.element(prs,FS)
avsel<-is.element(prs,AV)
dfcsel<-is.element(prs,dfc)
gpsel<-is.element(prs,gpgc)
relprs$relation[posel]<-"PO"
relprs$relation[fssel]<-"FS"
relprs$relation[avsel]<-"Av"
relprs$relation[dfcsel]<-"DFC"
relprs$relation[gpsel]<-"GpGc"
relprs$kinship[posel]<-.25
relprs$kinship[fssel]<-.25
relprs$kinship[avsel]<-.125
relprs$kinship[dfcsel]<-.125
relprs$kinship[gpsel]<-.125

relprs$family<-rep(1,nrow(pairs))
expected.result<-list("inbred.fam"=NULL,"inbred.KC"=NULL, relativeprs=relprs)
result<-pedigreePairwiseRelatedness(samp) # relationships check out
checkEquals(expected.result,result)
#####################################

family<-rep(1,12)
individ<-1:12
mother<-c(0,0,1,1,0,5,0,4,0,9,0,6)
father<-c(0,0,2,2,0,2,0,7,0,3,0,11)
sex<-c("F","M","M","F","F","F","M","F","F","F","M","F")
samp<-data.frame(family,individ,mother,father,sex,stringsAsFactors=FALSE)
pairs<-t(combn(samp$individ,2))
prs<-paste(pairs[,1],pairs[,2])
relprs<-data.frame("Individ1"=pairs[,1],"Individ2"=pairs[,2],stringsAsFactors=FALSE)
relprs$relation<-"U"
relprs$kinship<-0
po1<-c(1,1,2,2,2,5,3,9,6,11,4,7)
po2<-c(3,4,3,4,6,6,10,10,12,12,8,8)
po<-paste(po1,po2)
av<-c("4 10","3 8")
hav<-c("4 12","3 12","6 10","11 10","6 8", "11 8")
fc<-"8 10"
hfc<-c("10 12","8 12")
gpgc<-c("1 10","2 10","1 8","2 8","2 12", "5 12")
fs<-"3 4"
hs<-c("3 6", "4 6")
psel<-is.element(prs,po)
avsel<-is.element(prs,av)
havsel<-is.element(prs,hav)
fcsel<-is.element(prs,fc)
hfcsel<-is.element(prs,hfc)
gpsel<-is.element(prs,gpgc)
fssel<-is.element(prs,fs)
hssel<-is.element(prs,hs)

relprs$relation[psel]<-"PO"
relprs$relation[avsel]<-"Av"
relprs$relation[havsel]<-"HAv"
relprs$relation[fcsel]<-"FC"
relprs$relation[hfcsel]<-"HFC"
relprs$relation[gpsel]<-"GpGc"
relprs$relation[fssel]<-"FS"
relprs$relation[hssel]<-"HS"

relprs$kinship[psel]<-.25
relprs$kinship[avsel]<-.125
relprs$kinship[havsel]<-.0625
relprs$kinship[fcsel]<-.0625
relprs$kinship[hfcsel]<-.03125
relprs$kinship[gpsel]<-.125
relprs$kinship[fssel]<-.25
relprs$kinship[hssel]<-.125

relprs$family<-1

expected.result<-list(inbred.fam=NULL,inbred.KC=NULL,relativeprs=relprs)
result<-pedigreePairwiseRelatedness(samp)
checkEquals(expected.result,result)
########################
}
