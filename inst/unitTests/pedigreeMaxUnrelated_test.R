test_pedigreeMaxUnrelated <- function() {

# use checkTrue instead of checkEquals since order of results can be different
# checkEquals has problems if row names are different

########## Example 0 ###############
# 0a
# check more than one family
family <- c(1,1,1,1,2,2,2,2,2)
individ <- c(2,1,3,4,"A5","A6","A7","A8","A9")
mother <- c(3,3,0,0,0,0,"A5","A5",0)
father <- c(4,4,0,0,0,0,"A6","A9",0)
sex <- c("F","M","F","M","F","M","M","M","M")
samp <- data.frame(family, individ, mother, father, sex, stringsAsFactors=FALSE)
samp$selset<-1
samp$selset[samp$individ %in% c("A5",4)]<-0
# select founders 3,5 from family 1; A9,A6 from family 2 
fm<-c(1,2,2)
ex<-c(3,"A6","A9")
expected.result<-data.frame("family"=fm,"Individ"=ex,stringsAsFactors=FALSE)
result<-pedigreeMaxUnrelated(samp)
result2<-result[order(result$family,result$Individ),]
checkTrue(allequal(expected.result,result2))

# 0b: preference
family <- c(1,1,1,1,2,2,2,2,2)
individ <- c(2,1,3,4,"A5","A6","A7","A8","A9")
mother <- c(3,3,0,0,0,0,"A5","A5",0)
father <- c(4,4,0,0,0,0,"A6","A9",0)
sex <- c("F","M","F","M","F","M","M","M","M")
samp <- data.frame(family, individ, mother, father, sex, stringsAsFactors=FALSE)
samp$selset<-1
samp$selset[samp$individ %in% c("A5",4)]<-0
samp$pref<-2
samp$pref[is.element(samp$individ,c("A8","A7"))]<-1
# NOTE there is now no pref for family 1 (including no pref for founders) so 
# will select one unrelated from family 1: selects 2 since it is first listed in samp   

fm<-c(1,2,2)
ex<-c(2,"A6","A8")
expected.result<-data.frame("family"=fm,"Individ"=ex,stringsAsFactors=FALSE)
result<-pedigreeMaxUnrelated(samp,pref="pref")
result2<-result[order(result$family,result$Individ),]
checkTrue(allequal(expected.result,result2))

# 0c: preference change
family <- c(1,1,1,1,2,2,2,2,2)
individ <- c(2,1,3,4,"A5","A6","A7","A8","A9")
mother <- c(3,3,0,0,0,0,"A5","A5",0)
father <- c(4,4,0,0,0,0,"A6","A9",0)
sex <- c("F","M","F","M","F","M","M","M","M")
samp <- data.frame(family, individ, mother, father, sex, stringsAsFactors=FALSE)
samp$selset<-1
samp$selset[samp$individ %in% c("A5",4)]<-0
samp$pref<-2
samp$pref[is.element(samp$individ,"A7")]<-1   # changed
sel<-samp$family==1 & samp$mother==0 & samp$father==0 # prefer founders in family 1
samp$pref[sel]<-1
# will now select founder 3 from family 1 and will select pair A9, A7 from family 2
fm<-c(1,2,2)
ex<-c(3,"A7","A9")
expected.result<-data.frame("family"=fm,"Individ"=ex,stringsAsFactors=FALSE)
result<-pedigreeMaxUnrelated(samp,pref="pref")
result2<-result[order(result$family,result$Individ),]
checkTrue(allequal(expected.result,result2))

###########################

########### Example 1 #############
# 1a
family<-rep("A",8)
individ<-c("a","b","c","d","e","f","g","h")
mother<-c(0,"a","b",0,"f",0,0,"f")
father<-c(0,"d","e",0,"g",0,0,"g")
sex<-c(rep("F",3),"M","M","F","M","F")
samp <- data.frame(family, individ, mother, father, sex, stringsAsFactors=FALSE)
samp$selset<-1  # all genotyped - expect will choose founders

## preference default (i.e. choose founders first if possible)
fam<-rep("A",4)
fndr<-c("a","d","f","g")
expected.result<-data.frame("family"=fam,"Individ"=fndr,stringsAsFactors=FALSE)
result<-pedigreeMaxUnrelated(samp)
result2<-result[order(result$family,result$Individ),]
checkTrue(allequal(expected.result,result2))
#---------------   

# 1b
family<-rep("A",8)
individ<-c("a","b","c","d","e","f","g","h")
mother<-c(0,"a","b",0,"f",0,0,"f")
father<-c(0,"d","e",0,"g",0,0,"g")
sex<-c(rep("F",3),"M","M","F","M","F")
samp <- data.frame(family, individ, mother, father, sex, stringsAsFactors=FALSE)
samp$selset<-1
sel<-is.element(samp$individ,c("a","f","g"))
samp$selset[sel]<-0  #only one founder 'd' left as genotyped; default pref of founders

## should select 'd' and then first entry (as samp is written) not related - would be e
ex<-c("d","e")
expected.result<-data.frame("family"=rep("A",2),"Individ"=ex,stringsAsFactors=FALSE)
result<-pedigreeMaxUnrelated(samp)
result2<-result[order(result$family,result$Individ),]
checkTrue(allequal(expected.result,result2))
#----------------------

# 1c 
## preference
family<-rep("A",8)
individ<-c("a","b","c","d","e","f","g","h")
mother<-c(0,"a","b",0,"f",0,0,"f")
father<-c(0,"d","e",0,"g",0,0,"g")
sex<-c(rep("F",3),"M","M","F","M","F")
samp <- data.frame(family, individ, mother, father, sex, stringsAsFactors=FALSE)
samp$selset<-1
sel<-is.element(samp$individ,c("a","f","g"))
samp$selset[sel]<-0 
samp$pref<-2
sel2<-samp$individ %in% c("c","h") # preferred
## should select h (c related to everyone) and 'first' (as samp is written) entry not related to h which would be b
samp$pref[sel2]<-1
ex<-c("b","h")
expected.result<-data.frame("family"=rep("A",2),"Individ"=ex,stringsAsFactors=FALSE)
result<-pedigreeMaxUnrelated(samp,pref="pref")
result2<-result[order(result$family,result$Individ),]
checkTrue(allequal(expected.result,result2))
#----------------

## layered preference: secondary founders
family<-rep("A",8)
individ<-c("a","b","c","d","e","f","g","h")
mother<-c(0,"a","b",0,"f",0,0,"f")
father<-c(0,"d","e",0,"g",0,0,"g")
sex<-c(rep("F",3),"M","M","F","M","F")
samp <- data.frame(family, individ, mother, father, sex, stringsAsFactors=FALSE)
samp$selset<-1
sel<-is.element(samp$individ,c("a","f","g"))
samp$selset[sel]<-0 
samp$pref<-3
sel1<-samp$individ %in% c("c","h") # first pref
sel2<-samp$mother==0 & samp$father==0
samp$pref[sel1]<-1
samp$pref[sel2]<-2
## should select 'h' and then founder unrelated 'd'
#Note that other preferred 'c' is related to everyone so not chosen
ex<-c("d","h")
expected.result<-data.frame("family"=rep("A",2),"Individ"=ex,stringsAsFactors=FALSE)
result<-pedigreeMaxUnrelated(samp,pref="pref")
result2<-result[order(result$family,result$Individ),]
checkTrue(allequal(expected.result,result2))
#-----------------------------
###############################################

######## Example 2 #############
# note has inbreeding

#2a
family <- rep(2,7)
individ <- paste("I",c(1,2,3,4,5,6,7),sep="")
mother <- c(0,0,0,"I1","I1","I3","I5")
father <- c(0,0,0,"I2","I2","I4","I4")
sex <- c("F","M","F","M","F","F","F")
samp <- data.frame(family, individ, mother, father, sex, stringsAsFactors=FALSE)
samp$selset<-1  #all genotyped
## expect founders
ex<-c("I1","I2","I3")
expected.result<-data.frame("family"=rep(2,3),"Individ"=ex,stringsAsFactors=FALSE)
result<-pedigreeMaxUnrelated(samp)
result2<-result[order(result$family,result$Individ),]
checkTrue(allequal(expected.result,result2))
#-------------------------------

# 2b
family <- rep(2,7)
individ <- paste("I",c(1,2,3,4,5,6,7),sep="")
mother <- c(0,0,0,"I1","I1","I3","I5")
father <- c(0,0,0,"I2","I2","I4","I4")
sex <- c("F","M","F","M","F","F","F")
samp <- data.frame(family, individ, mother, father, sex, stringsAsFactors=FALSE)
samp$selset<-1
samp$selset[samp$individ=="I1"]<-0
ex<-c("I2","I3")
expected.result<-data.frame("family"=rep(2,2),"Individ"=ex,stringsAsFactors=FALSE)
result<-pedigreeMaxUnrelated(samp)
result2<-result[order(result$family,result$Individ),]
checkTrue(allequal(expected.result,result2))
#-------------------------------

# 2c
## first pref and then founder pref
family <- rep(2,7)
individ <- paste("I",c(1,2,3,4,5,6,7),sep="")
mother <- c(0,0,0,"I1","I1","I3","I5")
father <- c(0,0,0,"I2","I2","I4","I4")
sex <- c("F","M","F","M","F","F","F")
samp <- data.frame(family, individ, mother, father, sex, stringsAsFactors=FALSE)
samp$selset<-1
samp$selset[samp$individ=="I1"]<-0
samp$pref<-3
samp$pref[samp$mother==0 & samp$father==0]<-2
samp$pref[samp$individ=="I4"]<-1
## Note here that I3 is individual with most unrelated (3) so it gets chosen
#   then I2, I4, and I5 unrelated to only 1 other
#   preference given to 4 which gets chosen
ex<-c("I3","I4")
expected.result<-data.frame("family"=rep(2,2),"Individ"=ex,stringsAsFactors=FALSE)
result<-pedigreeMaxUnrelated(samp,pref="pref")
result2<-result[order(result$family,result$Individ),]
checkTrue(allequal(expected.result,result2))
#---------------------

# 2d
family <- rep(2,7)
individ <- paste("I",c(1,2,3,4,5,6,7),sep="")
mother <- c(0,0,0,"I1","I1","I3","I5")
father <- c(0,0,0,"I2","I2","I4","I4")
sex <- c("F","M","F","M","F","F","F")
samp <- data.frame(family, individ, mother, father, sex, stringsAsFactors=FALSE)
samp$selset<-1
samp$selset[samp$mother==0 & samp$father==0]<-0  # no founders genotyped; everyone else is related 
samp$pref<-2
samp$pref[samp$individ == "I6"]<-1
# all related so should choose one person - the preferred
ex<-"I6"
expected.result<-data.frame("family"=2,"Individ"=ex,stringsAsFactors=FALSE)
result<-pedigreeMaxUnrelated(samp,pref="pref")
result2<-result[order(result$family,result$Individ),]
checkTrue(allequal(expected.result,result2))
#--------------
###########################

######### Example 3 #############

# 3a
family<-rep(1,10)
individ<-paste("I",1:10,sep="")
mother<-c(0,0,"I1","I1",0,0,"I6","I6","I4","I7")
father<-c(0,0,"I2","I2",0,0,"I5","I5","I8","I3")
sex<-c("F","M","M","F","M","F","F","M","F","F")
samp<-data.frame(family,individ,mother,father,sex,stringsAsFactors=F)
samp$selset<-1
# should choose founders
ex<-c("I1","I2","I5","I6")
expected.result<-data.frame("family"=rep(1,4),"Individ"=ex,stringsAsFactors=FALSE)
result<-pedigreeMaxUnrelated(samp)
result2<-result[order(result$family,result$Individ),]
checkTrue(allequal(expected.result,result2))
#-----------------

# 3b
family<-rep(1,10)
individ<-paste("I",1:10,sep="")
mother<-c(0,0,"I1","I1",0,0,"I6","I6","I4","I7")
father<-c(0,0,"I2","I2",0,0,"I5","I5","I8","I3")
sex<-c("F","M","M","F","M","F","F","M","F","F")
samp<-data.frame(family,individ,mother,father,sex,stringsAsFactors=F)
samp$selset<-1
samp$selset[is.element(samp$individ,c("I1","I2","I5","I6"))]<-0  #founders not genotyped
# should choose 2 unrelated in order (as samp is written) I3, I7
ex<-c("I3","I7")
expected.result<-data.frame("family"=rep(1,2),"Individ"=ex,stringsAsFactors=FALSE)
result<-pedigreeMaxUnrelated(samp)
result2<-result[order(result$family,result$Individ),]
checkTrue(allequal(expected.result,result2))
#-------------


# 3c
family<-rep(1,10)
individ<-paste("I",1:10,sep="")
mother<-c(0,0,"I1","I1",0,0,"I6","I6","I4","I7")
father<-c(0,0,"I2","I2",0,0,"I5","I5","I8","I3")
sex<-c("F","M","M","F","M","F","F","M","F","F")
samp<-data.frame(family,individ,mother,father,sex,stringsAsFactors=F)
samp$selset<-1
samp$selset[is.element(samp$individ,c("I1","I2","I5","I6"))]<-0  
samp$pref<-2
samp$pref[is.element(samp$individ,c("I4","I8"))]<-1
# should choose 2 unrelated which will be the preferred in this case
ex<-c("I4","I8")
expected.result<-data.frame("family"=rep(1,2),"Individ"=ex,stringsAsFactors=FALSE)
result<-pedigreeMaxUnrelated(samp,pref="pref")
result2<-result[order(result$family,result$Individ),]
checkTrue(allequal(expected.result,result2))
#----------------

# 3d and 3e
family<-rep(1,10)
individ<-paste("I",1:10,sep="")
mother<-c(0,0,"I1","I1",0,0,"I6","I6","I4","I7")
father<-c(0,0,"I2","I2",0,0,"I5","I5","I8","I3")
sex<-c("F","M","M","F","M","F","F","M","F","F")
samp<-data.frame(family,individ,mother,father,sex,stringsAsFactors=F)
samp$selset<-1
samp$selset[is.element(samp$individ,c("I1","I2","I5","I6"))]<-0  
samp$pref<-2
samp$pref[is.element(samp$individ,c("I4","I8"))]<-1
samp$selset[samp$individ=="I1"]<-1  # so one founder is now genotyped
# default would choose founder - the first one in order and then another unrelated in order
ex<-c("I1","I7")
expected.result<-data.frame("family"=rep(1,2),"Individ"=ex,stringsAsFactors=FALSE)
result<-pedigreeMaxUnrelated(samp)
result2<-result[order(result$family,result$Individ),]
checkTrue(allequal(expected.result,result2))

# using the preference will choose 2 unrelated with preference above other choices
ex<-c("I4","I8")
expected.result<-data.frame("family"=rep(1,2),"Individ"=ex,stringsAsFactors=FALSE)
result<-pedigreeMaxUnrelated(samp,pref="pref")
result2<-result[order(result$family,result$Individ),]
checkTrue(allequal(expected.result,result2))
##############################

############# Example 4 ##############
# 4a
family<-rep(1,12)
individ<-1:12
mother<-c(0,0,1,1,0,5,0,4,0,9,0,6)
father<-c(0,0,2,2,0,2,0,7,0,3,0,11)
sex<-c("F","M","M","F","F","F","M","F","F","F","M","F")
samp<-data.frame(family,individ,mother,father,sex,stringsAsFactors=FALSE)
samp$selset<-1
samp$selset[samp$individ %in% c(7,9,11)]<-0
## default should give the three other founders
ex<-c(1,2,5)
expected.result<-data.frame("family"=rep(1,3),"Individ"=ex,stringsAsFactors=FALSE)
result<-pedigreeMaxUnrelated(samp)
result2<-result[order(result$family,result$Individ),]
checkTrue(allequal(expected.result,result2))
#-----------------

# 4b
family<-rep(1,12)
individ<-1:12
mother<-c(0,0,1,1,0,5,0,4,0,9,0,6)
father<-c(0,0,2,2,0,2,0,7,0,3,0,11)
sex<-c("F","M","M","F","F","F","M","F","F","F","M","F")
samp<-data.frame(family,individ,mother,father,sex,stringsAsFactors=FALSE)
samp$selset<-1
samp$selset[samp$individ %in% c(7,1,2,5,9,11)]<-0  # no founders genotyped
# max set contains only 1 person, here we get 3 (first person listed in samp genotyped)
ex<-3
expected.result<-data.frame("family"=1,"Individ"=ex,stringsAsFactors=FALSE)
result<-pedigreeMaxUnrelated(samp)
result2<-result[order(result$family,result$Individ),]
checkTrue(allequal(expected.result,result2))
#------------------

# 4c
family<-rep(1,12)
individ<-1:12
mother<-c(0,0,1,1,0,5,0,4,0,9,0,6)
father<-c(0,0,2,2,0,2,0,7,0,3,0,11)
sex<-c("F","M","M","F","F","F","M","F","F","F","M","F")
samp<-data.frame(family,individ,mother,father,sex,stringsAsFactors=FALSE)
samp$selset<-1
samp$selset[samp$individ==5]<-0
samp$pref<-2
samp$pref[samp$individ==6]<-1

# without pref - other founders
ex<-c(1,2,7,9,11)
expected.result<-data.frame("family"=rep(1,5),"Individ"=ex,stringsAsFactors=FALSE)
result<-pedigreeMaxUnrelated(samp)
result2<-result[order(result$family,result$Individ),]
checkTrue(allequal(expected.result,result2))

# with pref: choose 6 over 2
ex<-c(1,6,7,9,11)
expected.result<-data.frame("family"=rep(1,5),"Individ"=ex,stringsAsFactors=FALSE)
result<-pedigreeMaxUnrelated(samp,pref="pref")
result2<-result[order(result$family,result$Individ),]
checkTrue(allequal(expected.result,result2))
#---------------

# 4d
family<-rep(1,12)
individ<-1:12
mother<-c(0,0,1,1,0,5,0,4,0,9,0,6)
father<-c(0,0,2,2,0,2,0,7,0,3,0,11)
sex<-c("F","M","M","F","F","F","M","F","F","F","M","F")
samp<-data.frame(family,individ,mother,father,sex,stringsAsFactors=FALSE)
samp$selset<-1
samp$selset[c(1,2 ,7)]<-0
samp$pref<-2
samp$pref[12]<-1
# here using the pref will not choose 12 since it is more highly related than 11
ex<-c(3,5,9,11)
expected.result<-data.frame("family"=rep(1,4),"Individ"=ex,stringsAsFactors=FALSE)
result<-pedigreeMaxUnrelated(samp,pref="pref")
result2<-result[order(result$family,result$Individ),]
checkTrue(allequal(expected.result,result2))
#-----------------------

#4e
family<-rep(1,12)
individ<-1:12
mother<-c(0,0,1,1,0,5,0,4,0,9,0,6)
father<-c(0,0,2,2,0,2,0,7,0,3,0,11)
sex<-c("F","M","M","F","F","F","M","F","F","F","M","F")
samp<-data.frame(family,individ,mother,father,sex,stringsAsFactors=FALSE)
samp$selset<-1
samp$selset[c(1,2,5,7,11)]<-0 
# remove more founders; note there will be max set of two 
samp$pref<-2
samp$pref[12]<-1
ex<-c(9,12)
expected.result<-data.frame("family"=rep(1,2),"Individ"=ex,stringsAsFactors=FALSE)
result<-pedigreeMaxUnrelated(samp,pref="pref")
result2<-result[order(result$family,result$Individ),]
checkTrue(allequal(expected.result,result2))

# if don't use pref
ex<-c(3,9)
expected.result<-data.frame("family"=rep(1,2),"Individ"=ex,stringsAsFactors=FALSE)
result<-pedigreeMaxUnrelated(samp)
result2<-result[order(result$family,result$Individ),]
checkTrue(allequal(expected.result,result2))

}
