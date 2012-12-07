test_pedigreeCheck <- function() {
#test data for initial basic cleaning
family<-c("a","a","a","b","b","c","")
individ<-c("A","B","C","A","B",0,"")
mother<-c("B","C",0,0,0,NA,0)
father<-c("C","D",0,0,"",0,"D")
sex<-c("F","2","M","F","F","M","F")
samp<-data.frame(family, individ, mother,father,sex,stringsAsFactors=F)
expected.result<-list(family.missing.rows=7,individ.missing_or_0.rows = c(6,7),
    father.missing.rows=5,mother.missing.rows = 6, sexcode.error.rows=2)
result<-pedigreeCheck(samp)
checkEquals(expected.result, result) 
###################################
family<-c(rep(1,6),"2a","2a","2a",3,NA)
individ<-c(1,2,0,1,"",3,"0",2,3,1,1)
mother<-c(0,1,0,0,NA,0,2,0,0,1,0)
father<-c(0,3,0,0,3,1,3,"","0",2,2)
sex<-c("F","2","M","F","F","M","F","F","M","M","M")
samp<-data.frame(family, individ, mother,father,sex,stringsAsFactors=FALSE)
pid<-data.frame("row.num"=10L,"family"="3","no_individ_entry"="father","parentID"="2",stringsAsFactors=FALSE)
ofam<-data.frame("family"="3","founder"=FALSE,stringsAsFactors = FALSE)
expected.result<-list(family.missing.rows=11,individ.missing_or_0.rows =c(3,5,7), father.missing.rows=8,mother.missing.rows=5,sexcode.error.rows=2, parent.no.individ.entry=pid,one.person.fams=ofam)
result<-pedigreeCheck(samp)
checkEquals(expected.result,result)

## NOTE: family 1 has other issues but was not checked since had other earlier errors 
#########################################

family<-c("a","a","a","b","b","b","c")
individ<-c("A","B","C","A","B","C","B")
mother<-c("B",0,0,"B",0,0,0)
father<-c("C",0,0,"C",0,0,0)
sex<-c("F","F","M","F","F","M","F")
samp<-data.frame(family, individ, mother,father,sex,stringsAsFactors=FALSE)
pedigreeCheck(samp)
## should be ok for basic clean but has one person fam
ofam<-data.frame("family"="c","founder"=TRUE, stringsAsFactors=FALSE)
expected.result<-list(one.person.fams=ofam)
result<-pedigreeCheck(samp) 
checkEquals(expected.result,result)
#-----------------------------

#Test data unknown parents
family<-c(rep(1,4),2,3,3,4,4,4,5)
individ<-c(1,"a2",3,4,1,1,"a2",1,"a2","a3",1)
mother<-c("a2",0,0,"a2",2,"a2",0,"a2",0,0,3)
father<-c(3,0,0,0,3,3,0,"a3",0,0,0)
sex<-c("F","M","M","M","M","F","M","M","F","M","F")
samp<-data.frame(family, individ, mother,father,sex,stringsAsFactors=F)
   ## note here that samp$family is 'numeric' so output will be 'numeric'
pid<-data.frame("row.num"=c(5L,6L,11L),"family"=c(2,3,5),"no_individ_entry"=c("both","father","mother"),"parentID"=c("2;3","3","3"),stringsAsFactors=FALSE)
ukp<-data.frame("row.num"=c(4L,11L),"family"=c(1,5),stringsAsFactors=FALSE)
ofam<-data.frame("family"=c(2,5),founder=c(FALSE,FALSE),stringsAsFactors=FALSE)
expected.result<-list(parent.no.individ.entry=pid,one.person.fams=ofam,unknown.parent.rows=ukp)
result<-pedigreeCheck(samp)
checkEquals(expected.result,result)
##########################
family<-rep(1L,10) 
individ<-1:10
class(individ)<-"integer"
mother<-c(2,0,2,0,0,7,0,2,12,7)
father<-c(4,0,5,0,0,5,0,4,11,0)
sex<-c("F","F","M","M","M","F","F","F","M","F")
samp<-data.frame(family,individ,mother,father,sex,stringsAsFactors=FALSE)
pid<-data.frame("row.num"=9L,"family"=1L,"no_individ_entry"="both","parentID"="12;11",stringsAsFactors=FALSE)
ukp<-data.frame("row.num"=10L,"family"=1L,stringsAsFactors=FALSE)
expected.result<-list(parent.no.individ.entry=pid,unknown.parent.rows=ukp)
result<-pedigreeCheck(samp)
checkEquals(expected.result,result)
##################################
## 'duplicates', 'both.mother.father', 'parent.no.individ.entry'
family<-c("b","b","b","b","c","c",rep("d",5))
individ<-c("A","B","C","A","B","B",1:5)
mother<-c("B",0,0,"D",0,0,0,0,1,2,1)
father<-c("C",0,0,"C",0,0,0,0,2,1,2)
sex<-c("F","F","M","M","F","F","F","M","F","F","M")
samp<-data.frame(family, individ, mother,father,sex,stringsAsFactors=FALSE)
pid<-data.frame("row.num"=4L,"family"="b","no_individ_entry"="mother","parentID"="D",stringsAsFactors=FALSE)
dp<-data.frame("family"=c("b","c"), "individ"=c("A","B"),"copies"=c(2L,2L),"match"=c(FALSE,TRUE),stringsAsFactors=FALSE)
bmf<-data.frame("family"=rep("d",2),"parentID"=c("1","2"),"mother.row"=c("9;11","10"),"father.row"=c("10","9;11"),stringsAsFactors=FALSE)
expected.result<-list("both.mother.father"=bmf,"duplicates"=dp,"parent.no.individ.entry"=pid)
result<-pedigreeCheck(samp)
checkEquals(expected.result,result)
# there are other problems (such as mismatch.sex but not investigated 
#     directly because already had both.mother.father inconsistency)
###################################
#Test data for impossible related, sex mismatch, others
family<-c(1,1,1,2,2,2,3,4,4,4,5,5,5,5,6,6,6)
individ<-c(1,2,3,1,2,3,1,1,3,2,1,2,3,4,1,2,3)
mother<-c(2,0,1,2,1,0,1,2,0,2,2,4,0,0,2,1,0)
father<-c(3,0,3,0,3,0,2,3,1,0,3,1,0,0,3,3,0)
sex<-c("F","F","M","F","F","M","F","F","F","F","M","F","M","F","F","M","F")
samp<-data.frame(family, individ,mother,father,sex,stringsAsFactors=FALSE)
pid<-data.frame("row.num"=7,"family"=3,"no_individ_entry"="father","parentID"=2,stringsAsFactors=FALSE)
ukp<-data.frame("row.num"=c(4L,9L,10L),"family"=c(2,4,4),stringsAsFactors=FALSE)
ofam<-data.frame("family"=3,"founder"=FALSE)
msex<-data.frame("family"=c(6,6),"individ"=c(2,3),stringsAsFactors=FALSE)
ir<-vector("list",3)
names(ir)<-c(1,5,6)
ir[[1]]<-c(1L,3L)
ir[[2]]<-c(11L,12L)
ir[[3]]<-c(15L,16L)
expected.result<-list(parent.no.individ.entry=pid,one.person.fams=ofam,unknown.parent.rows=ukp,
   mismatch.sex=msex,impossible.related.rows=ir)
result<-pedigreeCheck(samp)
checkEquals(expected.result,result)

## NOTE: mismatch.sex and impossible.related.rows are investigated
#     only if no other errors
#    Hence even though family 2 has impossible relationship,
#        it had previous error of unknown parent so impossible relationship not investigated
#     similarly for one person family 3 - 'one person' designation is previous error
#    similarly for family 4 - had mismatch sex and impossible relationship but also unknown parent
###############################

 family<-c(rep(1,5),rep(4,3))
 individ<-c("R1","R2","R3","R4","R5","a1","a2","a3")
 mother<-c("R2",0,"R1",0,"R5","a2",0,0)
 father<-c("R3",0,"R4",0,"R4","a3",0,0)
 sex<-c("F","F","M","M","F","F","F","F")
 samp<-data.frame(family, individ,mother,father,sex,stringsAsFactors=FALSE)
 msex<-data.frame("family" =4, "individ" = "a3",stringsAsFactors=FALSE)
 ir<-vector("list",1)
 names(ir)<-"1"
 ir[[1]]<-c(1,3,5)
 expected.result<-list(mismatch.sex=msex,impossible.related.rows=ir)
 result<-pedigreeCheck(samp)
checkEquals(expected.result,result)
 # note that for family 1, impossible.related.rows actually 'contain' two impossible relationships
###############################
##Test data for checking relatedness
family<-rep(1,13)
individ<-c(1:13)
mother<-c(rep(0,7),2,2,4,4,6,9)
father<-c(rep(0,7),1,1,3,5,7,10)
sex<-c("M","F","M","F","M","F","M","M","F","M","F","F","F")
samp<-data.frame(family,individ,mother,father,sex,stringsAsFactors=FALSE)
fm<-rep(1,13)
sfm<-c(rep(1,10),rep(2,3))
id<-c(1:5,8:11,13,6,7,12)
sub<-data.frame("family"=fm,"subfamily"=sfm,"individ"=id,stringsAsFactors=FALSE)
#expected.result<-list(subfamilies.ident=sub)
result<-pedigreeCheck(samp)
checkTrue(allequal(sub,result[[1]]))
##################################
family<-rep(1,12)
individ<-1:12
mother<-c(0,0,2,2,0,0,5,0,7,0,0,10)
father<-c(0,0,1,1,0,0,6,0,8,0,0,11)
sex<-c("M",rep("F",4),"M","F","M","M","F","M","M")
samp<-data.frame(family,individ,mother,father,sex,stringsAsFactors=FALSE)
fm<-family
sfm<-c(rep(1,4),rep(2,5),rep(3,3))
id<-c(1:4,5:9,10:12)
sub<-data.frame("family"=fm,"subfamily"=sfm,"individ"=id,stringsAsFactors=FALSE)
#expected.result<-list(subfamilies.ident=sub)
result<-pedigreeCheck(samp)
checkTrue(allequal(sub,result[[1]]))
###################################
# connect two of the subfamilies
family<-rep(1,13)
individ<-1:13
mother<-c(0,0,2,2,0,0,5,0,7,0,0,10,4)
father<-c(0,0,1,1,0,0,6,0,8,0,0,11,9)
sex<-c("M",rep("F",4),"M","F","M","M","F","M","M","F")
samp<-data.frame(family,individ,mother,father,sex,stringsAsFactors=FALSE)
fm<-family
sfm<-c(rep(1,10),rep(2,3))
id<-c(1:9,13,10:12)
sub<-data.frame("family"=fm,"subfamily"=sfm,"individ"=id,stringsAsFactors=FALSE)
result<-pedigreeCheck(samp)
checkTrue(allequal(sub,result[[1]]))
#################################
## connect all subfamilies
family<-rep(1,14)
individ<-1:14
mother<-c(0,0,2,2,0,0,5,0,7,0,0,10,4,12)
father<-c(0,0,1,1,0,0,6,0,8,0,0,11,9,9)
sex<-c("M",rep("F",4),"M","F","M","M","F","M","F","F","F")
samp<-data.frame(family,individ,mother,father,sex,stringsAsFactors=FALSE)
expected.result<-NULL
result<-pedigreeCheck(samp) # NULL all subfams now connected
checkEquals(expected.result,result)

#############

}
