### LOH raw determination ###
########## SUB-functions ##########

### Find RUNS ##########
LOHruns<-function(index,geno,eleft,eright,run.size,inter.size) {
## index are eligible snp indices for a sample/chrom
## geno is genotypes corresponding to eligible snp's
## eleft and eright:endpoints of interval of interest given as snp indices
## run.size is minimum number of wanted markers to declare a run
## inter.size is number of unwanted markers allowed to interrupt a run
 
if(eright<eleft) stop("non-existent interval given to LOHruns")
indices<-index>=eleft & index<=eright
if(sum(indices)<run.size){out<-NULL; return(out)}

index1<-index[indices]
geno1<-geno[indices]

w<-rle(as.vector(geno1))
#w<-rle(  ) has w[[1]] = lengths of runs, w[[2]] corresponding values
vals<-w[[2]];lngs<-w[[1]]
r0<-vals

rlen<-length(r0)
if(rlen==1){ if(r0==0){left<-index1[1];right<-index1[length(index1)]
 out<-data.frame(left,right)
 names(out)<-c("left","right")
 return(out)}  else {out<-NULL;return(out)}  }

##establish initial positions of alternating runs
in.pos<-c(1,rep(NA,rlen-1))
for(i in 2:rlen) in.pos[i]<-in.pos[i-1]+lngs[i-1]

##merging homoz intervals if separated by inter.size no. of hets
# assuming sum of lengths of homo intervals on either side meets run.size criterion

tpos<-which(r0==1&lngs<=inter.size) #identify small runs of hets
smf<-which(r0==0 & lngs<run.size) #identify 'small' runs of homos
  if(length(tpos)!=0){
if(tpos[1]==1) {
   if(lngs[2]>=run.size){r0[1]<-0}
   tpos<-tpos[-1]}
if(length(tpos)!=0){
if(tpos[length(tpos)]==rlen) {tpos<-tpos[-length(tpos)]; if(lngs[rlen-1]>=run.size) {r0[rlen]<-0}}
 if(length(tpos)!=0){
for(k in tpos){if((lngs[k-1]+lngs[k+1])>=run.size) {r0[k]<-0}}
  }}  }

##want smaller runs of 0s to become runs of 'hets' but not if they
 #are part of a combined run
run2<-rle(r0)
vals2<-run2[[2]];lngs2<-run2[[1]]
w0<-which(vals2==0 & lngs2>1)
if(length(w0)==0){r0[smf]<-1} else {
index.set<-NULL
for(j in 1:length(w0)){ if(w0[j]==1){start<-1} else {start<-sum(lngs2[1:(w0[j]-1)])+1}
end<-sum(lngs2[1:w0[j]])
index.set<-c(index.set,start:end)}
smf.use<-setdiff(smf,intersect(smf,index.set))
r0[smf.use]<-1}

##after merging some of the initial runs, we get modified listing of r0 vals
## look for runs here; e.g. if run of two 0s that means that initially
#we had a run of homo with small run of hets - rle length of 2 then indicates
#putting those two original runs together as one run of homo 
new.rle<-rle(r0)
nvals<-new.rle[[2]]
nlens<-new.rle[[1]]  ##indicates how many of original runs to put together
if(length(nvals)==1){
  if(nvals!=0){out<-NULL;return(out)} else {
left<-index1[1];right<-index1[length(index1)]
 out<-data.frame(left,right)
 names(out)<-c("left","right")
 return(out)} }

newt<-which(nvals==0) #newt could be empty if originally there were no long homo runs
if(length(newt)==0){out<-NULL;return(out)}
left<-NULL
right<-NULL
### if newt indicates runs of 0s, change initial and end positions of 
#runs of homozy accordingly
if(newt[1]==1){left<-c(left,1);right<-c(right,in.pos[nlens[1]+1]-1);newt<-newt[-1]}
for(k in newt){ind<-sum(nlens[1:(k-1)]);left<-c(left,in.pos[ind+1])
kl<-length(newt);kk<-newt[kl]
if((ind+1+nlens[kk])<=length(in.pos)){
right<-c(right,in.pos[ind+1+nlens[k]]-1)} else {right<-c(right,length(index1))}}
##right and left positions are indices of geno1, indices of index1
if(length(right)==0|length(left)==0){out<-NULL;return(out)}
out<-data.frame(index1[left],index1[right])
 names(out)<-c("left","right")
return(out)
 }#end function
################

########## Find BASE info #############
LOHbase<-function(anoms,index,lrr,min.lrr.num) {
## anoms are anomalies found (BAF as well as potential LOH)  
# for a given sample/chrom with left, right given as snp indices
## index,lrr,pos are for a given sample/chrom
#min.lrr.num - minimum number of lrr pts to be considered/processed

if(dim(anoms)[1]==0){newright<-index[length(index)]
  newleft<-index[1]} else {anoms<-anoms[order(anoms$left),]
   newright<-c(anoms$left,index[length(index)])
   newleft<-c(index[1],anoms$right)
 }

lrr.pts<-NULL
for(k in 1:length(newright)){ 

 int.ind<-index>=newleft[k]& index<=newright[k] 
 if(sum(int.ind)<min.lrr.num) next
 
 y<-lrr[int.ind]
 lrr.pts<-c(lrr.pts,y)
   }
if(length(lrr.pts)<min.lrr.num) {out<-data.frame(NA,NA,NA,NA)} else{

 chmad<-mad(lrr.pts,na.rm=TRUE)
 chmed<-median(lrr.pts,na.rm=TRUE)
 chmn<-mean(lrr.pts,na.rm=TRUE)
 chsd<-sd(lrr.pts,na.rm=TRUE)
 out<-data.frame(chmad,chmed,chmn,chsd)
}
names(out)<-c("chrom.nonanom.mad","chrom.nonanom.median","chrom.nonanom.mean","chrom.nonanom.sd")
return(out)  }
################################## 
########### MAIN FUNCTION ##############

LOHfind<-function(snum,ch,geno,index,lrr,chr,segs,
  ansch,run.size=50,inter.size=4,min.lrr.num=20 ) {

#snum - current sample number
#ch - current chromosome 
#geno - genotypes for selected snps for that chrom
#index - snp indices for selected snps for that chrom
#lrr - LogRRatio values for selected snps for that chrom
#chr - vector corresponding to given ch (same length as index, lrr)
#segs - segments for given snum/ch from DNAcopy
#ansch - previously found anomalies (e.g. BAF) to remove before finding runs
#run.size - number of markers needed to declare a run
#inter.size - number of consecutive hets allowed in a run
#min.lrr.num - minimum number of lrr pts to be considered

if(!is.element(class(ansch),"data.frame")) stop(paste("known anom input for sample ",snum," chromosome ",ch," to LOHfind needs to be data.frame",sep=""))
if(!is.element(class(segs),"data.frame")) stop(paste("segs input to LOHfind for sample ",snum," chromosome ",ch," needs to be data.frame",sep=""))
if(any(ansch$left>=ansch$right))stop(paste("some left >= right for known anoms for sample ",snum," chromosome ",ch,sep=""))
###
num.segs<-dim(segs)[1]
annum<-dim(ansch)[1]

####### find RUNS ###########
## find homozygous runs in intervals where there are no known anomalies
## (i.e. take out intervals determined by known anomalies)
RUNS<-NULL
if(annum!=0){
 LF<-c(index[1],ansch$right)
 RT<-c(ansch$left,index[length(index)])} else {
    LF<-index[1];RT<-index[length(index)]    }
    
N<-length(LF)
for(j in 1:N){ eleft<-LF[j];eright<-RT[j]
 outrun<-GWASTools:::LOHruns(index,geno,eleft,eright,run.size,inter.size)  #looking for runs in regions outside anoms
 RUNS<-rbind(RUNS,outrun)  }  #RUNS$right and $left are snp indices
#################### find base info ########
if(is.null(RUNS)) { 
base<-GWASTools:::LOHbase(ansch,index,lrr,min.lrr.num)
base$num.runs<-0
base$num.segs<-num.segs
base$scanID<-snum
base$chrom<-ch

if(!is.na(base$chrom.nonanom.mean) & !is.na(base$chrom.nonanom.sd)) {
    segs$sd.fac<-(segs$seg.mean-base$chrom.nonanom.mean)/base$chrom.nonanom.sd } else {segs$sd.fac<-NA}
rout<-list(RUNS,base,segs)
 names(rout)<-c("RUNS","base.info","segments")
 return(rout)  }
###########
num.runs<-dim(RUNS)[1]
RUNS$scanID<-snum
RUNS$chrom<-ch

## Base line info: chrom.mad, etc are computed for 'non.anom' region####
##  'non.anom' region excludes BAF anomalies and homoz runs (since these are potential anoms)####
anoms1<-ansch[,c("left","right")]
anoms2<-RUNS[,c("left","right")]
anoms<-rbind(anoms1,anoms2)
anoms<-anoms[order(anoms$left),]
base<-GWASTools:::LOHbase(anoms,index,lrr,min.lrr.num)
base$num.runs<-num.runs
base$num.segs<-num.segs
base$scanID<-snum
base$chrom<-ch
if(!is.na(base$chrom.nonanom.mean) & !is.na(base$chrom.nonanom.sd)) {
    segs$sd.fac<-(segs$seg.mean-base$chrom.nonanom.mean)/base$chrom.nonanom.sd } else {segs$sd.fac<-NA}
rout<-list(RUNS,base,segs)
 names(rout)<-c("RUNS","base.info","segments")
 return(rout)  

} #end of function



