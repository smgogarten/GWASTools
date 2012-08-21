############ Find anomalous segments from homozygous runs ##############
#### subfunctions to merge consecutive intervals passing filter #########
mmerge<-function(w,segs){
N<-length(w)
  if(N<=1) {flag<-0;tmp2<-segs[w,c("left","right")]
out<-list(flag,tmp2);
  names(out)<-c("flag","anoms"); return(out)}
wdiff<-sapply(2:N,function(w,i){w[i]-w[i-1]},w=w)
 rn<-rle(wdiff)
vals<-rn[[2]]
lng<-rn[[1]]
wone<-which(vals==1)
if(length(wone)==0) { flag<-0;tmp2<-segs[w,c("left","right")]; out<-list(flag,tmp2);
  names(out)<-c("flag","anoms"); return(out)}

rlen<-length(vals)
if(rlen<2){init.pos<-1} else{
init.pos<-c(1,rep(NA,rlen-1))
for(i in 2:rlen) init.pos[i]<-init.pos[i-1]+lng[i-1]  }

n1<-length(wone)
comb.list<-list()
for(j in 1:n1){
 comb.list[[j]]<-init.pos[wone[j]]:(init.pos[wone[j]]+lng[wone[j]]) 
## comb.list are positions of w
##we now need positions in list of segments
 comb.list[[j]]<-w[comb.list[[j]]]}

mod.ind<-unlist(comb.list)
an.comb<-NULL
for(j in 1:n1) {
  ind<-comb.list[[j]]
  left<-min(segs$left[ind])
  right<-max(segs$right[ind])
  tmp<-data.frame(left,right)
  names(tmp)<-c("left","right")
  an.comb<-rbind(an.comb,tmp)   }
ws<-setdiff(w,mod.ind)
an.nomod<-segs[ws,c("left","right")]
tmp2<-rbind(an.nomod,an.comb)
flag<-1;out<-list(flag,tmp2)
names(out)<-c("flag","anoms")
return(out)   } #end function mmerge
################### main merge ##############
runsMerge<-function(segs,sig) {
## output is new list of runs with merged endpoints ##

 wup<-which(segs$sd.fac>sig)
 wdown<-which(segs$sd.fac < (-sig) )
 del<-union(wup,wdown)
 if(length(del)==0){rest<-segs[,c("left","right")]} else {
    
 rest<-segs[-del,c("left","right")] }

### up ###
tmpup<-NULL
 if(length(wup)!=0){
 resup<-GWASTools:::mmerge(wup,segs)
tmpup<-resup$anoms  } # end if(wup!=0)
## down ## 
tmpdown<-NULL
 if(length(wdown)!=0){
 resdown<-GWASTools:::mmerge(wdown,segs)
  tmpdown<-resdown$anoms
 } # end if(wdown!=0)
#####
flag<-FALSE
if(length(wup)!=0) { if(resup$flag==1) flag<-TRUE}
if(length(wdown)!=0) { if (resdown$flag==1) flag<-TRUE}
 out<-rbind(rest,tmpup,tmpdown)
out<-out[order(out$left),]
rout<-list(out,flag)
names(rout)<-c("newsegs","flag")
return(rout)
  } #end runsMerge 
############################################################
############local mad #########################
LOHlocalMad<-function(select,index,nonanom.index,lrr,length.factor,min.lrr.num){
if(any(select$left>select$right))stop("some left endpts > right endpts in runs input dataframe")
tm5<-NULL
 for(j in 1:dim(select)[1]){
 lt<-select$left[j]; rt<-select$right[j]
 int<-index>=select$left[j] & index<=select$right[j] 
 xx<-lrr[int]
 xxm<-median(xx,na.rm=TRUE)

 zf<-length.factor*sum(int)
  left.bdy<-index[1]
  right.bdy<-index[length(index)]
  nleft<-lt-zf; nright<-rt+zf
  nleft<-max(nleft,left.bdy);nright<-min(nright,right.bdy)
int.test<-index>=nleft & index<=nright
 nonanom<-is.element(index,nonanom.index)
 int5<-int.test&nonanom
  yy<-lrr[int5]
 if(length(yy)< min.lrr.num){tm5<-c(tm5,NA)} else { 
  yym<-median(yy,na.rm=TRUE);yymad<-mad(yy,na.rm=TRUE)
  tm<-(xxm-yym)/yymad
  tm5<-c(tm5,tm)  }

 } #end of loop on j
out<-select
 out$local<-tm5
return(out)
  } #end function

############# MAIN ########################
## find breakpoints in a homozygous run ##
##############

LOHselectAnoms<-function(snum,ch,segs,RUNS,base,index,nonanom.index,lrr,
 homodel.min.num=10,homodel.thresh=10,small.num=20,small.thresh=2.25, 
medium.num=50,medium.thresh=2,long.num=100,long.thresh=1.5, small.na.thresh=2.5,
 length.factor=5,merge.fac=.85,min.lrr.num=20) {

#segs: DNAcopy segments for current sample/chromsome
#
if(!all(is.element(class(RUNS),"data.frame"))) stop("RUNS needs to be a data.frame")
if(!all(is.element(c("scanID","chrom","left","right"),names(segs)))) stop("some names of RUNS missing")
if(!all(is.element(class(segs),"data.frame"))) stop("segs needs to be a data.frame")
if(!all(is.element(c("scanID","chrom","left","right","sd.fac"),names(segs)))) stop("some names of segs missing")
if(!all(is.element(class(base),"data.frame"))) stop("base needs to be a data.frame")
if(!all(is.element(c("scanID","chrom","chrom.nonanom.median","chrom.nonanom.mad","chrom.nonanom.mean","chrom.nonanom.sd"),names(base)))) stop("some names of base missing")

if(!all(RUNS$scanID==snum & RUNS$chrom==ch)) stop("RUNS not from same sample/chrom")
if(!all(segs$scanID==snum & segs$chrom==ch)) stop("segs not from same sample/chrom")
if(!all(base$scanID==snum & base$chrom==ch)) stop("base not from same sample/chrom")

##### merge segs if appropriate #####
new<- GWASTools:::runsMerge(segs,merge.fac)
FLAG<-new$flag
select<-new$newsegs

#### find overlap of segs with RUNS ##########
if(dim(RUNS)[1]<1) stop("No runs to test")
if(dim(segs)[1]<1)stop("No segment info")
runsegs<-NULL
for(i in 1:dim(RUNS)[1]) {
  for(k in 1:dim(select)[1] ){
   mxL<-max(c(RUNS$left[i],select$left[k]))
   mnR<-min(c(RUNS$right[i],select$right[k]))
   over<-mnR-mxL
   
   if(over<=0)next  #no overlap
     tp<-data.frame(mxL,mnR)
     names(tp)<-c("left","right")
     runsegs<-rbind(runsegs,tp)    
  } 
}

## find stats for the new runs
sd.fac<-NULL;mad.fac<-NULL;seg.median<-NULL
seg.median<-NULL;seg.mean<-NULL;nm<-NULL
 for(k in 1:dim(runsegs)[1]){
   subind<-index>=runsegs$left[k] & index<=runsegs$right[k]
   nm<-c(nm,sum(subind))
   sublrr<-lrr[subind]
   seg.med<-median(sublrr,na.rm=TRUE)
   seg.median<-c(seg.median,seg.med)
   sgm<-mean(sublrr,na.rm=TRUE)
   seg.mean<-c(seg.mean,sgm)
   sdf<-(sgm-base$chrom.nonanom.mean)/base$chrom.nonanom.sd
   mdf<-(seg.med-base$chrom.nonanom.median)/base$chrom.nonanom.mad
   sd.fac<-c(sd.fac,sdf); mad.fac<-c(mad.fac,mdf) 
 }#end of k loop
runsegs$num.mark<-nm
runsegs$seg.median<-seg.median
runsegs$seg.mean<-seg.mean
runsegs$mad.fac<-mad.fac 
runsegs$sd.fac<-sd.fac
runsegs$scanID<-snum
runsegs$chrom<-ch
runsegs$chrom.nonanom.median<-base$chrom.nonanom.median
runsegs$chrom.nonanom.mad<-base$chrom.nonanom.mad
runsegs$chrom.nonanom.mean<-base$chrom.nonanom.mean
runsegs$chrom.nonanom.sd<-base$chrom.nonanom.sd
### select ##

num.segs<-dim(segs)[1]
 
## LOHlocalMad computes and adds the variable "local"
select2<-GWASTools:::LOHlocalMad(runsegs,index,nonanom.index,lrr,length.factor,min.lrr.num)
select2$num.segs<-num.segs

lenrest<-dim(select2)[1]
#will have exited before if there is nothing in select2

## apply filters
jind<-NULL
for(j in 1:lenrest){ 
   
   if(is.na(select2$local[j])) {md<-abs(select2$mad.fac[j]);minv<-md} else {
     md<-mean(c(abs(select2$mad.fac[j]),abs(select2$local[j])))
     minv<-min(c(abs(select2$mad.fac[j]),abs(select2$local[j])))  
   }
   nm<-select2$num.mark[j]
   if(nm>=homodel.min.num & abs(select2$mad.fac[j])>homodel.thresh){ jind<-c(jind,j);next}
   if(!is.na(select2$local[j])& nm>small.num & nm<=medium.num &md>small.thresh & minv>=medium.thresh){ jind<-c(jind,j);next}
   if(is.na(select2$local[j])& nm>small.num & nm<=medium.num & md>small.na.thresh) {jind<-c(jind,j);next}
   if(nm>medium.num & nm<=long.num &md>medium.thresh & minv>=long.thresh) {jind<-c(jind,j);next}
   if(nm>long.num & md>long.thresh & signif(minv,2)>=long.thresh) {jind<-c(jind,j);next}
   
} #end j loop
if(length(jind)!=0) {selt<-select2[jind,];selt<-selt[order(selt$left),]} else selt<-NULL
outt<-list(select2,selt,FLAG)
names(outt)<-c("raw.adjusted","filtered","merge.flag")
return(outt)
  } #end of function
