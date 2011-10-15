## covariance matrix, inverse, eigenvalues - needed for full sib prediction  ellipse
## also include variance for half-sibs, cousins - for prediction interval
## from Hill, W.G. and B.S. Weir, Variation in actual relationship 
#     as a consequence of Mendelian sampling and linkage, Genet. Res., Camb. (2011), 93, 47--64
## uses map length data from Kong, X., Murphy, K., Rag, T., He, C., White, P.S. and
#  Matise, T.C., A combined physical-linkage map of the human genome. (2004) American
#  Journal of Human Genetics 75, 1143--1148 

## sub-function phi
phi<-function(len,n) {
  #len = map length
  #n: integer indicating level of function
coeff<-1/(2*len^2) * .25^n
out<-0
for(r in 1:n){
  s<-choose(n,r)*((2*r*len-1+exp(-2*r*len))/r^2)
  out<-out+s
}
out<-coeff*out
return(out)
}

##################################
## compute variance, covariance across each chromosome
## independent so can "add" - to get autosomal variance, use weighted average
#1/(sum of lengths)^1[length(i)^2*variance(i)]

## map lengths in cM from Kong, et al for chromosomes 1 - 22
chrom.len<-c(286.5, 263.3, 225.1, 212.2, 208.2, 192.2, 189.0, 173.3, 168.7, 173.5, 163.8, 174.2,
  128.9 ,123.8, 130.2, 134.2, 137.5, 124.2, 112.2, 102.5,  68.5,  86.1)

L<-sum(chrom.len)/100
## see formulas in Box 1 of Hill/Weir

### FULL SIB
vR<-0; vk1<-0; vk0<-0;cv10<-0 
for(i in 1:22){
  len<-chrom.len[i]/100  
  vrR<-2*phi(len,2)-phi(len,1)
  vrk0<-16*phi(len,4)-16*phi(len,3)+8*phi(len,2)-2*phi(len,1)
  vrk1<-4*vrk0-4*vrR
  cvr10<- -2*vrk0+2*vrR
  vR<-vR+len^2*vrR
  vk1<-vk1+len^2*vrk1
  vk0<-vk0+len^2*vrk0
  cv10<-cv10+len^2*cvr10
}
vR<-1/L^2*vR 
vk1<-1/L^2*vk1
vk0<-1/L^2*vk0
cv10<-1/L^2*cv10

mean.vec<-c(.25,.5)
sigma<-matrix(c(vk0,cv10,cv10,vk1),2,2)
sigma.inv<-solve(sigma)
eg<-eigen(sigma.inv)

FullSibs<-list(mean.vec, sigma,sigma.inv,eg$values,eg$vectors)
names(FullSibs)<-c("mean","cov","invCov","eigvals","eigvectors")
###################

## Half Sib & First Cousins ###
mean.vec.hs<-c(.5,.5)
mean.vec.fc<-c(.75,.25)
vhs<-0
vfc<-0
for(i in 1:22){
  len<-chrom.len[i]/100  
  vrhs<-4*phi(len,2)-2*phi(len,1)
  vrfc<-8*phi(len,4)-4*phi(len,3)+1.5*phi(len,2)-.5*phi(len,1)
  vhs<-vhs+len^2*vrhs
  vfc<-vfc+len^2*vrfc
}
vhs<-vhs/L^2
vfc<-vfc/L^2

HalfSibs<-list(mean.vec.hs,vhs)
names(HalfSibs)<-c("mean","var")

FirstCousins<-list(mean.vec.fc,vfc)
names(FirstCousins)<-c("mean","var")

MeanVarRelations<-list(FullSibs,HalfSibs,FirstCousins)
names(MeanVarRelations)<-c("FullSibs","HalfSibs","FirstCousins")

##saved as a data object relationsMeanVar.RData in svn trunk/data/







