#####
# Identify relatives from ibd coefficient estimates
#####

ibdAssignRelatedness <- function(
  k0,				# vector of k0 estimates
  k1,				# parallel vector of k1 estimates
  alpha=.05,          # significance level - finds 100(1-alpha)% predicition intervals/ellipse   	
  m = 0.04,		# width of rectangle along diagaonal line (distance from diagonal)
  po.w = 0.1,		# width of parent-offspring rectangle
  po.h = 0.1,		# height of parent-offspring rectangle
  dup.w = 0.1,	# width of duplicate rectangle
  dup.h = 0.1,		# height of duplicate rectangle
  un.w = 0.25,	# width of unrelated rectangle
  un.h = 0.25		# height of unrelated rectangle
){
  # returns a vector of assignments to "PO", "FS", "HS", "FC","U" (unrelated) and "Q" for everything else, for each pair of (k0,k1)
  
  # if(length(k0)!=length(k1)) stop("k0 and k1 must be parallel vectors of the same length")	
  
  # to identify PO within the rectangle
  po.sel <- k0 < po.w & k1 > 1-po.h
  # to identify duplicates within rectangle
  dup.sel <- k0 < dup.w & k1 < dup.h
  
  # identify full sibs, half sibs, first cousins
  
  # to identify full sibs within the ellipse - 
  
  FS<-GWASTools::relationsMeanVar$FullSibs
  mean.vec<-FS$mean  #vector
  sig.inv<-FS$invCov
  n0<-length(k0)
  X<-matrix(c(k0,k1),2,n0,byrow=TRUE)
  Mn<-matrix(rep(mean.vec,n0),2,n0)
  diff<-X-Mn
  tmp<-sig.inv%*%diff
  prob<-rep(NA,n0)
  for(i in 1:n0){
    prob[i]<-1-pchisq(diff[,i]%*%tmp[,i],2)
  }
  fs.sel<-prob>= alpha #inside ellipse
  
  sdm<-abs(qnorm(alpha/2))  # sd multiplier for determining half-sib, first cousin rectangles
  # to identify half-sibs within rectangle parallel to diagonal
  # ends of the rectangle on the diagonal
  HS<-GWASTools::relationsMeanVar$HalfSibs
  d<-sdm*sqrt(HS$var) # +/- d from k1-mean gives 100(1-alpha)% prediction interval for k1
  hsm<-HS$mean[2]
  y1<-hsm-d;x1<-1-y1
  y0<-hsm+d; x0<-1-y0
  # find the nearest point on the diagonal line
  x <- (1 + k0 - k1)/2; y <- (1 + k1 - k0)/2
  # is it within the diagonal segment
  chk1 <- x0 < x & x < x1 & y1 < y & y < y0  
  # is the point within perpendicular distance m of the diagonal
  chk2 <- (k0 - x)^2 + (k1 - y)^2 < m^2
  hs.sel <- chk1 & chk2
  
  # to identify first cousins within the rectangle parallel to diagonal
  # ends of the rectangle on the diagonal
  C<-GWASTools::relationsMeanVar$FirstCousins
  d<-sdm*sqrt(C$var)
  fcm<-C$mean[2]
  y1<-fcm-d;x1<-1-y1
  y0<-fcm+d; x0<-1-y0
  # find the nearest point on the diagonal line
  x <- (1 + k0 - k1)/2; y <- (1 + k1 - k0)/2
  # is it within the diagonal segment
  chk1 <- x0 < x & x < x1 & y1 < y & y < y0  
  # is the point within perpendicular distance m of the diagonal
  chk2 <- (k0 - x)^2 + (k1 - y)^2 < m^2
  fc.sel <- chk1 & chk2
  # check for overlap
  rels <- po.sel & dup.sel & fs.sel & hs.sel & fc.sel
  if(any(rels)) stop("one or more pairs assigned to more than one relationship")
  # unrelated
  rels <- po.sel | dup.sel | fs.sel | hs.sel | fc.sel
  un.sel <- !rels & k0 > 1 - un.w & k1 < un.h
  # combine logical vectors to get vector of assignments
  asnmt <- rep("Q", length(k0))
  asnmt[po.sel] <- "PO"
  asnmt[dup.sel] <- "Dup"
  asnmt[fs.sel] <- "FS"
  asnmt[hs.sel] <- "Deg2"
  asnmt[fc.sel] <- "Deg3"
  asnmt[un.sel] <- "U"
  return(asnmt)
}



ibdAssignRelatednessKing <- function(
  ibs0,  		# vector of ibs0 estimates
  kc,  			# vector of kinship coefficient estimates
  cut.kc.dup=1/(2^(3/2)),  # kinship coefficient threshold for duplicates
  cut.kc.fs=1/(2^(5/2)), # kc threshold for deg 1 relatives
  cut.kc.deg2=1/(2^(7/2)), # kc threshold for deg2 relatives
  cut.kc.deg3=1/(2^(9/2)), # kc threshold for deg3 relatives
  cut.ibs0.err=0.003 # should be 0 for PO, but sometimes is greater due to genotyping error.
){
  
  ## returns a vector of assignments to "PO", "FS", "HS", "FC","U" (unrelated) and "Q" for everything else, for each pair of (k0,k1)
  
  ## default thresholds for assigning relationships use kinship coefficients in table 1 of Manichaikul (2010) - KING paper
  
  # if(length(ibs0)!=length(kc)) stop("ibs0 and kc must be parallel vectors of the same length")  
  
  
  dup.sel <- kc > cut.kc.dup
  po.sel <- kc <= cut.kc.dup & kc > cut.kc.fs & ibs0 <= cut.ibs0.err
  fs.sel <- kc <= cut.kc.dup & kc > cut.kc.fs & ibs0 > cut.ibs0.err
  d2.sel <- kc <= cut.kc.fs & kc > cut.kc.deg2
  d3.sel <- kc <= cut.kc.deg2 & kc > cut.kc.deg3
  un.sel <- kc <= cut.kc.deg3
  
  # check for overlap - should be none
  rels <- po.sel & dup.sel & fs.sel & d2.sel & d3.sel
  if (any(rels)) stop("one or more pairs assigned to more than one relationship")
  
  rels <- po.sel | dup.sel | fs.sel | d2.sel | d3.sel
  un.sel <- !rels 
  
  asnmt <- rep("Q", length(kc))
  asnmt[dup.sel] <- "Dup"
  asnmt[po.sel] <- "PO"
  asnmt[fs.sel] <- "FS"
  asnmt[d2.sel] <- "Deg2"
  asnmt[d3.sel] <- "Deg3"
  asnmt[un.sel] <- "U"
  
  return(asnmt)
  
}
