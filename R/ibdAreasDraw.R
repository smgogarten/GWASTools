#####
# Draw relationship areas for ibd analysis
#####


ibdAreasDraw <- function(
      alpha=.05,         # significance level - finds 100(1-alpha)% prediction intervals/ellipse
	m = 0.04,		# width of rectangle along diagaonal line
	po.w = 0.1,		# width of parent-offspring rectangle
	po.h = 0.1,		# height of parent-offspring rectangle
	dup.w = 0.1,	# width of duplicate rectangle
	dup.h = 0.1,		# height of duplicate rectangle
	un.w = 0.25,	# width of unrelated rectangle
	un.h = 0.25,		# height of unrelated rectangle
        rel.lwd = 2,
	xcol = c("cyan","red","blue","lightgreen","magenta","black")	# colors for parent-offspring, full-sib, half-sib, first cousin, dup & unrelated areas
)
{
      sdm<-abs(qnorm(alpha/2)) #sd multiplier for rectangles for 100(1-alpha)% prediction interval

      # draw full-sib ellipse
      FS<-GWASTools::relationsMeanVar$FullSibs
      mean.vec<-FS$mean  #vector
      sig.inv<-FS$invCov            
      eg.vals<-FS$eigvals  #relates to length of ellipse axes
      eg.vec<-FS$eigvectors  #gives direction of ellipse axes
      c2<-qchisq(1-alpha,2)
      h<-sqrt(c2/eg.vals)   # gives lengths of ellipse axes
      ang <- seq(0, 2 * pi, by=.01) # angles for parameterization of ellipse
      
      mean.vec2<-matrix(mean.vec,2,1) # make column vector
      pts<-NULL
      for(a in ang){
        p1<-h[1]*cos(a) ;p2<-h[2]*sin(a)
        tp<-t(p1*eg.vec[,1]+p2*eg.vec[,2]+mean.vec2)
        pts<-rbind(pts,tp)
      }
      points(pts, type="l", lwd=rel.lwd, col=xcol[2])

      # draw rectangle for half-sibs
      HS<-GWASTools::relationsMeanVar$HalfSibs
      d<-sdm*sqrt(HS$var) # +/- d from k1-mean gives 100(1-alpha)% prediction interval for k1
      hsm<-HS$mean[2]
      s2<-sqrt(2)
      y1<-hsm-d;x1<-1-y1
      y0<-hsm+d; x0<-1-y0
      a<-x0-m/s2;b<-y0-m/s2
      c<-x1-m/s2;d<-y1-m/s2
      segments(a,b,x0,y0,col=xcol[3], lwd=rel.lwd)
      segments(c,d,x1,y1,col=xcol[3], lwd=rel.lwd)
      segments(a,b,c,d,col=xcol[3], lwd=rel.lwd)
      segments(x0,y0,x1,y1,col=xcol[3], lwd=rel.lwd)

      # draw rectangle for first cousins
      C<-GWASTools::relationsMeanVar$FirstCousins
      d<-sdm*sqrt(C$var)
      fcm<-C$mean[2]
      s2<-sqrt(2)
      y1<-fcm-d;x1<-1-y1
      y0<-fcm+d; x0<-1-y0
      a<-x0-m/s2;b<-y0-m/s2
      c<-x1-m/s2;d<-y1-m/s2
      segments(a,b,x0,y0,col=xcol[4], lwd=rel.lwd)
      segments(c,d,x1,y1,col=xcol[4], lwd=rel.lwd)
      segments(a,b,c,d,col=xcol[4], lwd=rel.lwd)
      segments(x0,y0,x1,y1,col=xcol[4], lwd=rel.lwd)

      # draw rectangle for parent-offspring
      rect(0,1-po.h,po.w,1, border=xcol[1], lwd=rel.lwd)
      # draw rectangle for duplicates
      rect(0,0,dup.w,dup.h, border=xcol[5], lwd=rel.lwd)
      # draw rectangle for unrelated 
      rect(1-un.w, 0, 1, un.h, border=xcol[6], lwd=rel.lwd)
}



