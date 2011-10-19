ibdPlot <- function(k0, k1, alpha=.05, relation=NULL, color=NULL, rel.lwd=2, ...) {

  stopifnot(length(k0) == length(k1))
  if (!is.null(relation)) stopifnot(length(relation) == length(k0))
  if (!is.null(color)) stopifnot(length(color) == length(k0))
  
  n <- length(k0)
  if (is.null(color) & is.null(relation)) {
    color <- rep("black", n)
  } else if (is.null(color)) {
    stopifnot(length(relation) == length(k0))
    stopifnot(all(relation %in% c("Dup", "PO", "FS", "HS", "FC", "U", "Q")))
    color <- rep("black", n)
    color[relation == "Dup"] <- "magenta"
    color[relation == "PO"] <- "cyan"
    color[relation == "FS"] <- "red"
    color[relation == "HS"] <- "blue"
    color[relation == "FC"] <- "lightgreen"
    color[relation == "Q"] <- "darkgreen"
  }
  
  # main plot - no data yet
  plot(k0, k1, xlim=c(0,1), ylim=c(0,1), xlab="k0", ylab="k1", type="n", ...)

  # light grey background boxes
  abline(h=seq(0,1,0.25), lty=2, col="gray")
  abline(v=seq(0,1,0.25), lty=2, col="gray")

  # k0 + k1 = 1
  abline(1, -1, lty=2)

  # add points
  points(k0, k1, col=color, ...)

  # delineate the expected values
  env = new.env(parent=emptyenv())
  data(relationsMeanVar, envir=env)
  sdm<-abs(qnorm(alpha/2))
  
  # draw full-sib ellipse
  FS = env[["relationsMeanVar"]]$FullSibs
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
  points(pts, type="l", lwd=rel.lwd, col="orange")

  HS<-env[["relationsMeanVar"]]$HalfSibs
  d<-sdm*sqrt(HS$var) # +/- d from k1-mean gives 100(1-alpha)% prediction interval for k1
  hsm<-HS$mean[2]
  s2<-sqrt(2)
  y1<-hsm-d;x1<-1-y1
  y0<-hsm+d; x0<-1-y0
  segments(x0,y0,x1,y1,lwd=rel.lwd,col="orange")
 
  C<-env[["relationsMeanVar"]]$FirstCousins
  d<-sdm*sqrt(C$var)
  fcm<-C$mean[2]
  s2<-sqrt(2)
  y1<-fcm-d;x1<-1-y1
  y0<-fcm+d; x0<-1-y0
  segments(x0,y0,x1,y1,lwd=rel.lwd,col="orange")
 
  # legend
  if (!is.null(relation)) {
    legend("topright", unique(relation), col=unique(color),
           pch=rep(1, length(unique(relation))))
  }

}
