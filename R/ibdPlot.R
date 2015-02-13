ibdPlot <- function(k0, k1, alpha=0.05, relation=NULL, color=NULL, rel.lwd=2,
                    rel.draw=c("FS", "Deg2", "Deg3"), ...) {

  stopifnot(length(k0) == length(k1))
  if (!is.null(relation)) stopifnot(length(relation) == length(k0))
  if (!is.null(color)) stopifnot(length(color) == length(k0))
  stopifnot(is.null(rel.draw) | all(rel.draw %in% c("FS", "Deg2", "Deg3")))
  
  n <- length(k0)
  ordered.cols <- FALSE
  if (is.null(color) & is.null(relation)) {
    color <- rep("black", n)
  } else if (is.null(color)) {
    stopifnot(length(relation) == length(k0))
    relation[relation %in% c("HS", "HSr", "HSFC", "Av", "GpGc", "DFC")] <- "Deg2"
    relation[relation %in% c("FC", "HAv", "OAv", "OC")] <- "Deg3"
    relation[!(relation %in% c("Dup", "PO", "FS", "Deg2", "Deg3", "Q"))] <- "U"
    
    rels <- c("Dup", "PO", "FS", "Deg2", "Deg3", "Q", "U")
    cols <- c("magenta", "cyan", "red", "blue", "lightgreen", "darkgreen", "black")
    names(cols) <- rels
    color <- cols[relation]
    ordered.cols <- TRUE
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
  sdm<-abs(qnorm(alpha/2))
  
  # draw full-sib ellipse
  if ("FS" %in% rel.draw) {
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
    points(pts, type="l", lwd=rel.lwd, col="orange")
  }

  if ("Deg2" %in% rel.draw) {
    HS<-GWASTools::relationsMeanVar$HalfSibs
    d<-sdm*sqrt(HS$var) # +/- d from k1-mean gives 100(1-alpha)% prediction interval for k1
    hsm<-HS$mean[2]
    s2<-sqrt(2)
    y1<-hsm-d;x1<-1-y1
    y0<-hsm+d; x0<-1-y0
    segments(x0,y0,x1,y1,lwd=rel.lwd,col="orange")
  }

  if ("Deg3" %in% rel.draw) {
    C<-GWASTools::relationsMeanVar$FirstCousins
    d<-sdm*sqrt(C$var)
    fcm<-C$mean[2]
    s2<-sqrt(2)
    y1<-fcm-d;x1<-1-y1
    y0<-fcm+d; x0<-1-y0
    segments(x0,y0,x1,y1,lwd=rel.lwd,col="orange")
  }
 
  # legend
  if (!is.null(relation)) {
    if (ordered.cols) {
      rel <- rels[rels %in% unique(relation)]
      col <- cols[rel]
    } else {
      rel <- unique(relation)
      col <- unique(color)
      ord <- order(rel)
      rel <- rel[ord]
      col <- col[ord]
    }
    legend("topright", legend=rel, col=col, pch=1)
  }
}
