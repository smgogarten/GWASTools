# qqPlot <- function(pval, truncate = FALSE, ylim = NULL, ...)
# {
# 	# pvalue is a vector of p-values
# 	# truncate is T/F whether to truncate the y-axis (observed) to the same limits as the x-axis (expected)
# 
# 	pval <- -log10(sort(pval)) # sort() removes NAs
# 	n <- length(pval)
#         a <- 1:n
#         b <- a/n
#         upper <- qbeta(0.025, a, rev(a))
#         lower <- qbeta(0.975, a, rev(a))
# 	x <- -log10(b)
# 
# 
#         char <- rep(1,n)
#         if(!truncate){
#           ylm <- ylim
#           ylb <- expression(paste(-log[10], "(observed P)"))
#         }else{
#           maxx <- max(x)+2
#           ylm <- c(0,maxx)
#           ylb <- expression(paste(-log[10], "(observed P) - truncated"))
#           nx <- length(which(pval > maxx))
#           if(nx > 0){
#             pval[1:nx] <- maxx
#             char[1:nx] <- 2
#           }
#         }
#         plot(x, pval, type = "n", ylim = ylm, ylab = ylb,
#              xlab = expression(paste(-log[10], "(expected P)")), ...)
#         polygon(-log10(c(b,rev(b))), -log10(c(upper, rev(lower))), density=NA, col="gray")
#         points(x, pval, pch = char, ...)
# 	abline(0,1,col="red")
# }


qqPlot <- function(pval, truncate = FALSE, ylim=NULL, thinThreshold=NULL, ...) 
{
  # pvalue is a vector of p-values
  # truncate is T/F whether to truncate the y-axis (observed) to the same limits as the x-axis (expected)
  # thin is whether to thin insignificant p-values
  
  pval <- -log10(sort(pval)) # sort() removes NAs
  n <- length(pval)
  a <- 1:n
  b <- a/n
  x <- -log10(b)
  
  if (!is.null(thinThreshold)){
    breaks <- seq(from=0, to=thinThreshold, length.out=11)
    pval.cut <- cut(pval, breaks=breaks, right=FALSE)
    #quant <- quantile(pval[pval < thinThreshold], probs=1:10/10)
    #pval.cut <- cut(pval, breaks=c(0, quant), right=FALSE)
    ind.list <- list()
    for (level in levels(pval.cut)){
      sel <- which(pval.cut == level)
      ind.list[[level]] <- sel[sample.int(length(sel), min(1000, length(sel)))]
    }
    ind.list[["max"]] <- which(pval >= thinThreshold)
    ind <- unlist(ind.list, use.names=FALSE)
    ind <- sort(ind) # sorting necessary for polygon, below
  } else {
    ind <- 1:n
  }
  
  upper <- qbeta(0.025, a, rev(a))
  lower <- qbeta(0.975, a, rev(a))
  
  char <- rep(1,n)
  if(!is.logical(truncate) | truncate){
    if (is.logical(truncate)){
      maxx <- max(x)+2  
    } else {
      maxx <- min(truncate, max(pval))
    }
    
    ylm <- c(0,maxx)
    ylb <- expression(paste(-log[10], "(observed P) - truncated"))
    nx <- length(which(pval > maxx))
    if(nx > 0){
      pval[1:nx] <- maxx
      char[1:nx] <- 2
    }
    
  } else {
    ylm <- ylim
    ylb <- expression(paste(-log[10], "(observed P)"))
  }
  plot(x[ind], pval[ind], type = "n", ylim = ylm, ylab = ylb,
       xlab = expression(paste(-log[10], "(expected P)")), ...)
  # upper and lower have already been subset
  polygon(-log10(c(b[ind], rev(b[ind]))), -log10(c(upper[ind], rev(lower[ind]))), density=NA, col="gray")
  points(x[ind], pval[ind], pch = char, ...)  
  abline(0,1,col="red")  
}
