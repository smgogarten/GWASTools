qqPlot <- function(pval, truncate = FALSE, ylim = NULL, ...)
{
	# pvalue is a vector of p-values
	# truncate is T/F whether to truncate the y-axis (observed) to the same limits as the x-axis (expected)

	pval <- -log10(sort(pval)) # sort() removes NAs
	n <- length(pval)
        a <- 1:n
        b <- a/n
        upper <- qbeta(0.025, a, rev(a))
        lower <- qbeta(0.975, a, rev(a))
	x <- -log10(b)


        char <- rep(1,n)
        if(!truncate){
          ylm <- ylim
          ylb <- expression(paste(-log[10], "(observed P)"))
        }else{
          maxx <- max(x)+2
          ylm <- c(0,maxx)
          ylb <- expression(paste(-log[10], "(observed P) - truncated"))
          nx <- length(which(pval > maxx))
          if(nx > 0){
            pval[1:nx] <- maxx
            char[1:nx] <- 2
          }
        }
        plot(x, pval, type = "n", ylim = ylm, ylab = ylb,
             xlab = expression(paste(-log[10], "(expected P)")), ...)
        polygon(-log10(c(b,rev(b))), -log10(c(upper, rev(lower))), density=NA, col="gray")
        points(x, pval, pch = char, ...)
	abline(0,1,col="red")
}

