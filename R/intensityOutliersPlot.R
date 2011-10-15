#####
# Function to plot mean intensity for chromosome i vs mean of intensities for autosomes (excluding i) and highlight outliers
#####

intensityOutliersPlot <- 
function(
	mean.intensities, 		# sample x chromosome matrix of mean intensities
	sex,		
	outliers,	
	sep = FALSE,
	label,
	...)
{
#	mean.intensities, 	# sample x chromosome matrix of mean intensities
#	sex,			# vector with values of "M" or "F" corresponding to samples in the rows of mean.intensities
#	outliers,		# list of outliers, each member corresponds to a chromosome (member "X" is itself a list of female and male outliers)
#	sep=F,			# plot outliers within a chromosome separately (T) or together (F)
#	label			# list of plot labels (to be positioned below X axis) corresponding to list of outliers
	
	
	dat<- mean.intensities

	# v2 adds label
	# v3 corrects error in reversing axis labels on x and y axes

	# check args
	if(!all(is.element(colnames(dat), c(1:22, "X", "XY", "Y", "M", "U")))) stop("column names should be 1-22, X, XY, Y, M and (optionally) U")
	if(!all(!is.na(sex) & is.element(sex, c("M","F")))) stop("all elements of sex should be either 'M' or 'F'")
	if(!length(sex)==dim(dat)[1]) stop("length of sex must match the number of rows of 'mean.intensities'")
	if(!is.list(outliers)) stop("outliers should be a list")
	if(!all(is.element(names(outliers), colnames(dat)))) stop("names of 'outliers' should be in column names of ' mean.intensities'")

	# prepare to plot

	chr <- names(outliers)
	data <- dat[, is.element(colnames(dat),1:22)]
	
	# plot for each member of 'outliers'
	for(i in chr)
	{
		ylab <- paste("chromosome",i); xlab="other autosomes"
		y <- mean.intensities[, i]

		dato <- data[, !is.element(colnames(data),i), drop = FALSE]
		x <- apply(dato, 1, mean)
		out <- outliers[[i]]
		# non-sex chromosomes
		if(!is.element(i, c("X","Y")))
		{
			ni <- length(out)
			coef <- lm(y~x)$coefficients
			if (!sep) 
			{ 
				plot(x, y, main=paste("chr",i),xlab=xlab, ylab=ylab, ...)
				abline(coef, col="red")
				sel <- is.element(rownames(dat),out)  # select outliers
				points(x[sel], y[sel], col="red", pch=18)
			} else
			{
				ni <- length(out)
				for(j in 1:ni){
					plot(x, y, main=paste("chr",i, "sample", out[j]), xlab=xlab, ylab=ylab, sub=label[[i]][j], ...)
					abline(coef, col="red")
					sel <- is.element(rownames(dat),out[j])  # select outliers
					points(x[sel], y[sel], col="red", pch=18)
				}
			}
		}
		# X chromosome
		if(is.element(i, "X"))
		{
			for(k in c("F","M"))
			{
				out2 <- out[[k]]
				ni <- length(out2)
				x2 <- x[sex==k]; y2 <- y[sex==k]
				if(!sep) 
				{ 
					plot(x2, y2, main=paste("chr",i,"- sex",k), xlab=xlab, ylab=ylab, ...)
					abline(lm(y2~x2)$coefficients, col="red")
					sel <- is.element(rownames(dat),out2)  # select outliers
					points(x[sel], y[sel], col="red", pch=18)
				} else 
				{
					ni <- length(out2)
					for(j in 1:ni)
					{
						plot(x2, y2, main=paste("chr",i,"- sex",k, "- sample", out2[j]), xlab=xlab, ylab=ylab, sub=label[[i]][[k]][j], ...)
						abline(lm(y2~x2)$coefficients, col="red")
						sel <- is.element(rownames(dat),out2[j])  # select outliers
						points(x[sel], y[sel], col="red", pch=18)
					}
				}
			} # end for k
		}  # end if
		# Y chromosome
		if(is.element(i,"Y"))
		{
			ni <- length(out)
			coef <- lm(y[sex=="M"]~x[sex=="M"])$coefficients
			if(!sep) 
			{ 
				plot(x[sex=="M"], y[sex=="M"], main=paste("chr",i,"- sex M"),xlab=xlab, ylab=ylab, ...)
				abline(coef, col="red")
				sel <- is.element(rownames(dat),out)  # select outliers
				points(x[sel], y[sel], col="red", pch=18, ...)
			} else 
			{
				ni <- length(out)
				for(j in 1:ni)
				{
					plot(x[sex=="M"], y[sex=="M"], main=paste("chr",i,"- sex M", "- sample", out[j]), xlab=xlab, ylab=ylab, sub=label[[i]][j], ...)
					abline(coef, col="red")
					sel <- is.element(rownames(dat),out[j])  # select outliers
					points(x[sel], y[sel], col="red", pch=18)
				}
			}
		} # end if
	} # end for i

}


		










