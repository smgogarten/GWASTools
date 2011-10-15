########################################################################
#
# To randomly select snps among which no pair is closer than min.dist apart
#
#  chromosome -- chromosome as integer
#	anything but 1-26 is missing
#  position -- genetic position of snps
#  min.dist -- the distance (basepair)
#  init.sel -- a logical vector indicating the initial snps
#  max.n.chromosomes -- the max number of snps for each chromosome,
#	"-1" means no number limit.
#  verbose -- show information when running
#
#  return a logical vector indicating the selection.
#
########################################################################
apartSnpSelection <- function(
	chromosome, 
	position, 
	min.dist = 100000,
	init.sel = NULL,
	max.n.chromosomes = -1,
	verbose = TRUE)
{
	# check
	if (!(length(chromosome) == length(position)))
		stop("Lengths of chomosome.code and position do not match")
	
	if (is.null(init.sel))
		init.sel <- rep(TRUE, length(chromosome))
	
	if (!is.logical(init.sel))
		stop("init.sel must be of type logical")

	# init
	rv <- rep(FALSE, length(chromosome))

	for (chr in 1:26)
	{
		b <- (chromosome == chr & !is.na(chromosome))
		sel <- which(init.sel[b])
		pos <- position[b]
		flag <- rep(FALSE, length(pos))
		iter.num <- 0
		show.flag <- TRUE

		repeat
		{
			if (length(sel) == 0) break
			if (iter.num == max.n.chromosomes) break
			iter.num <- iter.num + 1

			# new position
			if (length(sel) > 1)
				p.i <- sample(sel, 1)
			else
				p.i <- sel[1]
			flag[p.i] <- TRUE
			sel <- sel[abs(pos[sel] - pos[p.i]) >= min.dist]

			if (verbose & (iter.num %% 1000 == 0))
			{
				show.flag <- FALSE
				message(paste(date(), "Chromosome", chr, "num.snp", iter.num))
			}
		}

		if (verbose & show.flag)
			message(paste(date(), "Chromosome", chr))

		rv[b] <- flag
	}

	if (verbose)
		message(paste(date(), "Total # of SNPs selected:", sum(rv)))

	return(rv)
}

