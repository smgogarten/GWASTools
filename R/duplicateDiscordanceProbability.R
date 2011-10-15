# returns table of discordance probabilities for combined error rates of
# 1e-5, 1e-4, 1e-3, 1e-2, given npair observed pairs
duplicateDiscordanceProbability <- function(npair,
                                              error.rate = c(1e-5, 1e-4, 1e-3, 1e-2),
                                              max.disc = 7) {
  
  # since there are three genotypes, one possible call is correct and the other
  # two are erroneous, so theoretically we have two error rates, a and b
  # we assume a=b (previous work shows it doesn't matter)
  # error.rate argument to this function is the total error rate a+b

  #prob of alternate genotype 1
  a <- error.rate/2
  #prob of alternate genotype 2
  b <- a
  nerrs <- length(error.rate)
  P <- rep(NA, nerrs) # prob of discordance for one pair of samples at one SNP

  # get prob discordance for each pair of errors
  for(i in 1:nerrs){
    probs <- c(1-a[i]-b[i], a[i], b[i])
    perror <- probs %*% t(probs)
    pcon <- sum(diag(perror))
    #same as pcon = sum(probs^2)
    P[i] <- 1-pcon
  }
  #P is probability of discordance for each error rate

  d1 <- paste("dis>", 0:max.disc, sep="")
  d2 <- paste("error=", error.rate, sep="")
  n1 <- length(d1)
  n2 <- length(d2)
  dis <- matrix(nrow=n1, ncol=n2, dimnames=list(d1,d2))

  for(i in 1:nerrs){
    zero <- dbinom(0, npair, P[i]) # prob 0 disc
    dis[1,i] <- 1-zero # prob >0 disc
    if (max.disc > 0) {
      for (j in 1:max.disc) {
        obs <- dbinom(j, npair, P[i]) # prob j disc
        dis[j+1,i] <- dis[j,i] - obs # prob >i disc
      }
    }
  }

  return(dis)
}
