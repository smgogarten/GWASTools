test_genotypeToCharacter <- function() {
  g <- c(0,1,2,NA)
  A <- rep("G",4)
  B <- rep("C",4)
  checkIdentical(c("C/C","G/C","G/G",NA),
                 genotypeToCharacter(g, A, B, sort=FALSE))
  checkIdentical(c("C/C","C/G","G/G",NA),
                 genotypeToCharacter(g, A, B, sort=TRUE))
  checkIdentical(c("B/B","A/B","A/A",NA),
                 genotypeToCharacter(g, NULL, NULL, sort=FALSE))
  checkIdentical(matrix(c("C/C","G/C","G/G",NA), nrow=2),
                 genotypeToCharacter(matrix(g, nrow=2), A[1:2], B[1:2], sort=FALSE))
  checkIdentical(c("C/C","G/C","G/G",NA),
                 genotypeToCharacter(g, "G", "C", sort=FALSE))
}

test_genotypeToCharacter_NA <- function() {
  A <- c("A", NA, "A", NA)
  B <- c("C", "C", NA, NA)
  g <- rep(0,4)
  checkIdentical(c("C/C","C/C",NA,NA),
                 genotypeToCharacter(g, A, B))
  g <- rep(1,4)
  checkIdentical(c("A/C",NA,NA,NA),
                 genotypeToCharacter(g, A, B))
  g <- rep(2,4)
  checkIdentical(c("A/A",NA,"A/A",NA),
                 genotypeToCharacter(g, A, B))
}

test_pasteSorted <- function() {
    a <- c("a","b")
    b <- c("b","a")
    checkIdentical(c("a/b", "a/b"), pasteSorted(a,b))
    checkIdentical(c("ab", "ab"), pasteSorted(a,b, sep=""))
}
