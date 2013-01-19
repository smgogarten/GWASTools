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
