test_allequal <- function() {
  x <- c(1,2,NA,4)
  y <- x
  checkTrue(allequal(x,y))

  y <- c(1,NA,3,4)
  checkTrue(!allequal(x,y))

  checkTrue(!allequal(x,NULL))
}
