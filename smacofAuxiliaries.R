
makeMPInverseV <- function(theData) {
  nobj <- theData$nobj
  ndat <- theData$ndat
  wght <- matrix(0, nobj, nobj)
  for (k in 1:ndat) {
    i <- theData$iind[k]
    j <- theData$jind[k]
    wght[i, j] <- wght[j, i] <- theData$weights[k]
  }
  vmat <- -wght
  diag(vmat) <- -rowSums(vmat)
  vinv <- solve(vmat + (1 / nobj)) - (1 / nobj)
  return(as.vector(as.dist(vinv)))
}

smacofRandomConfiguration <- function(theData, ndim = 2) {
  nobj <- theData$nobj
  x <- matrix(rnorm(nobj * ndim), nobj, ndim)
  return(columnCenter(x))
}

columnCenter <- function(x) {
  apply(x, 2, function(x) x - mean(x))
}

matrixPrint <- function(x,
                        digits = 6,
                        width = 8,
                        format = "f",
                        flag = "+") {
  print(noquote(
    formatC(
      x,
      digits = digits,
      width = width,
      format = format,
      flag = flag
    )
  ))
}
