dyn.load("smacofSSWOEngine.so")

source("smacofUtils.R")
source("smacofSSData.R")

smacofSSWO <- function(theData,
                       ndim = 2,
                       xold = torgerson(theData, ndim),
                       ties = 1,
                       itmax = 1000,
                       eps = 1e-10,
                       digits = 10,
                       width = 15,
                       verbose = TRUE) {
  nobj <- nrow(xold)
  ndat <- nrow(theData)
  itel <- 1
  iind <- theData[, 1]
  jind <- theData[, 2]
  dhat <- theData[, 3]
  blks <- theData[, 4]
  wght <- theData[, 5]
  edis <- rep(0, ndat)
  for (k in 1:ndat) {
    i <- iind[k]
    j <- jind[k]
    edis[k] <- sqrt(sum((xold[i, ] - xold[j, ])^2))
  }
  sdd <- sum(wght * edis^2)
  sde <- sum(wght * dhat * edis)
  lbd <- sde / sdd
  edis <- lbd * edis
  xold <- lbd * xold
  sold <- sum(wght * (dhat - edis)^2) / wsum
  snew <- 0.0
  xold <- as.vector(xold)
  xnew <- rep(0, nobj * ndim)
  h <- .C(
    "smacofSSWOEngine",
    nobj = as.integer(nobj),
    ndim = as.integer(ndim),
    ndat = as.integer(ndat),
    itel = as.integer(itel),
    ties = as.integer(ties),
    itmax = as.integer(itmax),
    digits = as.integer(digits),
    width = as.integer(width),
    verbose = as.integer(verbose),
    wsum = as.double(wsum),
    sold = as.double(sold),
    snew = as.double(snew),
    eps = as.double(eps),
    iind = as.integer(iind),
    jind = as.integer(jind),
    blks = as.integer(blks),
    edis = as.double(edis),
    dhat = as.double(dhat),
    wght = as.double(wght),
    vinv = as.double(vinv),
    xold = as.double(xold),
    xnew = as.double(xnew)
  )
  return(
    list(
      data = theData,
      conf = h$xnew,
      loss = h$snew,
      edis = h$edis,
      dhat = h$dhat,
      itel = h$itel
    )
  )
}
