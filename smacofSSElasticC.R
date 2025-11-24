dyn.load("smacofSSElastic.so")

source("smacofAuxiliaries.R")
source("smacofDataUtilities.R")
source("smacofPlots.R")
source("smacofTorgerson.R")

smacofSSElasticC <- function(theData,
                     ndim = 2,
                     xinit = NULL,
                     ties = 1,
                     itmax = 1000,
                     eps = 1e-6,
                     digits = 10,
                     width = 15,
                     verbose = TRUE,
                     weighted = FALSE,
                     ordinal = FALSE) {
  if (is.null(xinit)) {
    xinit <- smacofTorgerson(theData, 2)$conf
  }
  xold <- xinit
  nobj <- theData$nobj
  ndat <- theData$ndat
  iind <- theData$iind
  jind <- theData$jind
  dhat <- theData$delta
  wght <- theData$weights
  if (!weighted) {
    wght <- rep(1, ndat)
  }
  blks <- theData$blocks
  dold <- rep(0, ndat)
  for (k in 1:ndat) {
    dold[k] <- sqrt(sum((xold[iind[k], ] - xold[jind[k], ])^2))
  }
  rold <- dold / dhat
  labd <- sum(wght * rold) / sum(wght * rold^2)
  dold <- labd * dold
  xold <- labd * xold
  rold <- dold / dhat
  sold <- sum(wght * (1 - rold)^2)
  snew <- 0.0
  xold <- as.vector(xold)
  xnew <- xold
  itel <- 1
  h <- .C(
    "smacofSSElasticEngine",
    nobj = as.integer(nobj),
    ndim = as.integer(ndim),
    ndat = as.integer(ndat),
    itel = as.integer(itel),
    ties = as.integer(ties),
    itmax = as.integer(itmax),
    digits = as.integer(digits),
    width = as.integer(width),
    verbose = as.integer(verbose),
    ordinal = as.integer(ordinal),
    weighted = as.integer(weighted),
    sold = as.double(sold),
    snew = as.double(snew),
    eps = as.double(eps),
    iind = as.integer(iind - 1),
    jind = as.integer(jind - 1),
    blks = as.integer(blks),
    wght = as.double(wght),
    edis = as.double(dold),
    dhat = as.double(dhat),
    xold = as.double(xold),
    xnew = as.double(xnew)
  )
  result <- list(
    delta = theData$delta,
    dhat = h$dhat,
    confdist = h$edis,
    conf = matrix(h$xnew, nobj, ndim),
    weightmat = h$wght,
    stress = h$snew,
    ndim = ndim,
    init = xinit,
    niter = h$itel,
    nobj = nobj,
    iind = h$iind,
    jind = h$jind,
    weighted = weighted,
    ordinal = ordinal
  )
  class(result) <- c("smacofSSResult", "smacofSSUOResult")
  return(result)
}
