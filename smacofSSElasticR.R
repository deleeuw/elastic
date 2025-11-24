source("smacofDataUtilities.R")
source("smacofTorgerson.R")
source("smacofPlots.R")

suppressPackageStartupMessages(library(isotone, quietly = TRUE))

smacofSSElasticR <- function(theData,
                      ndim = 2,
                      xinit = NULL,
                      ties = 1,
                      itmax = 1000,
                      eps = 1e-6,
                      verbose = TRUE,
                      ordinal = FALSE,
                      weighted = FALSE) {
  
  ties <- switch(ties, "primary", "secondary", "tertiary")
  ndat <- theData$ndat
  nobj <- theData$nobj
  iind <- theData$iind
  jind <- theData$jind
  dhat <- theData$delta
  wght <- theData$weights
  if (!weighted) {
    wght <- rep(1, ndat)
  }
  if (is.null(xinit)) {
    xinit <- smacofTorgerson(theData, ndim)$conf
  }
  xold <- xinit
  dold <- rep(0, ndat)
  for (k in 1:ndat) {
    dold[k] <- sqrt(sum((xold[iind[k], ] - xold[jind[k], ])^2))
  }
  rold <- dold / dhat
  labd <- sum(wght * rold) / sum(wght * rold^2)
  xold <- labd * xold
  dold <- labd * dold
  rold <- dold / dhat
  sold <- sum(wght * (1 - rold)^2)
  itel <- 1
  repeat {
    vmat <- matrix(0, nobj, nobj)
    bmat <- matrix(0, nobj, nobj)
    for (k in 1:ndat) {
      vmat[iind[k], jind[k]] <- -wght[k] / dhat[k]^2
      bmat[iind[k], jind[k]] <- -wght[k] / (dhat[k] * dold[k])
    }
    vmat <- vmat + t(vmat)
    bmat <- bmat + t(bmat)
    diag(vmat) <- -rowSums(vmat)
    diag(bmat) <- -rowSums(bmat)
    vinv <- solve(vmat + 1 / nobj) - 1 / nobj
    xmid <- vinv %*% bmat %*% xold
    dmid <- rep(0, ndat)
    for (k in 1:ndat) {
      dmid[k] <- sqrt(sum((xmid[iind[k], ] - xmid[jind[k], ])^2))
    }
    rmid <- dmid / dhat
    smid <- sum(wght * (1 - rmid)^2)
    if (ordinal) {
      dhat <- -1 / gpava(theData$delta, -1 / dmid, weights = wght * dmid^2, ties = ties)$x
      rnew <- dmid / dhat
      snew <- sum(wght * (1 - rnew)^2)
    } else {
      snew <- smid
    }
    if (verbose) {
      if (ordinal) {
        cat(
          "itel ",
          formatC(itel, format = "d"),
          "sold ",
          formatC(
            sold,
            digits = 10,
            width = 12,
            format = "f"
          ),
          "smid ",
          formatC(
            smid,
            digits = 10,
            width = 12,
            format = "f"
          ),
          "snew ",
          formatC(
            snew,
            digits = 10,
            width = 12,
            format = "f"
          ),
          "\n"
        )
      } else {
        cat(
          "itel ",
          formatC(itel, format = "d"),
          "sold ",
          formatC(
            sold,
            digits = 10,
            width = 12,
            format = "f"
          ),
          "snew ",
          formatC(
            snew,
            digits = 10,
            width = 12,
            format = "f"
          ),
          "\n"
        )
      }
    }
    if ((itel == itmax) || ((sold - snew) < eps)) {
      break
    }
    itel <- itel + 1
    sold <- snew
    xold <- xmid
    dold <- dmid
  }
  result <- list(
    delta = theData$delta,
    dhat = dhat,
    confdist = dmid,
    conf = xmid,
    weightmat = wght,
    stress = snew,
    ndim = ndim,
    niter = itel,
    nobj = nobj,
    iind = iind,
    jind = jind,
    ordinal = ordinal
  )
  class(result) <- c("smacofSSResult", "smacofSSUOResult")
  return(result)
}