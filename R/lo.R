"lo" <-
  function(..., span = 0.5, degree = 1)
{
  lodummy <- function(span = 0.5, degree = 1)
    list(span = span, degree = degree)
  vars <- list(...)
  locall <- sys.call()
  chcall <- deparse(locall)
  mcall <- match.call(expand = FALSE)
  mcall$... <- NULL
  nvars <- length(vars)
  if(nvars > 1) {
    scalars <- sapply(vars, length) == 1
    ## a bit of freedom in giving the span and degree
    if(any(scalars)) {
      nvars <- nvars - sum(scalars)
      mcall <- c(mcall, as.call(vars[scalars]))
      vars <- vars[!scalars]
    }
  }
  mcall[[1]] <- as.name("lodummy")
  m <- eval(mcall)
  degree <- m$degree
  span <- m$span
  if(degree > 2.)
    stop("degrees 1 or 2 are implemented")
  if(nvars == 1) {
    xvar <- as.matrix(vars[[1]])
    xnames <- deparse(locall[[2.]])
  }
  else {
    nobs <- length(vars[[1]])
    xvar <- matrix(0., nobs, nvars)
    xnames <- character(nvars)
    for(i in seq(nvars)) {
      tt <- vars[[i]]
      if(!is.null(dd <- dim(tt)) && dd[2.] > 1)
        stop("either call lo with a matrix argument, or else a comma separated list x1, x2"
             )
      exptt <- locall[[i + 1]]
      xnames[i] <- deparse(exptt)
      xvar[, i] <- as.numeric(tt)
    }
    dimnames(xvar) <- list(NULL, xnames)
  }
  ## for the moment we use polybasis from library(mda)
  polyx <- polylo(xvar, degree = degree)
  pd <- attr(polyx, "degree")
  opd <- order(pd)
  if(length(pd) > 1) {
    polyx <- polyx[, opd]
    p <- sum(pd == 1)
  }
  else p <- 1
  nobs <- dim(polyx)[1]
  nas <- is.na(polyx[, 1:p])
  if(any(nas)) {
    if(p > 1)
      nas <- nas %*% array(1, c(p, 1))
    attr(polyx, "NAs") <- seq(nobs)[nas > 0.]
  }
##  if(span * nobs < 1)
##    stop(paste("span is too small; the minimum is 1/n =", format(
##                                                                 round(1/nobs, 4.))))
  real.call <- substitute(gam.lo(data[[chcall]], z, w, span = span, 
                                 degree = degree, ncols = p), list(span = span, degree = degree,
                                                    chcall = chcall, p = p))
  atts <- c(attributes(polyx), list(span = span, degree = degree, ncols = 
                                   p, call = real.call))
  attributes(polyx) <- atts
  class(polyx) <- c("smooth", "matrix")
  polyx
}
