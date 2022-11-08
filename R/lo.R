#' Specify a loess fit in a GAM formula
#' 
#' A symbolic wrapper to indicate a smooth term in a formala argument to gam
#' 
#' A smoother in gam separates out the parametric part of the fit from the
#' non-parametric part. For local regression, the parametric part of the fit is
#' specified by the particular polynomial being fit locally. The workhorse
#' function \code{gam.lo} fits the local polynomial, then strips off this
#' parametric part. All the parametric pieces from all the terms in the
#' additive model are fit simultaneously in one operation for each loop of the
#' backfitting algorithm.
#' 
#' @aliases lo gam.lo
#' @param ...  the unspecified \code{\dots{}} can be a comma-separated list of
#' numeric vectors, numeric matrix, or expressions that evaluate to either of
#' these. If it is a list of vectors, they must all have the same length.
#' @param span the number of observations in a neighborhood. This is the
#' smoothing parameter for a \code{loess} fit. If specified, the full argument
#' name \code{span} must be written.
#' @param degree the degree of local polynomial to be fit; currently restricted
#' to be \code{1} or \code{2}. If specified, the full argument name
#' \code{degree} must be written.
#' @param x for \code{gam.lo}, the appropriate basis of polynomials generated
#' from the arguments to \code{lo}. These are also the variables that receive
#' linear coefficients in the GAM fit.
#' @param y a response variable passed to \code{gam.lo} during backfitting
#' @param w weights
#' @param ncols for \code{gam.lo} the number of columns in \code{x} used as the
#' smoothing inputs to local regression. For example, if \code{degree=2}, then
#' \code{x} has two columns defining a degree-2 polynomial basis. Both are
#' needed for the parameteric part of the fit, but \code{ncol=1} telling the
#' local regression routine that the first column is the actually smoothing
#' variable.
#' @param xeval If this argument is present, then \code{gam.lo} produces a
#' prediction at \code{xeval}.
#' @return \code{lo} returns a numeric matrix.  The simplest case is when there
#' is a single argument to \code{lo} and \code{degree=1}; a one-column matrix
#' is returned, consisting of a normalized version of the vector.  If
#' \code{degree=2} in this case, a two-column matrix is returned, consisting of
#' a degree-2 polynomial basis.  Similarly, if there are two arguments, or the
#' single argument is a two-column matrix, either a two-column matrix is
#' returned if \code{degree=1}, or a five-column matrix consisting of powers
#' and products up to degree \code{2}.  Any dimensional argument is allowed,
#' but typically one or two vectors are used in practice.
#' 
#' The matrix is endowed with a number of attributes; the matrix itself is used
#' in the construction of the model matrix, while the attributes are needed for
#' the backfitting algorithms \code{general.wam} (weighted additive model) or
#' \code{lo.wam} (currently not implemented). Local-linear curve or surface
#' fits reproduce linear responses, while local-quadratic fits reproduce
#' quadratic curves or surfaces. These parts of the \code{loess} fit are
#' computed exactly together with the other parametric linear parts
#' 
#' When two or more smoothing variables are given, the user should make sure
#' they are in a commensurable scale; \code{lo()} does no normalization. This
#' can make a difference, since \code{lo()} uses a spherical (isotropic)
#' neighborhood when establishing the nearest neighbors.
#' 
#' Note that \code{lo} itself does no smoothing; it simply sets things up for
#' \code{gam}; \code{gam.lo} does the actual smoothing. of the model.
#' 
#' One important attribute is named \code{call}. For example, \code{lo(x)} has
#' a call component \code{gam.lo(data[["lo(x)"]], z, w, span = 0.5, degree = 1,
#' ncols = 1)}. This is an expression that gets evaluated repeatedly in
#' \code{general.wam} (the backfitting algorithm).
#' 
#' \code{gam.lo} returns an object with components \item{residuals}{The
#' residuals from the smooth fit. Note that the smoother removes the parametric
#' part of the fit (using a linear fit with the columns in \code{x}), so these
#' residual represent the nonlinear part of the fit.} \item{nl.df}{the
#' nonlinear degrees of freedom} \item{var}{the pointwise variance for the
#' nonlinear fit}
#' 
#' When \code{gam.lo} is evaluated with an \code{xeval} argument, it returns a
#' matrix of predictions.
#' @author Written by Trevor Hastie, following closely the design in the
#' "Generalized Additive Models" chapter (Hastie, 1992) in Chambers and Hastie
#' (1992).
#' @seealso \code{\link{s}}, \code{\link{bs}}, \code{\link{ns}},
#' \code{\link{poly}}, \code{\link{loess}}
#' @references Hastie, T. J. (1992) \emph{Generalized additive models.} Chapter
#' 7 of \emph{Statistical Models in S} eds J. M. Chambers and T. J. Hastie,
#' Wadsworth & Brooks/Cole.
#' 
#' Hastie, T. and Tibshirani, R. (1990) \emph{Generalized Additive Models.}
#' London: Chapman and Hall.
#' @keywords models regression nonparametric smooth
#' @examples
#' 
#' y ~ Age + lo(Start)
#'      # fit Start using a loess smooth with a (default) span of 0.5.
#' y ~ lo(Age) + lo(Start, Number) 
#' y ~ lo(Age, span=0.3) # the argument name span cannot be abbreviated.
#' 
#' @export lo
lo <-
  function (..., span = 0.5, degree = 1) 
{
  vars <- list(...)
  locall <- sys.call()
  chcall <- deparse(locall)
  nvars <- length(vars)
  if (degree > 2) 
    stop("degrees 1 or 2 are implemented")
  if (nvars == 1) {
    xvar <- as.matrix(vars[[1]])
    xnames <- deparse(locall[[2]])
    if(is.null(dimnames(xvar)[[2]])){
      nc=ncol(xvar)
      dxnames=xnames
      if(nc>1)dxnames=paste(xnames,1:nc,sep=".")
      dimnames(xvar)=list(NULL,dxnames)
    }
  }
  else {
    nobs <- length(vars[[1]])
    xvar <- matrix(0, nobs, nvars)
    xnames <- character(nvars)
    for (i in seq(nvars)) {
      tt <- vars[[i]]
      if (!is.null(dd <- dim(tt)) && dd[2] > 1) 
        stop("either call lo with a matrix argument, or else a comma separated list x1, x2")
      exptt <- locall[[i + 1]]
      xnames[i] <- deparse(exptt)
      xvar[, i] <- as.numeric(tt)
    }
  dimnames(xvar) <- list(NULL, xnames)
  }
 
  polyx <- polylo(xvar, degree = degree)
  pd <- attr(polyx, "degree")
  opd <- order(pd)
  if (length(pd) > 1) {
    polyx <- polyx[, opd]
    p <- sum(pd == 1)
  }
  else p <- 1
  nobs <- dim(polyx)[1]
  nas <- is.na(polyx[, 1:p])
  if (any(nas)) {
    if (p > 1) 
      nas <- nas %*% array(1, c(p, 1))
    attr(polyx, "NAs") <- seq(nobs)[nas > 0]
  }
  real.call <- substitute(gam.lo(data[[chcall]], z, w, span = span, 
                                 degree = degree, ncols = p), list(span = span, degree = degree, 
                                                    chcall = chcall, p = p))
  atts <- c(attributes(polyx), list(span = span, degree = degree, 
                                    ncols = p, call = real.call))
  attributes(polyx) <- atts
  class(polyx) <- c("smooth", "matrix")
  polyx
}
