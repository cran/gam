#' Specify a Smoothing Spline Fit in a GAM Formula
#' 
#' A symbolic wrapper to indicate a smooth term in a formala argument to gam
#' 
#' 
#' @aliases s gam.s
#' @param x the univariate predictor, or expression, that evaluates to a
#' numeric vector.
#' @param df the target equivalent degrees of freedom, used as a smoothing
#' parameter. The real smoothing parameter (\code{spar} below) is found such
#' that \code{df=tr(S)-1}, where \code{S} is the implicit smoother matrix.
#' Values for \code{df} should be greater than \code{1}, with \code{df=1}
#' implying a linear fit. If both \code{df} and \code{spar} are supplied, the
#' former takes precedence. Note that \code{df} is not necessarily an integer.
#' @param spar can be used as smoothing parameter, with values typically in
#' \code{(0,1]}. See \code{\link{smooth.spline}} for more details.
#' @param y a response variable passed to \code{gam.s} during backfitting
#' @param w weights
#' @param xeval If this argument is present, then \code{gam.s} produces a
#' prediction at \code{xeval}.
#' @return
#' 
#' \code{s} returns the vector \code{x}, endowed with a number of attributes.
#' The vector itself is used in the construction of the model matrix, while the
#' attributes are needed for the backfitting algorithms \code{general.wam}
#' (weighted additive model) or \code{s.wam}. Since smoothing splines
#' reproduces linear fits, the linear part will be efficiently computed with
#' the other parametric linear parts of the model.
#' 
#' Note that \code{s} itself does no smoothing; it simply sets things up for
#' \code{gam}.
#' 
#' One important attribute is named \code{call}. For example, \code{s(x)} has a
#' call component \code{gam.s(data[["s(x)"]], z, w, spar = 1, df = 4)}. This is
#' an expression that gets evaluated repeatedly in \code{general.wam} (the
#' backfitting algorithm).
#' 
#' \code{gam.s} returns an object with components \item{residuals}{The
#' residuals from the smooth fit. Note that the smoother removes the parametric
#' part of the fit (using a linear fit in \code{x}), so these residual
#' represent the nonlinear part of the fit.} \item{nl.df}{the nonlinear degrees
#' of freedom} \item{var}{the pointwise variance for the nonlinear fit}
#' 
#' When \code{gam.s} is evaluated with an \code{xeval} argument, it returns a
#' vector of predictions.
#' @author Written by Trevor Hastie, following closely the design in the
#' "Generalized Additive Models" chapter (Hastie, 1992) in Chambers and Hastie
#' (1992).
#' @seealso \code{\link{lo}}, \code{\link{smooth.spline}}, \code{\link{bs}},
#' \code{\link{ns}}, \code{\link{poly}}
#' @references Hastie, T. J. (1992) \emph{Generalized additive models.} Chapter
#' 7 of \emph{Statistical Models in S} eds J. M. Chambers and T. J. Hastie,
#' Wadsworth & Brooks/Cole.
#' 
#' Hastie, T. and Tibshirani, R. (1990) \emph{Generalized Additive Models.}
#' London: Chapman and Hall.
#' 
#' Cantoni, E. and hastie, T. (2002) Degrees-of-freedom tests for smoothing
#' splines, \emph{Biometrika} 89(2), 251-263
#' @keywords models regression nonparametric smooth
#' @examples
#' 
#'      # fit Start using a smoothing spline with 4 df.
#'      y ~ Age + s(Start, 4)
#'      # fit log(Start) using a smoothing spline with 5 df.
#'      y ~ Age + s(log(Start), df=5)
#'  
#' @export
"gam.s" <-
  function(x, y, w = rep(1, length(x)), df = 4, spar = 1, xeval)
{
  storage.mode(x) <- storage.mode(y) <- storage.mode(w) <- storage.mode(
                                                                        spar) <- storage.mode(df) <- "double"
  n <- as.integer(length(x))
  x <- signif(x, 6)
  mat <- gam.match(x)
  omat <- mat$o
  nef <- mat$nef
  ##
  ## in rgam.r, splsm calls both splsm1 and splsm2.
  ## splsm2 needs (10+2*4)*(nef+2)+5*nef+n+15 doubles for work.
  ## splsm1 needs 3*nef+2*n+10.
  work.len <- max(3 * nef + 2 * n + 10, (10 + 2 * 4) * (nef + 2) + 5 *
                  nef + n + 15)
  fit <- .Fortran("splsm",
                  x,
                  y,
                  w,
                  n,
                  omat,
                  nef,
                  spar = spar,
                  df = df,
                  s = double(n),
                  s0 = double(1),
                  var = double(nef),
                  FALSE,
                  work = double(work.len),
                  PACKAGE="gam")
  if(missing(xeval))
    list(residuals = y - fit$s, nl.df = fit$df - 1, var = fit$
         var[omat])
  else {
    skn <- .Fortran("sknotl",
                    fit$work[seq(nef)],
                    nef,
                    knot = double(nef + 6),
                    k = integer(1),
                    PACKAGE="gam")
    smallest <- x[omat == 1][1]
    largest <- x[omat == nef][1]
    k <- skn$k
    gam.sp(xeval, skn$knot[seq(k)], k - 4, fit$work[seq(3 * nef +
                                                        n + 10, length = k - 4)], smallest, largest - smallest)
  }
}
