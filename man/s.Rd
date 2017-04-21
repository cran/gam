\name{s}
\alias{s}
\alias{gam.s}
\title{Specify a Smoothing Spline Fit in a GAM Formula}
\description{A symbolic wrapper to indicate a smooth term in a formala
 argument to gam}
\usage{
s(x, df=4, spar=1)
gam.s(x, y, w, df, spar, xeval)
}
\arguments{
\item{x}{
the univariate predictor, or expression, that evaluates to a numeric vector.
}

\item{df}{ the target equivalent degrees of freedom, used as a smoothing
parameter. The real smoothing parameter (\code{spar} below) is found
such that \code{df=tr(S)-1}, where \code{S} is the implicit smoother
matrix. Values for \code{df} should be greater than \code{1}, with
\code{df=1} implying a linear fit. If both \code{df} and \code{spar} are
supplied, the former takes precedence. Note that \code{df} is not
necessarily an integer.}
\item{spar}{
can be used as smoothing parameter, with values typically in
\code{(0,1]}. See \code{\link{smooth.spline}} for more details.}

\item{y}{a response variable passed to \code{gam.s} during backfitting}
\item{w}{weights}
\item{xeval}{If this argument is present, then \code{gam.s} produces a
  prediction at \code{xeval}.}
}
\value{

\code{s} returns the vector \code{x}, endowed with a number of
attributes. The vector itself is used in the construction of the model
matrix, while the attributes are needed for the backfitting algorithms
\code{general.wam} (weighted additive model) or \code{s.wam}.
Since smoothing splines reproduces linear fits, the linear
part will be efficiently computed with the other parametric linear parts
of the model.

Note that \code{s} itself does no smoothing; it simply sets things up
for \code{gam}.

One important attribute is named \code{call}. For example, \code{s(x)}
has a call component
\code{gam.s(data[["s(x)"]], z, w, spar = 1, df = 4)}.
This is an expression that gets evaluated repeatedly in \code{general.wam}
(the backfitting algorithm).

\code{gam.s} returns an object with components
\item{residuals}{The residuals from the smooth fit. Note that the
  smoother removes the parametric part of the fit (using a linear fit
  in \code{x}), so these residual represent the
  nonlinear part of the fit.}
\item{nl.df}{the nonlinear degrees of freedom}
\item{var}{the pointwise variance for the nonlinear fit}

When \code{gam.s} is evaluated with an \code{xeval} argument, it returns a
vector of predictions.
}
\seealso{
\code{\link{lo}}, \code{\link{smooth.spline}}, \code{\link{bs}}, \code{\link{ns}}, \code{\link{poly}}
}
\author{
  Written by Trevor Hastie, following closely the design in the
  "Generalized Additive Models" chapter (Hastie, 1992) in Chambers and
  Hastie (1992).
 }
\references{
  Hastie, T. J. (1992)
  \emph{Generalized additive models.}
  Chapter 7 of \emph{Statistical Models in S}
  eds J. M. Chambers and T. J. Hastie, Wadsworth \& Brooks/Cole.

  Hastie, T. and Tibshirani, R. (1990)
  \emph{Generalized Additive Models.}
  London: Chapman and Hall.

  Cantoni, E. and hastie, T. (2002)
  Degrees-of-freedom tests for smoothing splines,
  \emph{Biometrika} 89(2), 251-263
}
\examples{
     # fit Start using a smoothing spline with 4 df.
     y ~ Age + s(Start, 4)
     # fit log(Start) using a smoothing spline with 5 df.
     y ~ Age + s(log(Start), df=5)
 }
\keyword{models}
\keyword{regression}
\keyword{nonparametric}
\keyword{smooth}

