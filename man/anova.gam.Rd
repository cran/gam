\name{anova.gam}
\alias{anova.gam}
\alias{summary.gam}
\title{Analysis of Deviance for a Generalized Additive Model}
\description{Produces an ANODEV table for a set of GAM models, or else a
summary for a single GAM model}
\usage{
\method{anova}{gam}(object, \dots, test)
\method{summary}{gam}(object, dispersion=NULL,\dots)
}
\arguments{
  \item{object}{a fitted gam}
  \item{\dots}{other fitted gams for \code{anova}}
  \item{test}{a character string specifying the test statistic to be used.
          Can be one of '"F"', '"Chisq"' or '"Cp"', with partial
          matching allowed, or 'NULL' for no test.}
	\item{dispersion}{a dispersion parameter to be used in computing
	  standard errors}
	}
\details{
These are  methods for the functions \code{anova} or \code{summary} for objects inheriting from class `gam'.
See \code{\link{anova}} for the general behavior of this function
and for the interpretation of `test'.


When called with a single `gam' object, a special anova table for `gam'
models is returned. This gives a breakdown of the degrees of freedom for
all the terms in the model, separating the projection part and
nonparametric part of each. For example, a term specified by `s()' is
broken down into a single degree of freedom for its linear component,
and the remainder for the nonparametric component. In addition, a type
of score test is performed for each of the nonparametric terms. The
nonparametric component is set to zero, and the linear part is updated,
holding the other nonparametric terms fixed. This is done efficiently
and simulataneously for all terms.  }
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

  Venables, W. N. and Ripley, B. D. (2002)
  \emph{Modern Applied Statistics with S.}
  New York: Springer.
}

\examples{
data(gam.data)
gam.object <- gam(y~s(x,6)+z,data=gam.data)
anova(gam.object)
gam.object2 <- update(gam.object, ~.-z)
anova(gam.object, gam.object2, test="Chisq")
}
\keyword{models}
\keyword{regression}
\keyword{nonparametric}
\keyword{smooth}

