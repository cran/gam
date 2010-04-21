\name{predict.gam}
\alias{predict.gam}
\title{Predict method for GAM fits}
\description{Obtains predictions and optionally estimates standard errors of
     those predictions from a fitted generalized additive model object.}
\usage{
predict.gam(object, newdata, type, dispersion, se.fit = FALSE,na.action, terms,\dots)
}
\arguments{
\item{object}{
a fitted \code{gam} object, or one of its inheritants, such as a \code{glm} or \code{lm} object.
}
\item{newdata}{
a data frame containing the values at which predictions are required.
This argument can be missing, in which case predictions are made at the same values used to compute the object. 
Only those predictors, referred to in the right side of the formula in object need be present by name in \code{newdata}.
}
\item{type}{
type of predictions, with choices \code{"link"} (the default), \code{"response"}, or \code{"terms"}.
The default produces predictions on the scale of the additive predictors, and with \code{newdata} missing, \code{predict} is simply an extractor function for this component of a \code{gam} object. If \code{"response"} is selected, the predictions are on the scale of the response, and are monotone transformations of the additive predictors, using the inverse link function. If \code{type="terms"} is selected, a matrix of predictions is produced, one column for each term in the model.
}
\item{se.fit}{
if \code{TRUE}, pointwise standard errors are computed along with the predictions.
}
\item{dispersion}{the dispersion of the GLM fit to be assumed in computing
          the standard errors.  If omitted, that returned by 'summary'
          applied to the object is used}
\item{terms}{
if \code{type="terms"}, the \code{terms=} argument can be used to specify which terms should be included; the default is \code{labels(object)}. 
}
\item{na.action}{ function determining what should be done with missing values
  in 'newdata'.  The default is to predict 'NA'.}
\item{\dots}{Placemark for additional arguments to predict}
}
\value{
a vector or matrix of predictions, or a list consisting of the predictions and their standard errors if \code{se.fit = TRUE}. 
If \code{type="terms"}, a matrix of fitted terms is produced, with one column for each term in the model (or subset of these if the \code{terms=} argument is used). There is no column for the intercept, if present in the model, and each of the terms is centered so that their average over the original data is zero. 
The matrix of fitted terms has a \code{"constant"} attribute which, when added to the sum of these centered terms, gives the additive predictor.
See the documentation of \code{predict} for more details on the components returned.

When \code{newdata} are supplied, \code{predict.gam} simply invokes
inheritance and gets \code{predict.glm} to produce the parametric part
of the predictions. For each nonparametric term, \code{predict.gam} reconstructs the
partial residuals and weights from the final iteration of the local
scoring algorithm. The appropriate smoother is called for each term,
with the appropriate \code{xeval} argument (see \code{\link{s}} or
\code{\link{lo}}), and the prediction for that term is produced.

The standard errors are based on an approximation given in Hastie
(1992). Currently \code{predict.gam} does not produce standard errors for
predictions at \code{newdata}.

Warning: naive use of the generic
\code{predict} can produce incorrect predictions when the \code{newdata}
argument is used, if the formula in \code{object} involves
transformations such as \code{sqrt(Age - min(Age))}. 
}
\seealso{\code{\link{predict.glm}},
\code{\link{fitted}}, \code{\link{expand.grid}}
}
\author{
  Written by Trevor Hastie, following closely the design in the
  "Generalized Additive Models" chapter (Hastie, 1992) in Chambers and
  Hastie (1992). This version of \code{predict.gam} is adapted from the S
  version to match the corresponding predict methods for \code{glm} and
  \code{lm} objects in R. The \code{safe.predict.gam} function in S is
  no longer required, primarily because a safe prediction method is in
  place for functions like \code{ns}, \code{bs}, and \code{poly}.
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
gam.object <- gam(y ~ s(x,6) + z, data=gam.data)
predict(gam.object) # extract the additive predictors
data(gam.newdata)
predict(gam.object, gam.newdata, type="terms") 
}

\keyword{models}
\keyword{regression}
\keyword{nonparametric}
\keyword{smooth}
