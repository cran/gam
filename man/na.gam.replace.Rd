\name{na.gam.replace}
\alias{na.gam.replace}
\title{Missing Data Filter for GAMs}
\description{A method for dealing with missing values, friendly to GAM models.}
\usage{
na.gam.replace(frame)
}
\arguments{
\item{frame}{
a model or data frame
}}
\value{
a model or data frame is returned, with the missing observations (NAs) replaced.
The following rules are used. A factor with missing data is replaced by a new factor with one more level, labelled \code{"NA"}, which records the missing data. 
Ordered factors are treated similarly, except the result is an unordered factor.
A missing numeric vector has its missing entires replaced by the mean of the non-missing entries. Similarly, a matrix with missing entries has each missing entry replace by the mean of its column.
If \code{frame} is a model frame, the response variable can be identified, as can the weights (if present). Any rows for which the response or weight is missing are removed entirely from the model frame.


The word \code{"gam"} in the name is relevant, because \code{gam()} makes special use of this filter. All columns of a model frame that were created by a call to \code{lo()} or \code{s()} have an attribute names \code{"NAs"} if NAs are present in their columns. 
Despite the replacement by means, these attributes remain on the object,
and \code{gam()} takes appropriate action when smoothing against these
columns. See section 7.3.2 in Hastie (1992)
for more details.
}
\author{Trevor Hastie}

\references{
  Hastie, T. J. (1992)
  \emph{Generalized additive models.}
  Chapter 7 of \emph{Statistical Models in S}
  eds J. M. Chambers and T. J. Hastie, Wadsworth \& Brooks/Cole.
}

\examples{
data(airquality)
gam(Ozone^(1/3) ~ lo(Solar.R) + lo(Wind, Temp), data=airquality, na=na.gam.replace)
}
\seealso{
\code{\link{na.fail}}, \code{\link{na.omit}}, \code{\link{gam}}
}
\keyword{models}
\keyword{regression}
\keyword{nonparametric}
\keywords{smooth}

