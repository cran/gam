\name{gam.control}
\alias{gam.control}
\title{Auxilliary for controlling GAM fitting}
\description{Auxiliary function as user interface for 'gam' fitting. Typically
     only used when calling 'gam' or 'gam.fit'.}
\usage{
gam.control(epsilon=1e-07, bf.epsilon = 1e-07, maxit=30, bf.maxit = 30, trace=FALSE,\ldots)
}
\arguments{
\item{epsilon}{
convergence threshold for local scoring iterations
}
\item{bf.epsilon}{
convergence threshold for backfitting iterations
}
\item{maxit}{
maximum number of local scoring iterations
}
\item{bf.maxit}{
maximum number of backfitting iterations
}
\item{trace}{
should iteration details be printed while \code{gam} is fitting the model.
}
\item{\ldots}{}
}
\value{
a list is returned, consisting of the five parameters, conveniently packaged up to supply the \code{control} argument to \code{gam}. The values for \code{gam.control} can be supplied directly in a call to \code{gam}; these are then filtered through \code{gam.control} inside \code{gam}.
}
\references{
  Hastie, T. J. (1992)
  \emph{Generalized additive models.}
  Chapter 7 of \emph{Statistical Models in S}
  eds J. M. Chambers and T. J. Hastie, Wadsworth \& Brooks/Cole.
  }
\examples{
\dontrun{gam(formula, family, control = gam.control(bf.maxit=15))}
\dontrun{gam(formula, family, bf.maxit = 15) # these are equivalent}
}
\keyword{models}
\keyword{regression}
\keyword{nonparametric}
\keyword{smooth}
