\name{plot.Gam}
\alias{plot.Gam}
\alias{preplot.Gam}
\alias{plot.preplot.Gam}
\title{Plot Components of a GAM Object}
\description{A plot method for GAM objects, which can be used on GLM and
  LM
  objects as well. It focuses on terms (main-effects), and produces a
  suitable plot for terms of different types}
\usage{
\method{plot}{Gam}(x, residuals, rugplot, se, scale, ask = FALSE,terms,\dots)
\method{preplot}{Gam}(object, newdata, terms,\dots)
}
\arguments{
\item{x}{
a \code{Gam} object, or a \code{preplot.Gam} object. The first thing \code{plot.Gam()} does is check if \code{x} has a component called \code{preplot}; if not, it computes one using \code{preplot.Gam()}. Either way, it is this \code{preplot.Gam} object that is required for plotting a \code{Gam} object.
}
\item{object}{same as \code{x}}
\item{residuals}{
if \code{TRUE}, partial deviance residuals are plotted along with the fitted terms---default is \code{FALSE}. If \code{residuals} is a vector with the same length as each fitted term in \code{x}, then these are taken to be the overall residuals to be used for constructing the partial residuals.
}
\item{rugplot}{
if \code{TRUE} (the default), a univariate histogram or \code{rugplot} is displayed along the base of each plot, showing the occurrence of each `x'; ties are broken by jittering.
}
\item{se}{
if \code{TRUE}, upper and lower pointwise twice-standard-error curves are included for each plot. The default is \code{FALSE}.
}
\item{scale}{
a lower limit for the number of units covered by the limits on the `y' for each plot. The default is \code{scale=0}, in which case each plot uses the range of the functions being plotted to create their \code{ylim}. By setting \code{scale} to be the maximum value of \code{diff(ylim)} for all the plots, then all subsequent plots will produced in the same vertical units. This is essential for comparing the importance of fitted terms in additive models.
}
\item{ask}{
if \code{TRUE}, \code{plot.Gam()} operates in interactive mode. 
}
\item{newdata}{if supplied to \code{preplot.Gam}, the preplot object is
  based on them rather than the original.}
\item{terms}{subsets of the terms can be selected}
\item{\dots}{Additonal plotting arguments, not all of which will work
  (like xlim)}
}
\value{
a plot is produced for each of the terms in the object \code{x}. The function currently knows how to plot all main-effect functions of one or two predictors. So in particular, interactions are not plotted. An appropriate `x-y' is produced to display each of the terms, adorned with residuals, standard-error curves, and a rugplot, depending on the choice of options. The form of the plot is different, depending on whether the `x'-value for each plot is numeric, a factor, or a matrix.


When \code{ask=TRUE}, rather than produce each plot sequentially, \code{plot.Gam()} displays a menu listing all the terms that can be plotted, as well as switches for all the options. 


A \code{preplot.Gam} object is a list of precomputed terms. Each such
term (also a \code{preplot.Gam} object) is a list with components
\code{x}, \code{y} and others---the basic ingredients needed for each
term plot. These are in turn handed to the specialized plotting function
\code{gplot()}, which has methods for different classes of the leading
\code{x} argument. In particular, a different plot is produced if
\code{x} is numeric, a category or factor, a matrix, or a
list. Experienced users can extend this range by creating more
\code{gplot()} methods for other classes.  Graphical parameters (see
\code{\link{par}}) may also be supplied as arguments to this function.
This function is a method for the generic function \code{plot()} for
class \code{"Gam"}.

It can be invoked by calling \code{plot(x)} for an
object \code{x} of the appropriate class, or directly by
calling \code{plot.Gam(x)} regardless of the
class of the object.
}
\seealso{
\code{\link{preplot}}, \code{\link{predict.Gam}}
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
}
\examples{
data(gam.data)
Gam.object <- gam(y ~ s(x,6) + z,data=gam.data)
plot(Gam.object,se=TRUE)
data(gam.newdata)
preplot(Gam.object,newdata=gam.newdata)
}
\keyword{models}
\keyword{regression}
\keyword{nonparametric}
\keyword{smooth}

