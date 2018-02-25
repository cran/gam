\name{step.Gam}
\alias{step.Gam}
\title{Stepwise model builder for GAM}
\description{Builds a GAM model in a step-wise fashion. For each "term"
  there is an ordered list of alternatives, and the function traverses
  these in a greedy fashion. Note: this is NOT a method for \code{step},
which used to be a generic, so must be invoked with the full name.}
\usage{
step.Gam(object, scope, scale, direction, trace, keep, steps, parallel, \ldots)
}
\arguments{
\item{object}{
An object of class \code{Gam} or any of it's inheritants.
}
\item{scope}{
defines the range of models examined in the step-wise search. It is a list of formulas, with each formula corresponding to a term in the model. Each of these formulas specifies a "regimen" of candidate forms in which the particular term may enter the model. For example, a term formula might be 
\code{~1+ Income + log(Income) + s(Income)}.
This means that \code{Income} could either appear not at all,  linearly, linearly in its logarithm, or as a smooth function estimated nonparametrically. A \code{1} in the formula allows the additional option of leaving the term out of the model entirely. 
Every term in the model is described by such a term formula, and the
final model is built up by selecting a component from each
formula.

As an alternative more convenient for big models, each list can have
instead of a formula a
character vector corresponding to the candidates for that term. Thus we
could have \code{c("1","x","s(x,df=5")} rather than \code{~1+x+s(x,df=5)}.

The supplied model \code{object} is used as the starting model, and hence there is the requirement that one term from each of the term formulas be present in \code{formula(object)}. This also implies that any terms in \code{formula(object)} \emph{not} contained in any of the term formulas will be forced
to be present in every model considered.
The function \code{gam.scope} is helpful for generating the scope
argument for a large model.
}
\item{scale}{
an optional argument used in the definition of the AIC statistic used to evaluate models for selection. By default, the scaled Chi-squared statistic for the initial model is used, but if forward selection is to be performed, this is not necessarily a sound choice. 
}
\item{direction}{
The mode of step-wise search, can be one of \code{"both"}, \code{"backward"}, or \code{"forward"}, with a default of \code{"both"}. If \code{scope} is missing, the default for \code{direction} is "both". 
}
\item{trace}{
If \code{TRUE} (the default), information is printed during the running
of \code{step.Gam()}. This is an encouraging choice in general, since
\code{step.Gam()} can take some time to compute either for large models
or when called with an an extensive \code{scope=} argument. A simple one
line model summary is printed for each model selected. This argument can
also be given as the binary \code{0} or \code{1}. A value \code{trace=2}
gives a more verbose trace.
}
\item{keep}{
A filter function whose input is a fitted \code{Gam} object, and
anything else passed via \dots, and whose output is arbitrary. Typically \code{keep()} will select a subset of the components of the object and return them. The default is not to keep anything.
}
\item{steps}{
The maximum number of steps to be considered. The default is 1000 (essentially as
many as required). It is typically used to stop the process early.
}
\item{parallel}{If \code{TRUE}, use parallel \code{foreach} to fit each
  trial run.
      Must register parallel before hand, such as \code{doMC} or others.
    See the example below.}

\item{\dots}{Additional arguments to be passed on to \code{keep}}
}
\value{
The step-wise-selected model is returned, with up to two additional components.
There is an \code{"anova"} component corresponding to the steps taken in the search, as well as a \code{"keep"} component if the \code{keep=} argument was supplied in the call. 


We describe the most general setup, when \code{direction = "both"}.
At any stage there is a current model comprising a single term from each of the term formulas supplied in the \code{scope=} argument.
A series of models is fitted, each corrresponding to a formula obtained by moving each of the terms one step up or down in its regimen, relative to the formula of the current model.
If the current value for any term is at either of the extreme ends of its regimen, only one rather than two steps can be considered.
So if there are \code{p} term formulas, at most \code{2*p - 1}  models are considered.
A record is kept of all the models ever visited (hence the \code{-1} above), to avoid repetition.
Once each of these models has been fit, the "best" model
in terms of the AIC statistic is selected and defines the step.
The entire process is repeated until either the maximum number of steps has been used, or until the AIC criterion can not be decreased by any of the eligible steps.
}
\seealso{
\code{\link{gam.scope}},\code{\link{step}},\code{\link{glm}}, \code{\link{gam}}, \code{\link{drop1}}, \code{\link{add1}}, \code{\link{anova.Gam}}
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
Gam.object <- gam(y~x+z, data=gam.data)
step.object <-step.Gam(Gam.object, scope=list("x"=~1+x+s(x,4)+s(x,6)+s(x,12),"z"=~1+z+s(z,4)))
\dontrun{
# Parallel
require(doMC)
registerDoMC(cores=2)
step.Gam(Gam.object, scope=list("x"=~1+x+s(x,4)+s(x,6)+s(x,12),"z"=~1+z+s(z,4)),parallel=TRUE)
}

}
\keyword{models}
\keyword{regression}
\keyword{nonparametric}
\keyword{smooth}

