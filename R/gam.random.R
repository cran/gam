#' Specify a Random Effects Fit in a GAM Formula
#' 
#' A symbolic wrapper for a factor term, to specify a random effect term in a
#' formula argument to gam
#' 
#' This "smoother" takes a factor as input and returns a shrunken-mean fit.  If
#' \code{lambda=0}, it simply computes the mean of the response at each level
#' of \code{f}. With \code{lambda>0}, it returns a shrunken mean, where the
#' j'th level is shrunk by \code{nj/(nj+lambda)}, with \code{nj} being the
#' number of observations (or sum of their weights) at level \code{j}. Using
#' such smoother(s) in gam is formally equivalent to fitting a mixed-effect
#' model by generalized least squares.
#' 
#' @aliases random gam.random
#' @param f factor variable, or expression that evaluates to a factor.
#' @param y a response variable passed to \code{gam.random} during backfitting
#' @param w weights
#' @param df the target equivalent degrees of freedom, used as a smoothing
#' parameter. The real smoothing parameter (\code{lambda} below) is found such
#' that \code{df=tr(S)}, where \code{S} is the implicit smoother matrix. Values
#' for \code{df} should be greater than \code{0} and less than the number of
#' levels of \code{f}.  If both \code{df} and \code{lambda} are supplied, the
#' latter takes precedence. Note that \code{df} is not necessarily an integer.
#' @param lambda the non-negative penalty parameter. This is interpreted as a
#' variance ratio in a mixed effects model - namely the ratio of the noise
#' variance to the random-effect variance.
#' @param intercept if \code{intercept=TRUE} (the default) then the estimated
#' level effects are centered to average zero, otherwise they are left alone.
#' @param xeval If this argument is present, then \code{gam.random} produces a
#' prediction at \code{xeval}.
#' @return \code{random} returns the vector \code{f}, endowed with a number of
#' attributes. The vector itself is used in computing the means in backfitting,
#' while the attributes are needed for the backfitting algorithms
#' \code{general.wam}. Note that \code{random} itself does no smoothing; it
#' simply sets things up for \code{gam}.
#' 
#' One important attribute is named \code{call}. For example, \code{random(f,
#' lambda=2)} has a call component \code{gam.random(data[["random(f, lambda =
#' 2)"]], z, w, df = NULL, lambda = 2, intercept = TRUE)}. This is an
#' expression that gets evaluated repeatedly in \code{general.wam} (the
#' backfitting algorithm).
#' 
#' \code{gam.random} returns an object with components \item{residuals}{The
#' residuals from the smooth fit. } \item{nl.df}{the degrees of freedom}
#' \item{var}{the pointwise variance for the fit} \item{lambda}{the value of
#' \code{lambda} used in the fit} When \code{gam.random} is evaluated with an
#' \code{xeval} argument, it returns a vector of predictions.
#' @author Written by Trevor Hastie, following closely the design in the
#' "Generalized Additive Models" chapter (Hastie, 1992) in Chambers and Hastie
#' (1992).
#' @seealso \code{\link{lo}}, \code{\link{s}}, \code{\link{bs}},
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
#' @keywords models regression nonparametric smooth random effects mixed
#' effects
#' @examples
#' 
#' # fit a model with a linear term in Age and a random effect in the factor Level
#' y ~ Age + random(Level, lambda=1)
#'  
#' @export
"gam.random" <-
function(f, y, w, df = sum(non.zero), lambda = 0,intercept=TRUE, xeval)
{
	df.inv <- function(n, df, lambda = sum(n)/df - mean(n), iterations = 10
		)
	{
		if(df > length(n))
			return(0)
		current.df <- sum(n/(n + lambda))
		if(abs((df - current.df)/df) < 0.0001 | iterations == 1)
			lambda
		else {
			lambda <- exp(log(lambda) + (current.df - df)/(sum((
				n * lambda)/(n + lambda)^2)))
			Recall(n, df, lambda, iterations - 1)
		}
	}
        f=attr(f,"values")
	nw <- tapply(w, f, sum)
	non.zero <- !is.na(nw)
	if(is.null(df))
		df <- sum(non.zero)
	if(lambda == 0)
		lambda <- df.inv(nw[non.zero], df)
	df <- sum(nw[non.zero]/(nw[non.zero] + lambda))
	fit <- tapply(w * y, f, sum)/(nw + lambda)
        if(intercept)fit=fit-mean(fit)
	var <- as.vector(w/(nw[f] + lambda))
	residuals <- as.vector(y - fit[f])
	if(missing(xeval))
            list(x = seq(along = nw), y = fit, residuals = residuals, var = var,
                 nl.df = df, lambda = lambda)
        else fit[xeval]
}
