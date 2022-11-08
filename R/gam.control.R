#' Auxilliary for controlling GAM fitting
#' 
#' Auxiliary function as user interface for 'gam' fitting. Typically only used
#' when calling 'gam' or 'gam.fit'.
#' 
#' @param epsilon convergence threshold for local scoring iterations
#' @param bf.epsilon convergence threshold for backfitting iterations
#' @param maxit maximum number of local scoring iterations
#' @param bf.maxit maximum number of backfitting iterations
#' @param trace should iteration details be printed while \code{gam} is fitting
#' the model.
#' @param ... placemark for additional arguments
#' @return a list is returned, consisting of the five parameters, conveniently
#' packaged up to supply the \code{control} argument to \code{gam}. The values
#' for \code{gam.control} can be supplied directly in a call to \code{gam};
#' these are then filtered through \code{gam.control} inside \code{gam}.
#' @references Hastie, T. J. (1992) \emph{Generalized additive models.} Chapter
#' 7 of \emph{Statistical Models in S} eds J. M. Chambers and T. J. Hastie,
#' Wadsworth & Brooks/Cole.
#' @keywords models regression nonparametric smooth
#' @examples
#' 
#' \dontrun{gam(formula, family, control = gam.control(bf.maxit=15))}
#' \dontrun{gam(formula, family, bf.maxit = 15) # these are equivalent}
#' 
#' @export gam.control
"gam.control" <-
function(epsilon = 9.9999999999999995e-08, bf.epsilon = 9.9999999999999995e-08,
	maxit = 30, bf.maxit = 30, trace = FALSE, ...)
{
	if(epsilon <= 0) {
		warning("the value of epsilon supplied is zero or negative; the default value of 1e-7 was used instead"
			)
		epsilon <- 9.9999999999999995e-08
	}
	if(maxit < 1) {
		warning("the value of maxit supplied is too small; the default value of 30 was used instead"
			)
		maxit <- 30
	}
	if(bf.epsilon <= 0) {
		warning("the value of bf.epsilon supplied is zero or negative; the default value of 1e-7 was used instead"
			)
		bf.epsilon <- 9.9999999999999995e-08
	}
	if(bf.maxit < 1) {
		warning("the value of bf.maxit supplied is too small; the default value of 30 was used instead"
			)
		bf.maxit <- 30
	}
	list(epsilon = epsilon, maxit = maxit, bf.epsilon = bf.epsilon, 
		bf.maxit = bf.maxit, trace = as.logical(trace)[1])
}
