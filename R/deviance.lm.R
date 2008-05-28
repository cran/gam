"deviance.lm" <-
function(object, ...)
if(is.null(w <- object$weights)) sum(object$residuals^2.) else sum(w * object$
		residuals^2.)
