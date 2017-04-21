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
