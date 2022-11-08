#' Generate a scope for step.Gam
#' 
#' Given a data.frame as an argument, generate a scope list for use in
#' step.Gam, each element of which gives the candidates for that term.
#' 
#' This function creates a similar scope formula for each variable in the
#' frame. A column named "x" by default will generate a scope term
#' \code{~1+x+s(x)}. With \code{arg=c("df=4","df=6")} we get
#' \code{~1+x+s(x,df=4)+s(x,df=6)}. With form=FALSE, we would get the character
#' vector \code{c("1","x","s(x,df=4)","s(x,df=6")}.
#' 
#' @param frame a data.frame to be used in \code{step.Gam}. Apart from the
#' response column, all other columns will be used.
#' @param response The column in \code{frame} used as the response. Default is
#' 1.
#' @param smoother which smoother to use for the nonlinear terms; i.e. "s" or
#' "lo", or any other supplied smoother. Default is "s".
#' @param arg a character (vector), which is the argument to \code{smoother}.
#' For example, \code{arg="df=6"} would result in the expression
#' \code{s(x,df=6)} for a column named "x". This can be a vector, for example
#' \code{arg=c("df=4","df=6")}, which would result two smooth terms.
#' @param form if \code{TRUE}, each term is a formula, else a character vector.
#' @return a scope list is returned, with either a formula or a character
#' vector for each term, which describes the candidates for that term in the
#' Gam.
#' @author Written by Trevor Hastie, following closely the design in the
#' "Generalized Additive Models" chapter (Hastie, 1992) in Chambers and Hastie
#' (1992).  This version of \code{gam.scope} is adapted from the S version.
#' @seealso \code{\link{step.Gam}}
#' @references Hastie, T. J. (1991) \emph{Generalized additive models.} Chapter
#' 7 of \emph{Statistical Models in S} eds J. M. Chambers and T. J. Hastie,
#' Wadsworth & Brooks/Cole.
#' @keywords models regression nonparametric smooth
#' @examples
#' 
#' data(gam.data)
#' gdata=gam.data[,1:3]
#' gam.scope(gdata,2)
#' gam.scope(gdata,2,arg="df=5")
#' gam.scope(gdata,2,arg="df=5",form=FALSE)
#' gam.scope(gdata,2,arg=c("df=4","df=6"))
#' 
#' @export gam.scope
"gam.scope" <-
function(frame, response = 1, smoother = "s", arg = NULL, form = TRUE)
{
	vnames <- names(frame)
	vnames <- vnames[ - response]
	step.list <- as.list(vnames)
	names(step.list) <- vnames
	for(vname in vnames) {
		junk <- c("1", vname)
		if(is.vector(frame[[vname]]))
			junk <- c(junk, paste(smoother, "(", vname, if(is.null(
				arg)) ")" else paste(",", arg, ")", sep = ""),
				sep = ""))
		if(form)
			junk <- eval(parse(text = paste("~", paste(junk, 
				collapse = "+"))))
		step.list[[vname]] <- junk
	}
	step.list
}
