#' Fitting Generalized Additive Models
#' 
#' \code{gam} is used to fit generalized additive models, specified by giving a
#' symbolic description of the additive predictor and a description of the
#' error distribution. \code{gam} uses the \emph{backfitting algorithm} to
#' combine different smoothing or fitting methods. The methods currently
#' supported are local regression and smoothing splines.
#' 
#' The gam model is fit using the local scoring algorithm, which iteratively
#' fits weighted additive models by backfitting. The backfitting algorithm is a
#' Gauss-Seidel method for fitting additive models, by iteratively smoothing
#' partial residuals.  The algorithm separates the parametric from the
#' nonparametric part of the fit, and fits the parametric part using weighted
#' linear least squares within the backfitting algorithm. This version of
#' \code{gam} remains faithful to the philosophy of GAM models as outlined in
#' the references below.
#' 
#' An object \code{gam.slist} (currently set to \code{c("lo","s","random")})
#' lists the smoothers supported by \code{gam}. Corresponding to each of these
#' is a smoothing function \code{gam.lo}, \code{gam.s} etc that take particular
#' arguments and produce particular output, custom built to serve as building
#' blocks in the backfitting algorithm. This allows users to add their own
#' smoothing methods. See the documentation for these methods for further
#' information. In addition, the object \code{gam.wlist} (currently set to
#' \code{c("s","lo")}) lists the smoothers for which efficient backfitters are
#' provided. These are invoked if all the smoothing methods are of one kind
#' (either all \code{"lo"} or all \code{"s"}).
#' 
#' @aliases gam gam.fit
#' @param formula a formula expression as for other regression models, of the
#' form \code{response ~ predictors}. See the documentation of \code{lm} and
#' \code{formula} for details.  Built-in nonparametric smoothing terms are
#' indicated by \code{s} for smoothing splines or \code{lo} for \code{loess}
#' smooth terms.  See the documentation for \code{s} and \code{lo} for their
#' arguments. Additional smoothers can be added by creating the appropriate
#' interface functions. Interactions with nonparametric smooth terms are not
#' fully supported, but will not produce errors; they will simply produce the
#' usual parametric interaction.
#' @param family a description of the error distribution and link function to
#' be used in the model. This can be a character string naming a family
#' function, a family function or the result of a call to a family function.
#' (See \code{\link{family}} for details of family functions.)
#' @param data an optional data frame containing the variables in the model.
#' If not found in \code{data}, the variables are taken from
#' \code{environment(formula)}, typically the environment from which \code{gam}
#' is called.
#' @param weights an optional vector of weights to be used in the fitting
#' process.
#' @param subset an optional vector specifying a subset of observations to be
#' used in the fitting process.
#' @param na.action a function which indicates what should happen when the data
#' contain \code{NA}s.  The default is set by the \code{na.action} setting of
#' \code{\link{options}}, and is \code{\link{na.fail}} if that is unset.  The
#' \dQuote{factory-fresh} default is \code{\link{na.omit}}. A special method
#' \code{\link{na.gam.replace}} allows for mean-imputation of missing values
#' (assumes missing at random), and works gracefully with \code{gam}
#' @param start starting values for the parameters in the additive predictor.
#' @param etastart starting values for the additive predictor.
#' @param mustart starting values for the vector of means.
#' @param offset this can be used to specify an \emph{a priori} known component
#' to be included in the additive predictor during fitting.
#' @param control a list of parameters for controlling the fitting process.
#' See the documentation for \code{\link{gam.control}} for details. These can
#' also be set as arguments to \code{gam()} itself.
#' @param model a logical value indicating whether \emph{model frame} should be
#' included as a component of the returned value. Needed if \code{gam} is
#' called and predicted from inside a user function. Default is \code{TRUE}.
#' @param method the method to be used in fitting the parametric part of the
#' model.  The default method \code{"glm.fit"} uses iteratively reweighted
#' least squares (IWLS).  The only current alternative is \code{"model.frame"}
#' which returns the model frame and does no fitting.
#' @param x,y For \code{gam}: logical values indicating whether the response
#' vector and model matrix used in the fitting process should be returned as
#' components of the returned value.
#' 
#' For \code{gam.fit}: \code{x} is a model matrix of dimension \code{n * p},
#' and \code{y} is a vector of observations of length \code{n}.
#' @param smooth.frame for \code{gam.fit} only. This is essentially a subset of
#' the model frame corresponding to the smooth terms, and has the ingredients
#' needed for smoothing each variable in the backfitting algorithm. The
#' elements of this frame are produced by the formula functions \code{lo} and
#' \code{s}.
#' @param \dots further arguments passed to or from other methods.
#' @return \code{gam} returns an object of class \code{Gam}, which inherits
#' from both \code{glm} and \code{lm}.
#' 
#' Gam objects can be examined by \code{print}, \code{summary}, \code{plot},
#' and \code{anova}.  Components can be extracted using extractor functions
#' \code{predict}, \code{fitted}, \code{residuals}, \code{deviance},
#' \code{formula}, and \code{family}. Can be modified using \code{update}. It
#' has all the components of a \code{glm} object, with a few more. This also
#' means it can be queried, summarized etc by methods for \code{glm} and
#' \code{lm} objects. Other generic functions that have methods for \code{Gam}
#' objects are \code{step} and \code{preplot}.
#' 
#' The following components must be included in a legitimate `Gam' object. The
#' residuals, fitted values, coefficients and effects should be extracted by
#' the generic functions of the same name, rather than by the \code{"$"}
#' operator. The \code{family} function returns the entire family object used
#' in the fitting, and \code{deviance} can be used to extract the deviance of
#' the fit.
#' 
#' \item{coefficients}{ the coefficients of the parametric part of the
#' \code{additive.predictors}, which multiply the columns of the model matrix.
#' The names of the coefficients are the names of the single-degree-of-freedom
#' effects (the columns of the model matrix). If the model is overdetermined
#' there will be missing values in the coefficients corresponding to
#' inestimable coefficients. } \item{additive.predictors}{ the additive fit,
#' given by the product of the model matrix and the coefficients, plus the
#' columns of the \code{$smooth} component. } \item{fitted.values}{ the fitted
#' mean values, obtained by transforming the component
#' \code{additive.predictors} using the inverse link function. } \item{smooth,
#' nl.df, nl.chisq, var}{ these four characterize the nonparametric aspect of
#' the fit. \code{smooth} is a matrix of smooth terms, with a column
#' corresponding to each smooth term in the model; if no smooth terms are in
#' the \code{Gam} model, all these components will be missing. Each column
#' corresponds to the strictly nonparametric part of the term, while the
#' parametric part is obtained from the model matrix. \code{nl.df} is a vector
#' giving the approximate degrees of freedom for each column of \code{smooth}.
#' For smoothing splines specified by \code{s(x)}, the approximate \code{df}
#' will be the trace of the implicit smoother matrix minus 2. \code{nl.chisq}
#' is a vector containing a type of score test for the removal of each of the
#' columns of \code{smooth}. \code{var} is a matrix like \code{smooth},
#' containing the approximate pointwise variances for the columns of
#' \code{smooth}. } \item{smooth.frame}{This is essentially a subset of the
#' model frame corresponding to the smooth terms, and has the ingredients
#' needed for making predictions from a \code{Gam} object} \item{residuals}{
#' the residuals from the final weighted additive fit; also known as residuals,
#' these are typically not interpretable without rescaling by the weights. }
#' \item{deviance}{ up to a constant, minus twice the maximized log-likelihood.
#' Similar to the residual sum of squares. Where sensible, the constant is
#' chosen so that a saturated model has deviance zero. }
#' \item{null.deviance}{The deviance for the null model, comparable with
#' \code{deviance}. The null model will include the offset, and an intercept if
#' there is one in the model} \item{iter}{ the number of local scoring
#' iterations used to compute the estimates. } \item{bf.iter}{a vector of
#' length \code{iter} giving number of backfitting iterations used at each
#' inner loop.} \item{family}{ a three-element character vector giving the name
#' of the family, the link, and the variance function; mainly for printing
#' purposes. } \item{weights}{the \emph{working} weights, that is the weights
#' in the final iteration of the local scoring fit.} \item{prior.weights}{the
#' case weights initially supplied.} \item{df.residual}{the residual degrees of
#' freedom.} \item{df.null}{the residual degrees of freedom for the null
#' model.}
#' 
#' The object will also have the components of a \code{lm} object:
#' \code{coefficients}, \code{residuals}, \code{fitted.values}, \code{call},
#' \code{terms}, and some others involving the numerical fit.  See
#' \code{lm.object}.
#' @author Written by Trevor Hastie, following closely the design in the
#' "Generalized Additive Models" chapter (Hastie, 1992) in Chambers and Hastie
#' (1992), and the philosophy in Hastie and Tibshirani (1991).  This version of
#' \code{gam} is adapted from the S version to match the \code{glm} and
#' \code{lm} functions in R.
#' 
#' Note that this version of \code{gam} is different from the function with the
#' same name in the R library \code{mgcv}, which uses only smoothing splines
#' with a focus on automatic smoothing parameter selection via GCV. To avoid
#' issues with S3 method handling when both packages are loaded, the object
#' class in package "gam" is now "Gam".
#' 
#' @seealso \code{\link{glm}}, \code{\link{family}}, \code{\link{lm}}.
#' @references Hastie, T. J. (1991) \emph{Generalized additive models.} Chapter
#' 7 of \emph{Statistical Models in S} eds J. M. Chambers and T. J. Hastie,
#' Wadsworth & Brooks/Cole.
#' 
#' Hastie, T. and Tibshirani, R. (1990) \emph{Generalized Additive Models.}
#' London: Chapman and Hall.
#' 
#' Venables, W. N. and Ripley, B. D. (2002) \emph{Modern Applied Statistics
#' with S.} New York: Springer.
#' @keywords models regression nonparametric smooth
#' @examples
#' 
#' data(kyphosis)
#' gam(Kyphosis ~ s(Age,4) + Number, family = binomial, data=kyphosis,
#' trace=TRUE)
#' data(airquality)
#' gam(Ozone^(1/3) ~ lo(Solar.R) + lo(Wind, Temp), data=airquality, na=na.gam.replace)
#' gam(Kyphosis ~ poly(Age,2) + s(Start), data=kyphosis, family=binomial, subset=Number>2)
#' data(gam.data)
#' Gam.object <- gam(y ~ s(x,6) + z,data=gam.data)
#' summary(Gam.object)
#' plot(Gam.object,se=TRUE)
#' data(gam.newdata)
#' predict(Gam.object,type="terms",newdata=gam.newdata)
#' 
#' @export gam
"gam" <-
  function(formula, family = gaussian, data,
           weights, subset, na.action, start = NULL, etastart, mustart, control = gam.control(...),
           model = TRUE, method="glm.fit", x = FALSE, y = TRUE, ...)
{
  call <- match.call()
  if (is.character(family))
    family <- get(family, mode = "function", envir = parent.frame())
  if (is.function(family))
    family <- family()
  if (is.null(family$family)) {
    print(family)
    stop("`family' not recognized")
  }
  if (missing(data))
    data <- environment(formula)
  mf <- match.call(expand.dots = FALSE)
  ## m <- match(c("formula", "data", "subset", "weights", "na.action",
  ##              "etastart", "mustart", "offset"), names(mf), 0L)
 m <- match(c("formula", "data", "subset", "weights",
              "etastart", "mustart", "offset"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$na.action=quote(na.pass)## need to do this because model frame is not subsetting properly
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- quote(stats::model.frame)
  gam.slist <- gam.smoothers()$slist
  mt <- if(missing(data)) terms(formula, gam.slist) else terms(formula,gam.slist,data = data)
  mf$formula<-mt
  mf <- eval(mf, parent.frame())
  if(missing(na.action)){
      naa=getOption("na.action","na.fail")
      na.action=get(naa)
      }
  mf=na.action(mf)###because this was not done properly before
  mt=attributes(mf)[["terms"]]# the predvars are added here, while not before
    switch(method, model.frame = return(mf), glm.fit = 1, glm.fit.null = 1,
         stop("invalid `method': ", method))


  Y <- model.response(mf, "any")
  X <- if (!is.empty.model(mt))
###    model.matrix(mt, mf, contrasts) #not sure why the contrasts argument?
    model.matrix(mt, mf)
  else matrix(, NROW(Y), 0)
  weights <- model.weights(mf)
  offset <- model.offset(mf)
  if (!is.null(weights) && any(weights < 0))
    stop("Negative wts not allowed")
  if (!is.null(offset) && length(offset) != NROW(Y))
    stop("Number of offsets is ", length(offset), ", should equal ",
         NROW(Y), " (number of observations)")
  mustart <- model.extract(mf, "mustart")
  etastart <- model.extract(mf, "etastart")
fit<-gam.fit(x=X,y=Y,smooth.frame=mf,weights=weights,start=start,
             etastart=etastart,mustart=mustart,
             offset=offset,family=family,control=control)

### If both an offset and intercept are present, iterations are needed to
### compute the Null deviance; these are done here
###
  if(length(offset) && attr(mt, "intercept")>0) {
    fit$null.dev <- glm.fit(x = X[, "(Intercept)", drop = FALSE],
               y = Y, weights = weights, offset = offset, family = family,
               control = control[c("epsilon","maxit","trace")], intercept = TRUE)$deviance
  }
    if(model) fit$model <- mf
    fit$na.action <- attr(mf, "na.action")
    if(x) fit$x <- X
    if(!y) fit$y <- NULL
   fit <- c(fit, list(call = call, formula = formula,
		       terms = mt, data = data,
		       offset = offset, control = control, method = method,
		       contrasts = attr(X, "contrasts"),
                       xlevels = .getXlevels(mt, mf)))
    class(fit) <- c("Gam","glm", "lm")
  if(!is.null(fit$df.residual) && !(fit$df.residual > 0))
    warning("Residual degrees of freedom are negative or zero.  This occurs when the sum of the parametric and nonparametric degrees of freedom exceeds the number of observations.  The model is probably too complex for the amount of data available."
            )
  fit
}

