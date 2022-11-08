#' Generalized Additive Models
#'
#' This package provides functions for fitting and working with generalized additive models as described in chapter 7 of "Statistical Models in S" (Chambers and Hastie (eds), 1991) and "Generalized Additive Models" (Hastie and Tibshirani, 1990).
#' @name gam-package
#' @docType package
#' @author Trevor Hastie
#' @keywords models regression package
#' @useDynLib gam
#' @import methods stats splines foreach
#' @importFrom graphics axis lines mtext persp plot points rug segments
#' @importFrom utils head tail packageDescription menu assignInMyNamespace
NULL

#' Internal gam functions
#' 
#' @description Service functions and as yet undocumented functions for the gam library
#' @name gam-internal
#' @aliases .First.lib [.smooth general.wam anova.Gamlist as.anova
#' as.data.frame.lo.smooth assign.list Gamlist gam.match gam.nlchisq gam.sp
#' gplot gplot.default gplot.factor gplot.list gplot.matrix gplot.numeric
#' labels.Gam lo.wam newdata.predict.Gam polylo print.Gam print.Gamex
#' print.summary.Gam s.wam ylim.scale
#' @author Trevor Hastie
#' @keywords internal
NULL





#' Simulated dataset for gam
#' 
#' A simple simulated dataset, used to test out the gam functions
#' 
#' This dataset is artificial, and is used to test out some of the features of
#' gam.
#' 
#' @name gam.data
#' @aliases gam.data gam.newdata
#' @docType data
#' @format A data frame with 100 observations on the following 6 variables:
#' \describe{
#'   \item{x}{a numeric vector - predictor}
#'   \item{y}{a numeric vector - the response}
#'   \item{z}{a numeric vector - noise predictor}
#'   \item{f}{a numeric vector - true function}
#'   \item{probf}{a numeric vector - probability function}
#'   \item{ybin}{a numeric vector - binary response}
#' }
#' @keywords datasets
#' @examples
#' 
#' data(gam.data)
#' gam(y ~ s(x) + z, data=gam.data)
#' 
NULL





#' A classic example dataset for GAMs
#' 
#' Data on the results of a spinal operation "laminectomy" on children, to
#' correct for a condition called "kyphosis"; see Hastie and Tibshirani (1990)
#' for details
#' 
#' 
#' @name kyphosis
#' @docType data
#' @usage data(kyphosis)
#' @format A data frame with 81 observations on the following 4 variables.
#' \describe{ \item{Kyphosis}{a response factor with levels \code{absent}
#' \code{present}.} \item{Age}{of child in months, a numeric vector}
#' \item{Number}{of vertebra involved in the operation,a numeric vector}
#' \item{Start}{level of the operation, a numeric vector} }
#' @source Hastie, T. and Tibshirani, R. (1990) \emph{Generalized Additive
#' Models.} London: Chapman and Hall.
#' @keywords datasets
NULL


