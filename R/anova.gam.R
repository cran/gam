#' Analysis of Deviance for a Generalized Additive Model
#' 
#' Produces an ANODEV table for a set of GAM models, or else a summary for a single GAM model
#'
#' These are methods for the functions \code{anova} or \code{summary} for
#' objects inheriting from class `Gam`. See \code{\link{anova}} for the general
#' behavior of this function and for the interpretation of `test`.
#'
#' When called with a single `Gam` object, a special pair of anova tables for
#' `Gam` models is returned. This gives a breakdown of the degrees of freedom
#' for all the terms in the model, separating the projection part and
#' nonparametric part of each, and returned as a list of two anova objects. For
#' example, a term specified by `s()` is broken down into a single degree of
#' freedom for its linear component, and the remainder for the nonparametric
#' component. In addition, a type of score test is performed for each of the
#' nonparametric terms. The nonparametric component is set to zero, and the
#' linear part is updated, holding the other nonparametric terms fixed. This is
#' done efficiently and simulataneously for all terms.
#' 
#' @aliases anova.Gam summary.Gam
#' @param object a fitted Gam
#' @param ... other fitted Gams for \code{anova}
#' @param test a character string specifying the test statistic to be used.
#' Can be one of '"F"', '"Chisq"' or '"Cp"', with partial matching allowed, or
#' 'NULL' for no test.
#' @param dispersion a dispersion parameter to be used in computing standard
#' errors
#' @author Written by Trevor Hastie, following closely the design in the
#' "Generalized Additive Models" chapter (Hastie, 1992) in Chambers and Hastie
#' (1992).
#' @references Hastie, T. J. (1992) \emph{Generalized additive models.} Chapter
#' 7 of \emph{Statistical Models in S} eds J. M. Chambers and T. J. Hastie,
#' Wadsworth & Brooks/Cole.
#' 
#' Hastie, T. and Tibshirani, R. (1990) \emph{Generalized Additive Models.}
#' London: Chapman and Hall.
#' 
#' Venables, W. N. and Ripley, B. D. (2002) \emph{Modern Applied Statistics
#' with S.} New York: Springer.
#' @keywords models regression nonparametric smooth
#' @method anova Gam
#' @export
#' @export anova.Gam
#' @examples
#' 
#' data(gam.data)
#' Gam.object <- gam(y~s(x,6)+z,data=gam.data)
#' anova(Gam.object)
#' Gam.object2 <- update(Gam.object, ~.-z)
#' anova(Gam.object, Gam.object2, test="Chisq")
"anova.Gam" <-
  function(object, ..., test = c("Chisq", "F", "Cp"))
{
  test=match.arg(test)
  margs <- function(...)
    nargs()
  if(margs(...))
    anova(structure(list(object, ...),class="glmlist"), test = test)
  else summary.Gam(object)$anova
}
