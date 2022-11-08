#' @export
"gam.smooth.list"=list(
    slist=c("s","lo","random"),
    wlist=c("s","lo")
)


#' Smoothers available for backfitting
#' 
#' Auxiliary function as user interface for 'gam' fitting. Lists what smoothers
#' are implemented, and allows users to include new smoothers.
#' 
#' 
#' @aliases gam.smoothers gam.smooth.list
#' @param slist character vector giving names of smoothers available for
#' general backfitting. For every entry, eg "lo", there must exist a formula
#' function "lo()" that prepares the data, and a fitting function with the name
#' "gam.lo" which actually does the fitting. Look at "lo" and "s" as examples.
#' @param wlist character vector (subset of slist) giving names of smoothers
#' for which a special backfitting algorithm is available, when only that
#' smoother appears (multiple times) in the formula, along with other non
#' smooth terms.
#' @return a list is returned, consisting of the two named vectors. If the
#' function is called with no arguments, it gets the version of
#' "gam.smooth.list"' in the search path, by default from the package name
#' space. Once it is called with either of the arguments, it places a local
#' copy in the users namespace.
#' @references Hastie, T. J. (1992) \emph{Generalized additive models.} Chapter
#' 7 of \emph{Statistical Models in S} eds J. M. Chambers and T. J. Hastie,
#' Wadsworth & Brooks/Cole.
#' @keywords models regression nonparametric smooth
#' @examples
#' 
#' \dontrun{gam.smoothers()$slist # get the gam.smooth.list, and extract component slist}
#' \dontrun{gam.smoothers(slist=c("s","lo","random","tps") # add a new smoother "tps" to the list}
#' 
#' @export gam.smoothers
"gam.smoothers" <- function(slist=c("s","lo","random"), wlist=c("s","lo")){
    smooth.list=gam.smooth.list
    if(!missing(slist)){
        smooth.list$slist <- slist
        assignInMyNamespace("gam.smooth.list", smooth.list)
    }
    if(!missing(wlist)){
        smooth.list$wlist <- wlist
        assignInMyNamespace("gam.smooth.list", smooth.list)
    }
    smooth.list
    }

