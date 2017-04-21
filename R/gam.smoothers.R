"gam.smooth.list"=list(
    slist=c("s","lo","random"),
    wlist=c("s","lo")
)
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

