#' Stepwise model builder for GAM
#' 
#' Builds a GAM model in a step-wise fashion. For each "term" there is an
#' ordered list of alternatives, and the function traverses these in a greedy
#' fashion. Note: this is NOT a method for \code{step}, which used to be a
#' generic, so must be invoked with the full name.
#' 
#' @param object An object of class \code{Gam} or any of it's inheritants.
#' @param scope defines the range of models examined in the step-wise search.
#' It is a list of formulas, with each formula corresponding to a term in the
#' model. Each of these formulas specifies a "regimen" of candidate forms in
#' which the particular term may enter the model. For example, a term formula
#' might be \code{~1+ Income + log(Income) + s(Income)}. This means that
#' \code{Income} could either appear not at all, linearly, linearly in its
#' logarithm, or as a smooth function estimated nonparametrically. A \code{1}
#' in the formula allows the additional option of leaving the term out of the
#' model entirely.  Every term in the model is described by such a term
#' formula, and the final model is built up by selecting a component from each
#' formula.
#' 
#' As an alternative more convenient for big models, each list can have instead
#' of a formula a character vector corresponding to the candidates for that
#' term. Thus we could have \code{c("1","x","s(x,df=5")} rather than
#' \code{~1+x+s(x,df=5)}.
#' 
#' The supplied model \code{object} is used as the starting model, and hence
#' there is the requirement that one term from each of the term formulas be
#' present in \code{formula(object)}. This also implies that any terms in
#' \code{formula(object)} \emph{not} contained in any of the term formulas will
#' be forced to be present in every model considered. The function
#' \code{gam.scope} is helpful for generating the scope argument for a large
#' model.
#' @param scale an optional argument used in the definition of the AIC
#' statistic used to evaluate models for selection. By default, the scaled
#' Chi-squared statistic for the initial model is used, but if forward
#' selection is to be performed, this is not necessarily a sound choice.
#' @param direction The mode of step-wise search, can be one of \code{"both"},
#' \code{"backward"}, or \code{"forward"}, with a default of \code{"both"}. If
#' \code{scope} is missing, the default for \code{direction} is "both".
#' @param trace If \code{TRUE} (the default), information is printed during the
#' running of \code{step.Gam()}. This is an encouraging choice in general,
#' since \code{step.Gam()} can take some time to compute either for large
#' models or when called with an an extensive \code{scope=} argument. A simple
#' one line model summary is printed for each model selected. This argument can
#' also be given as the binary \code{0} or \code{1}. A value \code{trace=2}
#' gives a more verbose trace.
#' @param keep A filter function whose input is a fitted \code{Gam} object, and
#' anything else passed via \dots{}, and whose output is arbitrary. Typically
#' \code{keep()} will select a subset of the components of the object and
#' return them. The default is not to keep anything.
#' @param steps The maximum number of steps to be considered. The default is
#' 1000 (essentially as many as required). It is typically used to stop the
#' process early.
#' @param parallel If \code{TRUE}, use parallel \code{foreach} to fit each
#' trial run.  Must register parallel before hand, such as \code{doMC} or
#' others.  See the example below.
#' @param \dots Additional arguments to be passed on to \code{keep}
#' @return The step-wise-selected model is returned, with up to two additional
#' components. There is an \code{"anova"} component corresponding to the steps
#' taken in the search, as well as a \code{"keep"} component if the
#' \code{keep=} argument was supplied in the call.
#' 
#' We describe the most general setup, when \code{direction = "both"}. At any
#' stage there is a current model comprising a single term from each of the
#' term formulas supplied in the \code{scope=} argument. A series of models is
#' fitted, each corrresponding to a formula obtained by moving each of the
#' terms one step up or down in its regimen, relative to the formula of the
#' current model. If the current value for any term is at either of the extreme
#' ends of its regimen, only one rather than two steps can be considered. So if
#' there are \code{p} term formulas, at most \code{2*p - 1} models are
#' considered. A record is kept of all the models ever visited (hence the
#' \code{-1} above), to avoid repetition. Once each of these models has been
#' fit, the "best" model in terms of the AIC statistic is selected and defines
#' the step. The entire process is repeated until either the maximum number of
#' steps has been used, or until the AIC criterion can not be decreased by any
#' of the eligible steps.
#' @author Written by Trevor Hastie, following closely the design in the
#' "Generalized Additive Models" chapter (Hastie, 1992) in Chambers and Hastie
#' (1992).
#' @seealso \code{\link{gam.scope}},\code{\link{step}},\code{\link{glm}},
#' \code{\link{gam}}, \code{\link{drop1}}, \code{\link{add1}},
#' \code{\link{anova.Gam}}
#' @references Hastie, T. J. (1992) \emph{Generalized additive models.} Chapter
#' 7 of \emph{Statistical Models in S} eds J. M. Chambers and T. J. Hastie,
#' Wadsworth & Brooks/Cole.
#' 
#' Hastie, T. and Tibshirani, R. (1990) \emph{Generalized Additive Models.}
#' London: Chapman and Hall.
#' @keywords models regression nonparametric smooth
#' @examples
#' 
#' data(gam.data)
#' Gam.object <- gam(y~x+z, data=gam.data)
#' step.object <-step.Gam(Gam.object, scope=list("x"=~1+x+s(x,4)+s(x,6)+s(x,12),"z"=~1+z+s(z,4)))
#' \dontrun{
#' # Parallel
#' require(doMC)
#' registerDoMC(cores=2)
#' step.Gam(Gam.object, scope=list("x"=~1+x+s(x,4)+s(x,6)+s(x,12),"z"=~1+z+s(z,4)),parallel=TRUE)
#' }
#' 
#' 
#' @export step.Gam
step.Gam <-
  function (object, scope, scale, direction = c("both", "backward",
                                    "forward"), trace = TRUE, keep = NULL, steps = 1000, parallel=FALSE,...)
{
 trace=as.numeric(trace)
 get.visit <- function(trial, visited){
    match(paste(trial,collapse=""),apply(visited,2,paste,collapse=""),FALSE)
  }
deviancelm <-
function(object, ...)
if(is.null(w <- object$weights)) sum(object$residuals^2) else sum(w * object$
		residuals^2)
 scope.char <- function(formula) {
    formula = update(formula, ~-1 + .)
    tt <- terms(formula)
    tl <- attr(tt, "term.labels")
    if (attr(tt, "intercept"))
      c("1", tl)
    else tl
  }
  re.arrange <- function(keep) {
    namr <- names(k1 <- keep[[1]])
    namc <- names(keep)
    nc <- length(keep)
    nr <- length(k1)
    array(unlist(keep, recursive = FALSE), c(nr, nc), list(namr,
                         namc))
  }
  untangle.scope <- function(terms, regimens) {
    a <- attributes(terms)
    response <- deparse(a$variables[[2]])
    term.labels <- a$term.labels
    if (!is.null(a$offset)) {
      off1 <- deparse(a$variables[[a$offset]])
    }
    nt <- length(regimens)
    select <- integer(nt)
    for (i in seq(nt)) {
      j <- match(regimens[[i]], term.labels, 0)
      if (any(j)) {
        if (sum(j > 0) > 1)
          stop(paste("The elements of a regimen", i,
                     "appear more than once in the initial model",
                     sep = " "))
        select[i] <- seq(j)[j > 0]
        term.labels <- term.labels[-sum(j)]
      }
      else {
        if (!(j <- match("1", regimens[[i]], 0)))
          stop(paste("regimen", i, "does not appear in the initial model",
                     sep = " "))
        select[i] <- j
      }
    }
    if (length(term.labels))
      term.labels <- paste(term.labels, "+")
    if (!is.null(a$offset))
      term.labels <- paste(off1, term.labels, sep = " + ")
    return(list(response = paste(response, term.labels, sep = " ~ "),
                select = select))
  }
  make.step <- function(models, fit, scale, object) {
    chfrom <- sapply(models, "[[", "from")
    chfrom[chfrom == "1"] <- ""
    chto <- sapply(models, "[[", "to")
    chto[1]="<start>"
    chto[chto == "1"] <- ""
    dev <- sapply(models, "[[", "deviance")
    df <- sapply(models, "[[", "df.resid")
    ddev <- c(NA, diff(dev))
    ddf <- c(NA, diff(df))
    AIC <- sapply(models, "[[", "AIC")
    heading <- c("Stepwise Model Path \nAnalysis of Deviance Table",
                 "\nInitial Model:", deparse(as.vector(formula(object))),
                 "\nFinal Model:", deparse(as.vector(formula(fit))),
                 paste("\nScale: ", format(scale), "\n", sep = ""))
#    rowns=paste(chfrom,chto,sep=" -> ")
#    rowns[1]="<start>"
#    rowns=paste(seq(rowns)-1,rowns,sep=": ")
    aod <- data.frame(From=chfrom,To=chto, Df = ddf,
                      Deviance = ddev, `Resid. Df` = df, `Resid. Dev` = dev,
                      AIC = AIC, check.names = FALSE)
     aod <- as.anova(aod, heading)
     class(aod)=c("stepanova","data.frame")
     fit$anova=aod
    fit
  }
  direction <- match.arg(direction)
  if (missing(scope))
    stop("you must supply a scope argument to step.Gam(); the gam.scope() function might be useful")
  if (!is.character(scope[[1]]))
    scope <- lapply(scope, scope.char)
  response <- untangle.scope(object$terms, scope)
  form.y <- response$response
  backward <- direction == "both" | direction == "backward"
  forward <- direction == "both" | direction == "forward"
  items <- response$select
  family <- family(object)
  Call <- object$call
  term.lengths <- sapply(scope, length)
  n.items <- length(items)
  visited <- matrix(items)
  form.vector <- character(n.items)
  for (i in seq(n.items)) form.vector[i] <- scope[[i]][items[i]]
  form <- deparse(object$formula)
  if (trace>0)
    cat("Start: ", form)
  fit <- object
  n <- length(fit$fitted)
  if (missing(scale)) {
    famname <- family$family["name"]
    scale <- switch(famname, Poisson = 1, Binomial = 1, deviancelm(fit)/fit$df.resid)
  }
  else if (scale == 0)
    scale <- deviancelm(fit)/fit$df.resid
  bAIC <- fit$aic
  if (trace>0)
    cat("; AIC=", format(round(bAIC, 4)), "\n")
  models <- list(
                 list(deviance = deviance(fit), df.resid = fit$df.resid,
                      AIC = bAIC, from = "", to = "")
                 )
  if (!is.null(keep))   {
    keep.list <- list(keep(fit,...))
    keep.it=TRUE}
  else keep.it=FALSE
  AIC <- bAIC + 1
  stepnum=0
  while (bAIC < AIC & steps > 0) {
    steps <- steps - 1
    stepnum=stepnum+1
    AIC <- bAIC
    form.list=NULL
###First some prelimenaries to see what formulas to try
    for (i in seq(n.items)) {
      if (backward) {
        trial <- items
        trial[i] <- trial[i] - 1
        if (trial[i] > 0 && !get.visit(trial,visited)) {
          visited<-cbind(visited,trial)
          tform.vector <- form.vector
          tform.vector[i] <- scope[[i]][trial[i]]
          form.list=c(form.list,list(list(trial=trial, form.vector=tform.vector, which=i)))
        }
      }
      if (forward) {
        trial <- items
        trial[i] <- trial[i] + 1
        if (trial[i] <= term.lengths[i] && !get.visit(trial,visited)){
          visited<-cbind(visited,trial)
          tform.vector <- form.vector
          tform.vector[i] <- scope[[i]][trial[i]]
          form.list=c(form.list,list(list(trial=trial, form.vector=tform.vector, which=i)))
        }
      }
    }
    if(is.null(form.list))break
### Now we are ready for the expensive loop
### Parallel is set up
#if(parallel&&require(foreach)){
if(parallel){
#   step.list=foreach(i=1:length(form.list),.inorder=FALSE,.packages="gam",.verbose=trace>1)%dopar%
   step.list=foreach(i=1:length(form.list),.inorder=FALSE,.verbose=trace>1)%dopar%
    {
      tform=paste(form.y, paste(form.list[[i]]$form.vector, collapse = " + "))
      update(object, eval(parse(text = tform)),trace = FALSE, ...)
    }
  }
### No parallel
    else {
    step.list=as.list(sequence(length(form.list)))
    for(i in 1:length(form.list)){
      tform=paste(form.y, paste(form.list[[i]]$form.vector, collapse = " + "))
      step.list[[i]]=update(object, eval(parse(text = tform)),trace = FALSE, ...)
      if(trace>1)cat("Trial: ", tform,"; AIC=", format(round(step.list[[i]]$aic, 4)), "\n")
    }
}
### end expensive loop
    taic.vec=sapply(step.list,"[[","aic")
    if(keep.it)  keep.list=c(keep.list, lapply(step.list,keep,...))
    bAIC=min(taic.vec)
    if (bAIC >= AIC | steps == 0) {
      if (keep.it) fit$keep <- re.arrange(keep.list)
      return(make.step(models, fit, scale, object))
    }
    else {
      o1=order(taic.vec)[1]
      fit=step.list[[o1]]
      form.list=form.list[[o1]]
      bwhich=form.list$which
      bfrom=form.vector[bwhich]
      form.vector=form.list$form.vector #this is the new one
      bto=form.vector[bwhich]
      if (trace>0)
        cat(paste("Step:",stepnum,sep=""), deparse(fit$formula), "; AIC=",
            format(round(bAIC, 4)), "\n")
      items <- form.list$trial
      models <- c(models,list(list(deviance = deviance(fit), df.resid = fit$df.resid,
                           AIC = bAIC, from = bfrom, to = bto)))
    }
  }
}
