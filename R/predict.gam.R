"predict.gam" <-
  function(object, newdata, type = c("link", "response", "terms"), dispersion=NULL, se.fit = FALSE, na.action=na.pass, terms = labels(object),...)
{
  type <- match.arg(type)
  if(missing(newdata)) {
    if(inherits(object, "gam") && !is.null(object$smooth)) {
      if(se.fit)
        switch(type,
               response = {
                 out <- predict.gam(object,
                                    type = "link", se.fit
                                    = TRUE, ...)
                 famob <- family(object)
                 out$se.fit <- drop(out$se.fit*abs(famob$mu.eta(out$fit)))
                 out$fit <- fitted(object)
                 out
               }
               ,
               link = {
                 out <- NextMethod("predict")
                 out$fit <- object$additive.predictors
                 TS <- out$residual.scale^2
                 TT <- ncol(object$var)
                 out$se.fit <- sqrt(out$se.fit^
                                    2 + TS * object$var %*%
                                    rep(1, TT))
                 out
               }
               ,
               terms = {
                 out <- NextMethod("predict")
                 TT <- dimnames(s <- object$smooth)[[2]]
                 out$fit[, TT] <- out$fit[,
                                          TT] + s
                 TS <- out$residual.scale^2
                 out$se.fit[, TT] <- sqrt(out$
                                          se.fit[, TT]^2 + TS *
                                          object$var)
                 out
               }
               )
      else switch(type,
                  terms = {
                    out <- NextMethod("predict")
                    TT <- dimnames(s <- object$smooth)[[2]]
                    out[, TT] <- out[, TT] + s
                    out
                  }
                  ,
                  link = object$additive.predictors,

                  response = object$fitted)
    }
    else {
      if(inherits(object, "gam")) {
        if(type == "link" && !se.fit)
          object$additive.predictors
        else NextMethod("predict")
      }
      else UseMethod("predict")
    }
  }
  else newdata.predict.gam(object, newdata, type, dispersion,se.fit, na.action, terms, ...)
}
