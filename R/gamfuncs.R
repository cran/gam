".First.lib" <-
function (libname, pkgname, where) 
{
library.dynam("gam",pkgname,libname)
}
"[.smooth" <-
function(x, ..., drop = FALSE)
{
	cl <- oldClass(x)
	oldClass(x) <- NULL
	ats <- attributes(x)
	ats$dimnames <- NULL
	ats$dim <- NULL
	ats$names <- NULL
	y <- x[..., drop = drop]
	if(!is.null(nas <- ats$NAs)) {
		if(is.null(d <- dim(x)))
			d <- c(length(x), 1.)
		navec <- array(logical(d[1.]), d)
		navec[nas,  ] <- TRUE
		navec <- navec[...]
		nas <- if(is.null(dim(navec))) navec else navec[, 1.]
		nas <- seq(nas)[nas]
		ats$NAs <- nas
	}
	attributes(y) <- c(attributes(y), ats)
	oldClass(y) <- cl
	y
}
"all.wam" <-
  function(x, y, w, s, which, smooth.frame, maxit = 30, tol = 1e-7, trace = FALSE,
           se = TRUE, ...)
{
  if(inherits(smooth.frame, "data.frame")) {
    data <- smooth.frame
### Note; the lev component of the smooths is the diagonal hat matrix elements
### for the NONLINEAR part of the fit.
###The smoother can return both the linear and nonlinear parts, although only
### the nonlinear part is strictly necessary. 
###
    oldClass(data) <- NULL
    names.calls <- names(which)
    smooth.calls <- lapply(data[names.calls], attr, "call")
    names(smooth.calls) <- names.calls
    smooth.frame <- list(data = data, smooth.calls = smooth.calls)
  }
  else {
    data <- smooth.frame$data
    smooth.calls <- smooth.frame$smooth.calls
  }
  names.calls <- names(smooth.calls)
  y <- as.vector(y)
  residuals <- as.vector(y - s %*% rep(1., ncol(s)))
  n <- length(y)
  fit <- list(fitted.values = 0.)
  rss <- weighted.mean(residuals^2., w)
  rssold <- rss * 10.
  nit <- 0.
  df <- rep(NA, length(which))
  var <- s
  if(trace)
    cat("\nWAM   iter   rss/n     term\n")
  ndig <-  - log10(tol) + 1.
  RATIO <- tol + 1.
  while(RATIO > tol & nit < maxit) {
    rssold <- rss
    nit <- nit + 1.
    z <- residuals + fit$fitted.values
    fit <- lm.wfit(x, z, w, method = "qr", singular.ok = TRUE,
                   ...)
    residuals <- fit$residuals
    rss <- weighted.mean(residuals^2., w)
    if(trace)
      cat("\n         ", nit, "   ", format(round(rss, ndig)),
          "  Parametric -- lm.wfit\n", sep = "")
    deltaf <- 0.
    for(j in seq(names.calls)) {
      old <- s[, j]
      z <- residuals + s[, j]
      fit.call <- eval(smooth.calls[[j]])
      residuals <- as.double(fit.call$residuals)
      if(length(residuals) != n)
        stop(paste(names.calls[j], 
                   "returns a vector of the wrong length")
             )
      s[, j] <- z - residuals
      deltaf <- deltaf + weighted.mean((s[, j] - old)^2.,
                                       w)
      rss <- weighted.mean(residuals^2., w)
      if(trace) {
        cat("         ", nit, "   ", format(round(
                                                  rss, ndig)), "  Nonparametric -- ",
            names.calls[j], "\n", sep = "")
      }
      df[j] <- fit.call$nl.df
      if(se)
        var[, j] <- fit.call$var
    }
    RATIO <- sqrt(deltaf/sum(w * apply(s, 1., sum)^2.))
    if(trace)
      cat("Relative change in functions:", format(round(
                                                        RATIO, ndig)), "\n")
  }
  if((nit == maxit) & maxit > 1.)
    warning(paste("all.wam convergence not obtained in ", maxit,
                  " iterations"))
  names(df) <- names.calls
  if(trace)
    cat("\n")
  fit$fitted.values <- y - residuals
  rl <- c(fit, list(smooth = s, nl.df = df))
  rl$df.residual <- rl$df.residual - sum(df)
  if(se)
    rl <- c(rl, list(var = var))
  c(list(smooth.frame = smooth.frame), rl)
}
"anova.gam" <-
  function(object, ..., test = c("Chisq", "F", "Cp"))
{
  margs <- function(...)
    nargs()
  if(margs(...))
    anova.glmlist(list(object, ...), test = test)
  else summary.gam(object)$anova
}
"anova.gamlist" <-
function(object, ..., test = c("none", "Chisq", "F", "Cp"))
anova.glmlist(object, test = test)
"as.anova" <-
  function(df, heading)
{
  if(!inherits(df, "data.frame"))
    stop("df must be a data frame")
  attr(df, "heading") <- heading
                                        #if the "class" attribute of df already starts with "anova" return(df)
  if(inherits(df, "anova")) {
    dfClasses <- attr(df, "class")
    if(dfClasses[1] == "anova")
      return(df)
  }
  class(df) <- unique(c("anova", class(df)))
  df
}
"as.data.frame.lo.smooth" <-
function(x, row.names = NULL, optional = FALSE,...)
{
	d <- dim(x)
	nrows <- d[[1.]]
	dn <- dimnames(x)
	row.names <- dn[[1.]]
	value <- list(x)
	if(length(row.names)) {
		row.names <- as.character(row.names)
		if(length(row.names) != nrows)
			stop(paste("supplied", length(row.names), 
				"names for a data frame with", nrows, "rows"))
	}
	else if(optional)
		row.names <- character(nrows)
	else row.names <- as.character(seq(length = nrows))
	if(!optional)
		names(value) <- deparse(substitute(x))[[1.]]
	attr(value, "row.names") <- row.names
	oldClass(value) <- "data.frame"
	value
}
assign.list<-function(assignx,term.labels){
  ass<-as.list(seq(term.labels))
  names(ass)<-term.labels
  indexset<-seq(along=assignx)
  lapply(ass,function(i,indexset,assignx)indexset[assignx==i],indexset,assignx)
}
"deviance.default" <-
function(object, ...)
object$deviance
"deviance.glm" <-
function(object, ...)
object$deviance
"deviance.lm" <-
function(object, ...)
if(is.null(w <- object$weights)) sum(object$residuals^2.) else sum(w * object$
		residuals^2.)
"gam" <-
  function(formula, family = gaussian, data, 
           weights, subset, na.action, start = NULL, etastart, mustart, control = gam.control(...),
           model = FALSE, method="glm.fit", x = FALSE, y = TRUE, ...)
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
  m <- match(c("formula", "data", "subset", "weights", "na.action", 
               "etastart", "mustart", "offset"), names(mf), 0)
  mf <- mf[c(1, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1]] <- as.name("model.frame")
  mt <- if(missing(data)) terms(formula, gam.slist) else terms(formula,gam.slist,data = data)
  mf$formula<-mt                                                          
  mf <- eval(mf, parent.frame())
   switch(method, model.frame = return(mf), glm.fit = 1, glm.fit.null = 1, 
         stop("invalid `method': ", method))


  Y <- model.response(mf, "numeric")
  X <- if (!is.empty.model(mt)) 
    model.matrix(mt, mf, contrasts)
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
  if(any(offset) && attr(mt, "intercept")>0) {
    fit$null.dev <- glm.fit(x = X[, "(Intercept)", drop = FALSE], 
               y = Y, weights = weights, offset = offset, family = family, 
               control = control, intercept = TRUE)$deviance
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
    class(fit) <- c("gam","glm", "lm")
  if(!is.null(fit$df.residual) && !(fit$df.residual > 0))
    warning("Residual degrees of freedom are negative or zero.  This occurs when the sum of the parametric and nonparametric degrees of freedom exceeds the number of observations.  The model is probably too complex for the amount of data available."
            )
  fit
}

"gam.control" <-
function(epsilon = 9.9999999999999995e-08, bf.epsilon = 9.9999999999999995e-08,
	maxit = 30, bf.maxit = 30, trace = FALSE, ...)
{
	if(epsilon <= 0) {
		warning("the value of epsilon supplied is zero or negative; the default value of 1e-7 was used instead"
			)
		epsilon <- 9.9999999999999995e-08
	}
	if(maxit < 1) {
		warning("the value of maxit supplied is too small; the default value of 30 was used instead"
			)
		maxit <- 30
	}
	if(bf.epsilon <= 0) {
		warning("the value of bf.epsilon supplied is zero or negative; the default value of 1e-7 was used instead"
			)
		bf.epsilon <- 9.9999999999999995e-08
	}
	if(bf.maxit < 1) {
		warning("the value of bf.maxit supplied is too small; the default value of 30 was used instead"
			)
		bf.maxit <- 30
	}
	list(epsilon = epsilon, maxit = maxit, bf.epsilon = bf.epsilon, 
		bf.maxit = bf.maxit, trace = as.logical(trace)[1])
}
"gam.exact" <-
function(gam.obj)
### -----------------------------------------------------------------------------------
### gam.exact is a method for the gam class.
###  
### Computes the asymptotically exact variance-covariance matrix for the linear
### terms in the model (except for the intercept).
###
### Note: Use of lo in the model formula is not allowed.  
###
### Author: Aidan McDermott (AMcD)
### Date:   Mar 5, 2003
###
###         Mar 28, 2003
###         Fixed single linear term models -- thanks to Tim Ramsay
###         April 17, 2006
###         Modified to work in R by Trevor Hastie
###
### See:
### 
### [1] Issues in Semiparametric Regression: A Case Study of Time Series
###     Models in Air Pollution and Mortality,
###     Dominici F., McDermott A., Hastie T.J.,
###     Technical Report, Department of Biostatistics, Johns Hopkins University,
###     Baltimore, MD, USA. 
###  
### -----------------------------------------------------------------------------------
  {

    if ( is.na(match("gam",class(gam.obj))) ) {
      stop("not a gam object")
    }
    
    nl.df    <- gam.obj$nl.df
    terms    <- terms(gam.obj)
    at.terms <- attributes(terms)

    coef <- coef(gam.obj)
    
    w   <- gam.obj$weights
    mu  <- gam.obj$fitted.values
    eta <- gam.obj$additive.predictors
    y   <- as.matrix(gam.obj$y)

    family   <- family(gam.obj)
    mu.eta.val <- family$mu.eta(eta)
    z <- eta + (y - mu)/mu.eta.val

    
### Don't want lo in gam formula.
###    if ( length((at.terms$specials)$lo) > 0 ) {
###      stop("lo found in gam formula.")
###    }

    X   <- model.matrix(gam.obj)
    Y   <- as.matrix(gam.obj$y)

### only take terms that survived the original gam call
    names.coef <- names(coef)
    has.intercept <- match("(Intercept)",names.coef)
    if ( !is.na(has.intercept) ) names.coef <- names.coef[-has.intercept]
    X   <- X[,names.coef]
    tnames <- dimnames(X)[[2]]
    form   <- "y~"
    special.list <- c()
### Replace the df with the actual df returned by gam.
### Rewrite fromula to match names in X
    for ( k in 1:length(tnames) ) {
      if ( substring(tnames[k],1,2) == "s(" ) {
        s.call     <- match.call(s,parse(text=tnames[k]))
        this.name  <- as.name(paste("x",k,sep=""))

        which      <- match(tnames[k],names(nl.df))
        if ( is.na(which) ) stop(paste("can't find df for term",tnames[k]))
        this.df    <- nl.df[which]+1

        form <- paste(form,
                      "+s(",this.name,",df =",this.df,")")
        special.list <- c(special.list,k)
      }
      else if ( substring(tnames[k],1,3) == "lo(" ) {
        lname <- length(tnames[k])
        if ( substring(tnames[k],lname,lname) == "1" ) tnames[k] <- substring(tnames[k],1,(lname-1))
        if ( substring(tnames[k],lname,lname) == ")" ) {
        lo.call    <- match.call(lo,parse(text=tnames[k]))
        this.name  <- as.name(paste("x",k,sep=""))

        lo.call[[2]] <- this.name
        lo.call <- deparse(lo.call)
        form <- paste(form,"+",lo.call)
      }
        special.list <- c(special.list,k)
      }
      else form <- paste(form,"+x",k,sep="")
    }
    mydat <- data.frame(cbind(Y,X))
    names(mydat) <- c("y",paste("x",1:ncol(X),sep=""))

    XX <- X
    mydat[,"w"] <- w

    Control <- gam.obj$call$control
    if ( is.null(Control) ) {
      call      <- gam.obj$call
      call[[1]] <- as.name("gam.control")
      Control   <- eval(call,sys.parent())
    }

    for ( k in 1:length(tnames) ) {
      if ( substring(tnames[k],1,2) != "s("  & substring(tnames[k],1,3) != "lo(" ) {
        this.var <- paste("x",k,sep="")
        upd.form <- update(as.formula(form),paste(this.var,"~. -",this.var))

        XX[,k] <- gam(formula=upd.form,data=mydat,family=gaussian,weight=w,
                      control=eval(Control))$fitted
      }
    }

### Need to test we get some data
    if ( length(X) == 0 ) stop("nothing to do")
    
    X   <- X[,-special.list,drop=FALSE]
    sx  <- XX[,-special.list,drop=FALSE]
    swx <- w*sx
    
    if ( length(X) == 0 ) stop("no linear terms in the model -- nothing to do")
    
    A <- t(X) %*% ( w * X ) - t(X) %*% ( w * sx )
    B <- t(X*w) - t(swx)
    H <- solve(A) %*% B

    beta    <- H %*% z
    varbeta <- (H * (1/w)) %*% t(H) * as.vector(summary(gam.obj)$dispersion)
    se      <- sqrt(diag(varbeta))

    coef <- cbind(summary.glm(gam.obj)$coef,NA,NA,NA)
    tab <- cbind(beta,se,beta/se,2*(1-pnorm(beta/se)))
    coef[dimnames(tab)[[1]],c(5,6,7)] <- tab[,c(2,3,4)]

    dimnames(coef) <- list(dimnames(coef)[[1]],
                           c(dimnames(coef)[[2]][1:4],
                             "A-exact SE","A-exact Z","A-exact P")) 

    out.object <- list(coefficients=coef,covariance=varbeta)
    class(out.object) <- c("gamex")
    
    return(out.object)
  }

"gam.fit" <-
  function (x, y, smooth.frame, weights = rep(1, nobs), start = NULL, 
            etastart = NULL, mustart = NULL, offset = rep(0, nobs), family = gaussian(), 
            control = gam.control()) 
{
  ynames <- if (is.matrix(y)) 
    dimnames(y)[[1]]
  else names(y)
  xnames <- dimnames(x)[[2]]
  nobs <- NROW(y)
  nvars <- ncol(x)
  maxit <- control$maxit
  bf.maxit <- control$bf.maxit
  epsilon <- control$epsilon
  bf.epsilon <- control$bf.epsilon
  trace <- control$trace
  digits <- -log10(epsilon) + 1
  if (is.null(weights)) 
    weights <- rep.int(1, nobs)
  if (is.null(offset)) 
    offset <- rep.int(0, nobs)
  variance <- family$variance
  dev.resids <- family$dev.resids
  aic <- family$aic
  linkinv <- family$linkinv
  mu.eta <- family$mu.eta
  if (!is.function(variance) || !is.function(linkinv)) 
    stop("illegal `family' argument")
  valideta <- family$valideta
  if (is.null(valideta)) 
    valideta <- function(eta) TRUE
  validmu <- family$validmu
  if (is.null(validmu)) 
    validmu <- function(mu) TRUE
  eval(family$initialize)
  if (is.null(mustart)) {
    eval(family$initialize)
  }
  else {
    mukeep <- mustart
    eval(family$initialize)
    mustart <- mukeep
  }
  eta <- if (!is.null(etastart)) 
    etastart
  else if (!is.null(start)) 
    if (length(start) != nvars) 
      stop("Length of start should equal ", nvars, " and correspond to initial coefs for ", 
           deparse(xnames))
    else {
      coefold <- start
      offset + as.vector(if (NCOL(x) == 1) 
                         x * start
      else x %*% start)
    }
  else family$linkfun(mustart)
  mu <- linkinv(eta)
  if (!(validmu(mu) && valideta(eta))) 
    stop("Can't find valid starting values: please specify some")
  new.dev <- sum(dev.resids(y, mu, weights))
  a <- attributes(attr(smooth.frame, "terms"))
  smoothers <- a$specials
  if (length(smoothers) > 0) {
    smoothers <- smoothers[sapply(smoothers, length) > 0]
    for (i in seq(along = smoothers)) {
      tt <- smoothers[[i]]
      ff <- apply(a$factors[tt, , drop = FALSE], 2, any)
      smoothers[[i]] <- if (any(ff)) 
        seq(along = ff)[a$order == 1 & ff]
      else NULL
    }
  }
  if (length(smoothers) > 0) {
    smooth.labels <- a$term.labels[unlist(smoothers)]
    assignx <- attr(x, "assign")
    assignx <- assign.list(assignx, a$term.labels)
    which <- assignx[smooth.labels]
    if (length(smoothers) > 1) 
      bf <- "all.wam"
    else {
      sbf <- match(names(smoothers), gam.wlist, FALSE)
      bf <- if (sbf) 
        paste(gam.wlist[sbf], "wam", sep = ".")
      else "all.wam"
    }
    bf.call <- parse(text = paste(bf, "(x, z, wz, fit$smooth, which, fit$smooth.frame,bf.maxit,bf.epsilon, trace)", 
                       sep = ""))[[1]]
    s <- matrix(0, length(y), length(which))
    dimnames(s) <- list(names(y), names(which))
    fit <- list(smooth = s, smooth.frame = smooth.frame)
  }
  else {
    bf.call <- expression(lm.wfit(x, z, wz, method = "qr", 
        singular.ok = TRUE))
    bf <- "lm.wfit"
  }
  old.dev <- 10 * new.dev + 10
  for (iter in 1:maxit) {
    good <- weights > 0
    varmu <- variance(mu)
    if (any(is.na(varmu[good]))) 
      stop("NAs in V(mu)")
    if (any(varmu[good] == 0)) 
      stop("0s in V(mu)")
    mu.eta.val <- mu.eta(eta)
    if (any(is.na(mu.eta.val[good]))) 
      stop("NAs in d(mu)/d(eta)")
    good <- (weights > 0) & (mu.eta.val != 0)
    z <- eta - offset
    z[good] <- z[good] + (y - mu)[good]/mu.eta.val[good]
    wz <- weights
    wz[!good] <- 0
    wz[good] <- wz[good] * mu.eta.val[good]^2/varmu[good]
    fit <- eval(bf.call)
    eta <- fit$fitted.values + offset
    mu <- linkinv(eta)
    old.dev <- new.dev
    new.dev <- sum(dev.resids(y, mu, weights))
    if (trace) 
      cat("GAM ", bf, " loop ", iter, ": deviance = ", 
          format(round(new.dev, digits)), " \n", sep = "")
    if (is.na(new.dev)) {
      one.more <- FALSE
      warning("iterations terminated prematurely because of singularities")
    }
    else one.more <- abs(old.dev - new.dev)/(old.dev + 0.1) > 
      epsilon
    if (!one.more) 
      break
  }
  fitqr <- fit$qr
  xxnames <- xnames[fitqr$pivot]
  nr <- min(sum(good), nvars)
  if (nr < nvars) {
    Rmat <- diag(nvars)
    Rmat[1:nr, 1:nvars] <- fitqr$qr[1:nr, 1:nvars]
  }
  else Rmat <- fitqr$qr[1:nvars, 1:nvars]
  Rmat <- as.matrix(Rmat)
  Rmat[row(Rmat) > col(Rmat)] <- 0
  dimnames(Rmat) <- list(xxnames, xxnames)
  names(fit$residuals) <- ynames
  names(mu) <- ynames
  names(eta) <- ynames
  fit$additive.predictors <- eta
  fit$fitted.values <- mu
  names(fit$weights) <- ynames
  names(fit$effects) <- c(xxnames[seq(len = fitqr$rank)], rep.int("", 
                                        sum(good) - fitqr$rank))
  if (length(fit$smooth) > 0) 
    fit$smooth.frame <- smooth.frame[smooth.labels]
  wtdmu <- if (a$intercept) 
    sum(weights * y)/sum(weights)
  else linkinv(offset)
  nulldev <- sum(dev.resids(y, wtdmu, weights))
  n.ok <- nobs - sum(weights == 0)
  nulldf <- n.ok - as.integer(a$intercept)
  rank <- n.ok - fit$df.residual
  aic.model <- aic(y, n, mu, weights, new.dev) + 2 * rank
  if (!is.null(fit$smooth)) {
    nonzeroWt <- (wz > 0)
    nl.chisq <-  gam.nlchisq(fit$qr, fit$residuals, wz, fit$smooth)
  }
  else nl.chisq <- NULL
  fit <- c(fit, list(R = Rmat, rank = fitqr$rank, family = family, 
                     deviance = new.dev, aic = aic.model, null.deviance = nulldev, 
                     iter = iter, prior.weights = weights, y = y, df.null = nulldf, 
                     nl.chisq = nl.chisq))
  fit
}

"gamlist" <-
function(...)
{
	gl <- list(...)
	oldClass(gl) <- c("gamlist", "glmlist")
	gl
}
"gam.lo" <-
function(x, y, w = rep(1, length(y)), span = 0.5, degree = 1, ncols = p, xeval
	 = x)
{
	storage.mode(x) <- storage.mode(y) <- storage.mode(w) <- storage.mode(
		span) <- "double"
	storage.mode(degree) <- "integer"
	if(is.null(np <- dim(x))) {
		n <- as.integer(length(x))
		p <- as.integer(1)
	}
	else {
		np <- as.integer(np)
		n <- np[1]
		p <- np[2]
	}
	storage.mode(ncols) <- "integer"
	o <- gam.match(x)
	nef <- o$nef
	nvmax <- as.integer(200 + 300 * (1 - 1/log(max(c(nef - 200, 3)))))
	liv <- as.integer(50 + (2^ncols + 4) * nvmax + 2 * nef)
	lv <- as.integer(50 + (3 * ncols + 3) * nvmax + nef + (ifelse(degree ==
		2, ((ncols + 2) * (ncols + 1))/2, ncols + 1) + 2) * (nef * span +
		1))
	fit <- .Fortran("lo0",
		x,
		y,
		w,
		n,
		ncols,
		p,
		nvmax,
		span,
		degree,
		o$o,
		nef,
		df = double(1),
		s = double(n),
		var = double(n),
		beta = double(p + 1),
		iv = integer(liv),
		liv,
		lv,
		v = double(lv),
                integer(2*ncols),
		double(nef * (p + ncols + 8) + 2 * p + n + 9),
                        PACKAGE="gam")
	if(!missing(xeval)) {
		storage.mode(xeval) <- "double"
		m <- as.integer(dim(xeval)[1])
		if(length(m) == 0)
			m <- as.integer(length(xeval))
		.Fortran("lowese",
			fit$iv,
			liv,
			lv,
			fit$v,
			m,
			xeval,
			s = double(m),
                         PACKAGE="gam")$s - cbind(1, xeval) %*% fit$beta
	}
	else list(residuals = y - fit$s, var = fit$var, nl.df = fit$df)
}
"gam.match" <-
function(x)
{
	if(is.list(x)) {
		junk <- Recall(x[[1]])
		if((nvar <- length(x)) == 1)
			return(list(o = junk$o, nef = junk$nef))
		else {
			o <- matrix(junk$o, length(junk$o), nvar)
			nef <- rep(junk$nef, nvar)
			for(i in 2:nvar) {
				junk <- Recall(x[[i]])
				o[, i] <- junk$o
				nef[i] <- junk$nef
			}
			names(nef) <- nn <- names(x)
			dimnames(o) <- list(NULL, nn)
			return(list(o = o, nef = nef))
		}
	}
	if(is.matrix(x)) {
		ats <- attributes(x)
		a <- ats$NAs
		ncols <- ats$ncols
		d <- dim(x)
		if(is.null(ncols))
			ncols <- d[2]
		if(ncols == 1)
			return(Recall(structure(x[, 1, drop = TRUE], NAs = a)))
		if(is.null(a)) {
			o <- seq(d[1])
			nef <- d[1]
		}
		else {
			nef <- d[1] - length(a)
			o <- rep(nef + 1, d[1])
			o[ - a] <- seq(nef)
		}
		return(list(o = as.integer(o), nef = as.integer(nef)))
	}
	else {
		a <- attributes(x)$NAs
		if(!is.null(a))
			x[a] <- NA
		xr <- signif(as.vector(x), 6)
		sx <- unique(sort(xr))
		nef <- as.integer(length(sx))
		if(nef <= 3)
			stop("A smoothing variable encountered with 3 or less unique values; at least 4 needed"
				)
		o <- match(xr, sx, nef + 1)
		o[is.na(o)] <- nef + 1
		return(list(o = as.integer(o), nef = as.integer(nef)))
	}
}
"gam.nlchisq" <-
function(qr, resid, w, s)
{
	wt <- sqrt(w)
	s <- s * wt
	resid <- wt * resid
	Rsw <- qr.resid(qr, s)
	apply(Rsw^2 + 2 * s * resid, 2, sum)
}
"gam.random" <-
function(x, y, w, df = sum(non.zero), sigma = 0)
{
	df.inv <- function(n, df, sigma = sum(n)/df - mean(n), iterations = 10
		)
	{
		if(df > length(n))
			return(0)
		current.df <- sum(n/(n + sigma))
		if(abs((df - current.df)/df) < 0.0001 | iterations == 1)
			sigma
		else {
			sigma <- exp(log(sigma) + (current.df - df)/(sum((
				n * sigma)/(n + sigma)^2)))
			Recall(n, df, sigma, iterations - 1)
		}
	}
	nw <- tapply(w, x, sum)
	non.zero <- !is.na(nw)
	if(is.null(df))
		df <- sum(non.zero)
	if(sigma == 0)
		sigma <- df.inv(nw[non.zero], df)
	df <- sum(nw[non.zero]/(nw[non.zero] + sigma))
	fit <- tapply(w * y, x, sum)/(nw + sigma)
	var <- as.vector(w/(nw[x] + sigma))
	residuals <- as.vector(y - fit[x])
	list(x = seq(along = nw), y = fit, residuals = residuals, var = var,
		nl.df = df, sigma = sigma)
}
"gam.s" <-
  function(x, y, w = rep(1, length(x)), df = 4, spar = 1, xeval)
{
  storage.mode(x) <- storage.mode(y) <- storage.mode(w) <- storage.mode(
                                                                        spar) <- storage.mode(df) <- "double"
  n <- as.integer(length(x))
  x <- signif(x, 6)
  mat <- gam.match(x)
  omat <- mat$o
  nef <- mat$nef
  ##
  ## in rgam.r, splsm calls both splsm1 and splsm2.
  ## splsm2 needs (10+2*4)*(nef+2)+5*nef+n+15 doubles for work.
  ## splsm1 needs 3*nef+2*n+10.
  work.len <- max(3 * nef + 2 * n + 10, (10 + 2 * 4) * (nef + 2) + 5 *
                  nef + n + 15)
  fit <- .Fortran("splsm",
                  x,
                  y,
                  w,
                  n,
                  omat,
                  nef,
                  spar = spar,
                  df = df,
                  s = double(n),
                  s0 = double(1),
                  var = double(nef),
                  FALSE,
                  work = double(work.len),
                  PACKAGE="gam")
  if(missing(xeval))
    list(residuals = y - fit$s, nl.df = fit$df - 1, var = fit$
         var[omat])
  else {
    skn <- .Fortran("sknotl",
                    fit$work[seq(nef)],
                    nef,
                    knot = double(nef + 6),
                    k = integer(1),
                    PACKAGE="gam")
    smallest <- x[omat == 1][1]
    largest <- x[omat == nef][1]
    k <- skn$k
    gam.sp(xeval, skn$knot[seq(k)], k - 4, fit$work[seq(3 * nef +
                                                        n + 10, length = k - 4)], smallest, largest - smallest)
  }
}
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
"gam.slist" <-
c("s", "lo", "random")
"gam.sp" <-
function(x, knots, nknots, coef, smallest, scale)
{
	nas <- is.na(x)
	xs <- as.double((x[!nas] - smallest)/scale)
	bad.left <- xs < 0
	bad.right <- xs > 1
	good <- !(bad.left | bad.right)
	y <- xs
	if(any(good)) {
		junk <- .Fortran("bvalus",
			as.integer(sum(good)),
			knots,
			coef,
			as.integer(nknots),
			xs[good],
			s = double(sum(good)),
			as.integer(0),
                                 PACKAGE="gam")
		y[good] <- junk$s
	}
	if(any(!good)) {
		end.fit <- .Fortran("bvalus",
			as.integer(2),
			knots,
			coef,
			as.integer(nknots),
			as.double(c(0, 1)),
			s = double(2),
			as.integer(0),
                        PACKAGE="gam")$s
		end.slopes <- .Fortran("bvalus",
			as.integer(2),
			knots,
			coef,
			as.integer(nknots),
			as.double(c(0, 1)),
			s = double(2),
			as.integer(1),
                        PACKAGE="gam")$s
		if(any(bad.left))
			y[bad.left] <- end.fit[1] + end.slopes[1] * (xs[
				bad.left])
		if(any(bad.right))
			y[bad.right] <- end.fit[2] + end.slopes[2] * (xs[
				bad.right] - 1)
	}
	pred <- x * 0
	pred[!nas] <- y
	pred
}
"gam.wlist" <-
c("s","lo")
"gplot" <-
function(x, ...)
UseMethod("gplot")
"gplot.default" <-
function(x, y, se.y = NULL, xlab = "", ylab = "", residuals = NULL, rugplot = FALSE,
	scale = 0, se = FALSE, fit = TRUE, ...)
switch(data.class(x)[1],
       AsIs = { class(x)<-NULL
                gplot.default(x , y = y, se.y = se.y, xlab = xlab,
		ylab = ylab, residuals = residuals, rugplot = rugplot, scale = 
		scale, se = se, fit = fit, ...)
              },
	logical = gplot.factor(x = factor(x), y = y, se.y = se.y, xlab = xlab,
		ylab = ylab, residuals = residuals, rugplot = rugplot, scale = 
		scale, se = se, fit = fit, ...),
	list = gplot.list(x = x, y = y, se.y = se.y, xlab = xlab, ylab = ylab,
		residuals = residuals, rugplot = rugplot, scale = scale, se = 
		se, fit = fit, ...),
	if(is.numeric(x)) gplot.numeric(x = as.vector(x), y = y, se.y = se.y,
			xlab = xlab, ylab = ylab, residuals = residuals, 
			rugplot = rugplot, scale = scale, se = se, fit = fit,
			...) else warning(paste("The \"x\" component of \"",
			ylab, "\" has class \"", paste(class(x), collapse = 
			"\", \""), "\"; no gplot() methods available", sep = ""
			)))
 
"gplot.factor" <-
function(x, y, se.y = NULL, xlab, ylab, residuals = NULL, rugplot = FALSE, scale = 
	0, se = FALSE, xlim = NULL, ylim = NULL, fit = TRUE, ...)
{
	if(length(x) != length(y))
		stop("x and y do not have the same length; possibly a consequence of an na.action"
			)
	nn <- as.numeric(table(x))
	codex <- as.numeric(x)
	ucodex <- seq(nn)[nn > 0]
	o <- match(ucodex, codex, 0)
	uy <- as.numeric(y[o])
	ylim <- range(ylim, uy)
	xlim <- range(c(0, sum(nn), xlim))
	rightx <- cumsum(nn)
	leftx <- c(0, rightx[ - length(nn)])
	ux <- ((leftx + rightx)/2)
	delta <- (rightx - leftx)/8
	jx <- runif(length(codex), (ux - delta)[codex], (ux + delta)[codex])
	nnajx <- jx[!is.na(jx)]
	if(rugplot)
		xlim <- range(c(xlim, nnajx))
	if(se && !is.null(se.y)) {
		se.upper <- uy + 2 * se.y[o]
		se.lower <- uy - 2 * se.y[o]
		ylim <- range(c(ylim, se.upper, se.lower))
	}
	if(!is.null(residuals)) {
		if(length(residuals) == length(y)) {
			residuals <- y + residuals
			ylim <- range(c(ylim, residuals))
		}
		else {
			residuals <- NULL
			warning(paste("Residuals do not match x in \"", ylab,
				"\" preplot object", sep = ""))
		}
	}
	ylim <- ylim.scale(ylim, scale)
	Levels <- levels(x)
	if(!all(nn)) {
		keep <- nn > 0
		ux <- ux[keep]
		delta <- delta[keep]
		leftx <- leftx[keep]
		rightx <- rightx[keep]
		Levels <- Levels[keep]
	}
	plot(ux, uy, ylim = ylim, xlim = xlim, xlab = "", type = "n", ylab = 
		ylab, xaxt = "n", ...)
	mtext(xlab, 1, 2)
	axis(side = 3, at = ux - delta, labels = Levels, srt = 45, tick = FALSE,
		adj = 0)
	if(fit)
		segments(leftx + delta, uy, rightx - delta, uy)
	if(!is.null(residuals))
		points(jx, residuals)
	if(rugplot)
		rug(nnajx)
	if(se) {
		segments(ux + delta, se.upper, ux - delta, se.upper)
		segments(ux + delta, se.lower, ux - delta, se.lower)
		segments(ux, se.lower, ux, se.upper, lty = 2)
	}
	invisible(diff(ylim))
}
"gplot.list" <-
function(x, y, se.y = NULL, xlab, ylab, residuals = NULL, rugplot = FALSE, scale = 
	0, se = FALSE, fit = TRUE, ...)
{
	if(length(x) != 2) {
		warning(paste("A perspective plot was requested for \"", ylab,
			"\" but the \"x\" variable has dimension other than 2",
			sep = ""))
		invisible(return(0))
	}
	names(x) <- xlab
	x <- data.matrix(data.frame(x))
	#	UseMethod("gplot")
	gplot.matrix(x, y, se.y, xlab, ylab, residuals, rugplot, scale, se,
		fit, ...)
}
"gplot.matrix" <-
  function(x, y, se.y = NULL, xlab, ylab, residuals = NULL, rugplot = FALSE, scale = 
           0, se = FALSE, fit, ...)
{
  if(ncol(x) != 2) {
    warning(paste("A perspective plot was requested for \"", ylab,
                  "\" but the \"x\" variable has dimension other than 2",
                  sep = ""))
    invisible(return(0))
  }
  bivar.dup <- function(x)
    {
      if(is.null(dx <- dim(x)) || dx[2] > 2)
        stop("x must be bivariate")
      duplicated(x[, 1] + (1i) * x[, 2])
    }
  interp.loaded<-require("akima")
  if(!interp.loaded)stop("You need to install and load the package 'akima' from the R contributed libraries")
  xname <- dimnames(x)[[2]]
  dups <- bivar.dup(x)
  xyz <- interp(x[!dups, 1], x[!dups, 2], y[!dups])
  zmin <- min(xyz$z[!is.na(xyz$z)])
  z <- ifelse(is.na(xyz$z), zmin, xyz$z)
  scale2 <- diff(range(z))
                                        # Adjust scale
  scale <- max(scale, scale2)
                                        #	persp(xyz$x, xyz$y, (z - zmin)/scale, xlab = xname[1], ylab = xname[
                                        #		2], zlab = ylab, ...)
  persp(xyz$x, xyz$y, z, xlab = xname[1], ylab = xname[2], zlab = ylab,
        ...)
  invisible(scale)
}
"gplot.numeric" <-
function(x, y, se.y = NULL, xlab, ylab, residuals = NULL, rugplot = FALSE, scale = 
	0, se = FALSE, xlim = NULL, ylim = NULL, fit = TRUE, ...)
{
	if(length(x) != length(y))
		stop("x and y do not have the same length; possibly a consequence of an na.action"
			)
### Here we check if its a simple linear term; if so, we include a point at the mean of x
        if(se &&  !is.null(se.y) && ylab==paste("partial for",xlab)){
          x=c(x,mean(x))
          y=c(y,0)
          se.y=c(se.y,0)
                  }
	ux <- unique(sort(x))
	o <- match(ux, x)
	uy <- y[o]
	xlim <- range(xlim, ux)
	ylim <- range(ylim, uy)
	if(rugplot) {
		jx <- jitter(x[!is.na(x)])
		xlim <- range(c(xlim, jx))
	}
	if(se && !is.null(se.y)) {
		se.upper <- uy + 2 * se.y[o]
		se.lower <- uy - 2 * se.y[o]
		ylim <- range(c(ylim, se.upper, se.lower))
	}
	if(!is.null(residuals)) {
		if(length(residuals) == length(y)) {
			residuals <- y + residuals
			ylim <- range(c(ylim, residuals))
		}
		else {
			residuals <- NULL
			warning(paste("Residuals do not match x in \"", ylab,
				"\" preplot object", sep = ""))
		}
	}
	ylim <- ylim.scale(ylim, scale)
	if(!is.null(residuals)) {
		plot(x, residuals, xlim = xlim, ylim = ylim, xlab = xlab, ylab
			 = ylab, ...)
		if(fit)
			lines(ux, uy)
	}
	else {
		if(fit)
			plot(ux, uy, type = "l", xlim = xlim, ylim = ylim,
				xlab = xlab, ylab = ylab, ...)
	}
	if(rugplot)
		rug(jx)
	if(se) {
		lines(ux, se.upper, lty = 3)
		lines(ux, se.lower, lty = 3)
	}
	invisible(diff(ylim))
}
labels.gam<-function(object,...){
      attr(object$terms, "term.labels")
    }
"lo" <-
  function(..., span = 0.5, degree = 1)
{
  lodummy <- function(span = 0.5, degree = 1)
    list(span = span, degree = degree)
  vars <- list(...)
  locall <- sys.call()
  chcall <- deparse(locall)
  mcall <- match.call(expand = FALSE)
  mcall$... <- NULL
  nvars <- length(vars)
  if(nvars > 1) {
    scalars <- sapply(vars, length) == 1
    ## a bit of freedom in giving the span and degree
    if(any(scalars)) {
      nvars <- nvars - sum(scalars)
      mcall <- c(mcall, as.call(vars[scalars]))
      vars <- vars[!scalars]
    }
  }
  mcall[[1]] <- as.name("lodummy")
  m <- eval(mcall)
  degree <- m$degree
  span <- m$span
  if(degree > 2.)
    stop("degrees 1 or 2 are implemented")
  if(nvars == 1) {
    xvar <- as.matrix(vars[[1]])
    xnames <- deparse(locall[[2.]])
  }
  else {
    nobs <- length(vars[[1]])
    xvar <- matrix(0., nobs, nvars)
    xnames <- character(nvars)
    for(i in seq(nvars)) {
      tt <- vars[[i]]
      if(!is.null(dd <- dim(tt)) && dd[2.] > 1)
        stop("either call lo with a matrix argument, or else a comma separated list x1, x2"
             )
      exptt <- locall[[i + 1]]
      xnames[i] <- deparse(exptt)
      xvar[, i] <- as.numeric(tt)
    }
    dimnames(xvar) <- list(NULL, xnames)
  }
  ## for the moment we use polybasis from library(mda)
  polyx <- polylo(xvar, degree = degree)
  pd <- attr(polyx, "degree")
  opd <- order(pd)
  if(length(pd) > 1) {
    polyx <- polyx[, opd]
    p <- sum(pd == 1)
  }
  else p <- 1
  nobs <- dim(polyx)[1]
  nas <- is.na(polyx[, 1:p])
  if(any(nas)) {
    if(p > 1)
      nas <- nas %*% array(1, c(p, 1))
    attr(polyx, "NAs") <- seq(nobs)[nas > 0.]
  }
##  if(span * nobs < 1)
##    stop(paste("span is too small; the minimum is 1/n =", format(
##                                                                 round(1/nobs, 4.))))
  real.call <- substitute(gam.lo(data[[chcall]], z, w, span = span, 
                                 degree = degree, ncols = p), list(span = span, degree = degree,
                                                    chcall = chcall, p = p))
  atts <- c(attributes(polyx), list(span = span, degree = degree, ncols = 
                                   p, call = real.call))
  attributes(polyx) <- atts
  class(polyx) <- c("smooth", "matrix")
  polyx
}
"lo.wam" <-
function(x, y, w, s, which, smooth.frame, maxit = 30, tol = 1e-7, trace = FALSE,
	se = TRUE, ...)
{
	if(is.data.frame(smooth.frame)) {
		first <- TRUE
		# first call to wam; set up some things
		#on first entry, smooth.frame is a data frame with elements the terms to be
		#smoothed in which
		data <- smooth.frame[, names(which), drop = FALSE]
		smooth.frame <- gam.match(data)
		dx <- as.integer(dim(x))
		oldClass(data) <- NULL
		atts <- lapply(data, attributes)
		span <- sapply(atts, "[[", "span")
		degree <- sapply(atts, "[[", "degree")
		nvars <- sapply(atts, "[[", "ncols")
		ndim <- sapply(atts, "[[", "dim")[2.,  ]
		npetc <- as.integer(c(dx, length(which), 0., maxit, 0.))
		nef <- smooth.frame$nef
		nvmax <- 200. + 300. * (1. - 1./log(apply(cbind(nef - 200.,
			3.), 1., max)))
		nspar <- (nef * span + 1.)
		liv <- 50. + (2.^nvars + 4.) * nvmax + 2. * nef
		lv <- 50. + (3. * nvars + 3.) * nvmax + nef + (ifelse(degree ==
			2., ((nvars + 2.) * (nvars + 1.))/2., nvars + 1.) +
			2.) * nspar
		LL <- nspar * nvmax
		liv <- liv + LL
		lv <- lv + (nvars + 1.) * LL
		which <- sapply(which, "[", 1.)
		wddnfl <- cbind(unlist(which), nvars, ndim, degree, nef, liv,
			lv, nvmax)
		storage.mode(wddnfl) <- "integer"
		spatol <- as.double(c(span, tol))
		nwork <- 9. * dx[1.] + sum(nef * (nvars + ndim + 4.) + 5. +
			3. * ndim)
		liv <- sum(liv)
		lv <- sum(lv)
		smooth.frame <- c(smooth.frame, list(npetc = npetc, wddnfl = 
			wddnfl, spatol = spatol,niwork=2*sum(nvars), nwork = nwork, liv = liv,
			lv = lv))
	}
	else first <- FALSE
	storage.mode(y) <- "double"
	storage.mode(w) <- "double"
	n <- smooth.frame$npetc[1.]
	p <- smooth.frame$npetc[2.]
	q <- smooth.frame$npetc[3.]
	fit <- .Fortran("baklo",
		x,
		y = y,
		w = w,
		npetc = smooth.frame$npetc,
		smooth.frame$wddnfl,
		smooth.frame$spatol,
		smooth.frame$o,
		etal = double(n),
		s = s,
		eta = double(n),
		beta = double(p),
		var = s,
		df = double(q),
		qr = x,
		qraux = double(p),
		qpivot = as.integer(1.:p),
                effects=double(n),
		integer(smooth.frame$liv),
		double(smooth.frame$lv),
                integer(smooth.frame$niwork),
		double(smooth.frame$nwork),
                        PACKAGE="gam")

	nit <- fit$npetc[4.]
	qrank <- fit$npetc[6.]
	if((nit == maxit) & maxit > 1.)
		warning(paste("lo.wam convergence not obtained in ", maxit,
			" iterations"))
	names(fit$df) <- dimnames(s)[[2]]
	names(fit$beta) <- labels(x)[[2]]
                qrx <- structure(list(qr = fit$qr,qraux = fit$qraux,
                     rank = qrank, pivot = fit$qpivot,tol=1e-7),class="qr")
        effects<-fit$effects
        r1 <- seq(len = qrx$rank)
        dn <- colnames(x)
        if (is.null(dn)) 
          dn <- paste("x", 1:p, sep = "")
        names(effects) <- c(dn[qrx$pivot[r1]], rep.int("", n - qrx$rank))
	rl <- list(coefficients = fit$beta, residuals = fit$y - fit$eta, 
                   fitted.values = fit$eta,
                   effects=effects, weights=w, rank=qrank,
                   assign=attr(x,"assign"),
                   qr=qrx,
                   smooth = fit$s,
                   nl.df = fit$df
                   )
	rl$df.residual <- n - qrank - sum(rl$nl.df) - sum(fit$w == 0.)
	if(se)
		rl <- c(rl, list(var = fit$var))
	c(list(smooth.frame = smooth.frame), rl)
}
"na.gam.replace" <-
function(frame)
{
	vars <- names(frame)
	if(!is.null(resp <- attr(attr(frame, "terms"), "response"))) {
		vars <- vars[ - resp]
		x <- frame[[resp]]
		pos <- is.na(x)
		if(any(pos)) {
			frame <- frame[!pos,  , drop = FALSE]
			warning(paste(sum(pos), 
				"observations omitted due to missing values in the response"
				))
		}
	}
	for(j in vars) {
		x <- frame[[j]]
		pos <- is.na(x)
		if(any(pos)) {
			if(length(levels(x))) {
				xx <- as.character(x)
				xx[pos] <- "NA"
				x <- factor(xx, exclude = NULL)
			}
			else if(is.matrix(x)) {
				ats <- attributes(x)
				w <- !pos
				x[pos] <- 0
				n <- nrow(x)
				TT <- array(1, c(1, n))
				xbar <- (TT %*% x)/(TT %*% w)
				xbar <- t(TT) %*% xbar
				x[pos] <- xbar[pos]
				attributes(x) <- ats
			}
			else {
				ats <- attributes(x)
				x[pos] <- mean(x[!pos])
				attributes(x) <- ats
			}
			frame[[j]] <- x
		}
	}
	frame
}
"newdata.predict.gam" <-
  function(object, newdata, type = c("link", "response", "terms"), dispersion=NULL, se.fit = FALSE, na.action=na.pass,terms=labels(object), ...)
{
  out.attrs <- attr(newdata, "out.attrs")
  is.gam<-inherits(object, "gam") && !is.null(object$smooth)
 if(is.gam) {
   if(se.fit){
     se.fit<-FALSE
     warning("No standard errors (currently) for gam predictions with newdata")
   }
   ##First get the linear predictions
   type <- match.arg(type)
   local.type<-type
   if(type=="response")local.type<-"link"
   pred<-predict.glm(object,newdata,type=local.type,dispersion=dispersion,se.fit=FALSE,terms=terms)
   ##Build up the smooth.frame for the new data
   tt <- terms(object)
    Terms <- delete.response(tt)
    smooth.frame <- model.frame(Terms, newdata, na.action = na.action, 
                     xlev = object$xlevels)
   nrows<-nrow(smooth.frame)
   old.smooth<-object$smooth
   data<-object$smooth.frame # this was the old smooth frame
   smooth.labels<-names(data)
   n.smooths<-length(smooth.labels)
   if (!is.null(cl <- attr(Terms, "dataClasses"))) 
      .checkMFClasses(cl, smooth.frame)
    out.attrs <- attr(newdata, "out.attrs")
  

   w <- object$weights
   pred.s <- array(0, c(nrows, n.smooths), list(row.names(smooth.frame), 
                                                 smooth.labels))
   smooth.wanted <- smooth.labels[match(smooth.labels, terms,
                                         0) > 0]
   pred.s<-pred.s[,smooth.wanted,drop=FALSE]
    residuals <- object$residuals
    for(TT in smooth.wanted) {
      Call <- attr(data[[TT]], "call")
      Call$xeval <- substitute(smooth.frame[[TT]], list(TT = TT))
      z <- residuals + object$smooth[, TT]
       pred.s[, TT] <- eval(Call)
    }
    if(type == "terms")
      pred[, smooth.wanted] <- pred[, smooth.wanted] + pred.s[
                                                              , smooth.wanted]
    else pred <- pred + rowSums(pred.s)
   if(type == "response") {
     famob <- family(object)
     pred <- famob$linkinv(pred)
   }
  }
  else {
    pred<-predict.glm(object,newdata,type=type,dispersion=dispersion,se.fit=se.fit,terms=terms)
  }
  if(type != "terms" && !is.null(out.attrs)) {
    if(!is.null(out.attrs)) {
      if(se.fit) {
        attributes(pred$fit) <- out.attrs
        attributes(pred$se.fit) <- out.attrs
      }
      else attributes(pred) <- out.attrs
    }
  }
pred
}
"plot.gam" <-
  function(x,  residuals = NULL, rugplot = TRUE, se = FALSE, scale = 0, ask = FALSE,
terms=labels.gam(x), ...)
{
  
  if(!is.null(x$na.action))
    x$na.action <- NULL
  preplot.object <- x$preplot
  if(is.null(preplot.object))
    preplot.object <- preplot.gam(x,terms=terms)
  x$preplot <- preplot.object
  Residuals <- resid(x)
  if(!is.null(residuals)) {
    if(length(residuals) == 1)
      if(residuals)
        residuals <- Residuals
      else residuals <- NULL
    else Residuals <- residuals
  }
  if(!ask) {
    plot.preplot.gam(preplot.object, residuals = residuals, rugplot
                     = rugplot, scale = scale, se = se, fit = TRUE, ...)
    invisible(x)
  }
  else{
    nterms <- names(preplot.object)
    tterms <- substring(nterms, 1, 40)
                                        #truncate long names
    residualsmenu <- if(!is.null(residuals)) "residuals off" else 
    "residuals on"
    rugmenu <- if(rugplot) "rug off" else "rug on"
    semenu <- if(se) "se off" else "se on"
    scalemenu <- paste("scale (", round(scale, 1), ")", sep = "")
    scales <- numeric()
    tmenu <- c(paste("plot:", tterms), "plot all terms", residualsmenu,
               rugmenu, semenu, scalemenu)
    tnames <- character()
    pick <- 1
    while(pick > 0 && pick <= length(tmenu)) {
      pick <- menu(tmenu, title = 
                   "Make a plot selection (or 0 to exit):\n")
      if(pick > 0 && pick <= length(nterms)) {
        tscale <- plot.preplot.gam(preplot.object[[pick]],
                                   residuals = residuals, rugplot = rugplot, scale
                                   = scale, se = se, fit = TRUE, ...)
        names(tscale) <- nterms[pick]
        scales <- c(scales, tscale)
        cat("Plots performed:\n ")
        print(scales)
      }
      else switch(pick - length(nterms),
                  {
                    scales <- plot.preplot.gam(
                                               preplot.object, residuals = 
                                               residuals, rugplot = rugplot,
                                               scale = scale, se = se, fit = 
                                               TRUE, ...)
                    print(scales)
                  }
                  ,
                  {
                    residuals <- if(is.null(residuals)) 
                      Residuals else NULL
                    residualsmenu <- if(!is.null(residuals)
                                        ) "residuals off" else 
                    "residuals on"
                  }
                  ,
                  {
                    rugplot <- !rugplot
                    rugmenu <- if(rugplot) "rug off" else 
                    "rug on"
                  }
                  ,
                  {
                    se <- !se
                    semenu <- if(se) "se off" else "se on"
                  }
                  ,
                  {
                    cat("Type in a new scale\n")
                    scale <- eval(parse(n=1))
                    scalemenu <- paste("scale (", round(
                                                        scale, 1), ")", sep = "")
                  }
                  ,
                  invisible(return(x)))
      tmenu <- c(paste("plot:", tterms), "plot all terms", 
                 residualsmenu, rugmenu, semenu, scalemenu)
    }
    invisible(x)
  }
}
"plot.preplot.gam" <-
function(x, y = NULL, residuals = NULL, rugplot = TRUE, se = FALSE, scale = 0, fit = TRUE,
	...)
{
	listof <- inherits(x[[1]], "preplot.gam")
	if(listof) {
		TT <- names(x)
		scales <- rep(0, length(TT))
		names(scales) <- TT
		for(i in TT)
			scales[i] <- plot.preplot.gam(x[[i]], y = NULL, 
				residuals, rugplot, se, scale, fit, ...)
		#			scales[i] <- UseMethod("plot",x[[i]])
		invisible(scales)
	}
	else {
		dummy <- function(residuals = NULL, rugplot = TRUE, se = FALSE, scale
			 = 0, fit = TRUE, ...)
		c(list(residuals = residuals, rugplot = rugplot, se = se, scale
			 = scale, fit = fit), list(...))
		d <- dummy(residuals, rugplot, se, scale, fit, ...)
		uniq.comps <- unique(c(names(x), names(d)))
		Call <- c(as.name("gplot"), c(d, x)[uniq.comps])
		mode(Call) <- "call"
		invisible(eval(Call))
	}
}
"polylo" <-
  function (x, degree = 1, monomial = FALSE) 
{
  if (degree >= 4) 
    warning("This is not a smart polynomial routine. You may get numerical problems with a degree of 4 or more")
  x <- as.matrix(x)
  dn <- dimnames(x)
  dd <- dim(x)
  np <- dd[2]
  ad<-rep(1,ncol(x))
  if (np == 1) 
    monomial <- TRUE
  if (degree > 1) {
    if (monomial) {
      ad<-seq(degree)
      px <- x
      cc <- sapply(split(paste(diag(np)), rep(seq(np), 
                                              rep(np, np))), paste, collapse = "")
      tx <- x
      for (i in 2:degree) {
        px <- px * tx
        x <- cbind(x, px)
        cc <- c(cc, sapply(split(paste(diag(np) * i), 
                                 rep(seq(np), rep(np, np))), paste, collapse = ""))
      }
      
    }
    else {
      matarray <- array(x, c(dd, degree))
      for (i in 2:degree) matarray[, , i] <- x^i
      matarray <- aperm(matarray, c(1, 3, 2))
      x <- matarray[, , np,drop=FALSE]
      ad0 <- seq(degree)
      ad <- ad0
      ncol.mat0 <- degree
      ncol.x <- degree
      d0 <- paste(ad0)
      cc <- d0
      for (ii in seq(np - 1, 1)) {
        index0 <- rep(seq(ncol.mat0), ncol.x)
        index <- rep(seq(ncol.x), rep(ncol.mat0, ncol.x))
        newad <- ad0[index0] + ad[index]
        retain <- newad <= degree
        mat0 <- matarray[, , ii,drop=FALSE]
        browser()
        if (any(retain)) 
          newmat <- mat0[, index0[retain],, drop = FALSE] * 
            x[, index[retain], ,drop = FALSE]
        else newmat <- NULL
        ddn <- paste(d0[index0[retain]], cc[index[retain]], 
                     sep = "")
        zeros <- paste(rep(0, nchar(cc[1])), collapse = "")
        cc <- paste(0, cc, sep = "")
        d00 <- paste(d0, zeros, sep = "")
        x <- cbind(mat0, x, newmat)
        cc <- c(d00, cc, ddn)
        ad <- c(ad0, ad, newad[retain])
        ncol.x <- length(ad)
      }
    }
    if (!is.null(dn)) 
      dn[[2]] <- cc
    else dn <- list(NULL, cc)
    dimnames(x) <- dn
  }
  attr(x,"degree")<-ad
  x
}
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
"preplot.gam" <-
  function(object, newdata, terms = labels.gam(object),...)
{
  ## this labels.gam above is because there does not seem to be a label method for glms
  Terms <- object$terms
  a <- attributes(Terms)
  Call <- object$call
  all.terms <- labels(Terms)
  xvars <- parse(text=all.terms)
  names(xvars) <- all.terms
  terms <- sapply(terms,match.arg, all.terms)
  Interactions <- a$order > 1
 if(any(Interactions)) {
    all.terms <- all.terms[!Interactions]
    TM <- match(terms, all.terms, 0)
    if(!all(TM)) {
      terms <- terms[TM > 0]
      warning("No terms saved for \"a:b\" style interaction terms"
              )
    }
  }
   xvars <- xvars[terms]
  xnames <- as.list(terms)
  names(xnames) <- terms
  modes <- sapply(xvars, mode)
   for(term in terms[modes != "name"]) {
    evars <- all.names(xvars[term], functions = FALSE, unique = TRUE)
    if(!length(evars))
      next
    xnames[[term]] <- evars
    evars <- parse(text = evars)
    if(length(evars) == 1)
      evars <- evars[[1]]
    else {
      evars <- c(as.name("list"), evars)
      mode(evars) <- "call"
    }
    xvars[[term]] <- evars
  }
  xvars <- c(as.name("list"), xvars)
  mode(xvars) <- "call"

  if(!missing(newdata))
    xvars <- eval(xvars, newdata)
  else {
    if(!is.null(Call$subset) | !is.null(Call$na.action) | !is.null(
                                                                   options("na.action")[[1]])) {
      Rownames <- names(object$fitted)

      if(!(Rl <- length(Rownames)))
        stop("need to have names for fitted.values when call has a subset or na.action argument"
             )
      form<-paste("~",unlist(xnames),collapse="+")
      Mcall <- c(as.name("model.frame"), list(formula = 
                                              terms(as.formula(form)),
                                              subset = Rownames, na.action = function(x)
                                              x))
      mode(Mcall) <- "call"
      Mcall$data <- Call$data
      xvars <- eval(xvars, eval(Mcall))
    }
    else {
      ecall <- substitute(eval(expression(xvars)))
      ecall$local <- Call$data
      xvars <- eval(ecall)
    }
  }
  if(missing(newdata))
    pred <- predict(object, type = "terms", terms = terms,
			se.fit = TRUE)
  else pred <- predict(object, newdata, type = "terms", terms = terms,
                           se.fit = TRUE)
  fits <- pred$fit
  if(is.null(fits)) {
    fits <- pred
    se.fits <- NULL
  }
  else se.fits <- pred$se.fit
  gamplot <- xnames
  for(term in terms) {
    x <- xvars[[term]]
    ## oldClass(x) <- unique(c(oldClass(x), data.class(unclass(x))))
    xlab <- xnames[[term]]
    ## Fix ylab for linear terms:
    ylab <- if(length(xlab) == 1 && term == xlab) paste(
                                      "partial for", term) else term
    TT <- list(x = x, y = fits[, term], se.y = if(is.null(se.fits)
                                          ) NULL else se.fits[, term], xlab = xlab, ylab = ylab)
    oldClass(TT) <- "preplot.gam"
    gamplot[[term]] <- TT
  }
  oldClass(gamplot) <- "preplot.gam"
  gamplot
}
"print.gam" <-
  function(x, digits = 5, ...)
{
  if(!is.null(cl <- x$call)) {
    cat("Call:\n")
    dput(cl)
  }
  n <- x$df.null
  if(is.null(df.resid <- x$df.resid))
    df.resid <- n - sum(!is.na(x$coef)) - sum(x$nl.df)
  cat("\nDegrees of Freedom:", n, "total;", format(round(df.resid, digits
                                                         )), "Residual\n")
  if(!is.null(x$na.action))
    cat(naprint(x$na.action), "\n")
  cat("Residual Deviance:", format(round(x$deviance, digits)), "\n")
  invisible(x)
}
"print.gamex" <-
  function(x,...)
  {
    print(x$coefficients)
    invisible()
  }

"print.summary.gam" <-
  function(x,  digits = max(3, getOption("digits") - 3), quote = TRUE, prefix = "", ...)
{
  cat("\nCall: ")
  dput(x$call)
  dresid <- x$deviance.resid
  n <- length(dresid)
  rdf <- x$df[2]
  if(rdf > 5) {
    cat("Deviance Residuals:\n")
    rq <- quantile(as.vector(dresid))
    names(rq) <- c("Min", "1Q", "Median", "3Q", "Max")
    print(rq, digits = digits)
  }
  else if(rdf > 0) {
    cat("Deviance Residuals:\n")
    print(dresid, digits = digits)
  }
  cat(paste("\n(Dispersion Parameter for ", names(x$dispersion), 
            " family taken to be ", format(round(x$dispersion, digits)),
            ")\n",sep=""))
  int <- attr(x$terms, "intercept")
  cat("\n    Null Deviance:", format(round(x$null.deviance, digits)),
      "on", n - int, "degrees of freedom")
  cat("\nResidual Deviance:", format(round(x$deviance, digits)), "on",
      format(round(rdf, digits)), "degrees of freedom")
  cat("\nAIC:", format(round(x$aic, digits)),"\n")
  if(!is.null(x$na.action))
    cat(naprint(x$na.action), "\n")
  cat("\nNumber of Local Scoring Iterations:", format(trunc(x$iter)),
      "\n")
  print(x$anova)
}
"random" <-
  function(xvar, df = NULL, sigma = 0.)
{
  scall <- deparse(sys.call())
  if(!inherits(xvar, "factor"))
    stop("random() expects a factor or category as its first argument"
         )
  xvar <- C(xvar, rep(0., length(levels(xvar))), 1.)
  attr(xvar, "call") <- substitute(gam.random(data[[scall]], z, w, df = 
                                              df, sigma))
  oldClass(xvar) <- c("smooth", oldClass(xvar))
  xvar
}
"s" <-
  function(x, df = 4, spar = 1)
{
  scall <- deparse(sys.call())
  if(missing(df)){
    if(!missing(spar))df<-0
  }
    
  if(ncol(as.matrix(x)) > 1)
    stop(paste(
               "The default smoother is bivariate; you gave a matrix as an argument in ",
               scall, "\n"))
  if(!is.null(levels(x))) {
    if(inherits(x, "ordered"))
      x <- as.numeric(x)
    else stop("unordered factors cannot be used as smoothing variables"
              )
  }
  attr(x, "spar") <- spar
  attr(x, "df") <- df
  real.call <- substitute(gam.s(data[[scall]], z, w, spar = spar, df = df
                                ))
  attr(x, "call") <- real.call
  attr(x, "class") <- "smooth"
  a <- is.na(x)
  if(any(a))
    attr(x, "NAs") <- seq(along = x)[a]
  x
}
"s.wam" <-
function(x, y, w, s, which, smooth.frame, maxit = 30, tol = 1e-7, trace = FALSE,
	se = TRUE, ...)
{
	if(is.data.frame(smooth.frame)) {
		first <- TRUE
		# first call to wam; set up some things
		#on first entry, smooth.frame is a data frame with elements the terms to be
		#smoothed in which
		data <- smooth.frame[, names(which), drop = FALSE]
		smooth.frame <- gam.match(data)
		dx <- as.integer(dim(x))
		smooth.frame$n <- dx[1]
		smooth.frame$p <- dx[2]
		oldClass(data) <- NULL
		smooth.frame$spar <- unlist(lapply(data, attr, "spar"))
		smooth.frame$df <- unlist(lapply(data, attr, "df"))
	}
	else first <- FALSE
	storage.mode(tol) <- "double"
	storage.mode(maxit) <- "integer"
	which <- unlist(which)
	storage.mode(which) <- "integer"
	storage.mode(y) <- "double"
	storage.mode(w) <- "double"
	p <- smooth.frame$p
	n <- smooth.frame$n
	fit <- .Fortran("bakfit",
		x,
		npetc = as.integer(c(n, p, length(which), se, 0, maxit, 0)),

			y = y,
		w = w,
		which,
		spar = as.double(smooth.frame$spar),
		df = as.double(smooth.frame$df),
		as.integer(smooth.frame$o),
		as.integer(smooth.frame$nef),
		etal = double(n),
		s = s,
		eta = double(n),
		beta = double(p),
		var = s,
		tol,
		qr = x,
		qraux = double(p),
		qpivot = as.integer(1:p),
                effects=double(n),        
		double((10 + 2 * 4 + 5) * (max(smooth.frame$nef) + 2) + 15 *
			n + 15 + length(which)),
                        PACKAGE="gam")
	nit <- fit$npetc[5]
	qrank <- fit$npetc[7]
	if((nit == maxit) & maxit > 1)
		warning(paste("s.wam convergence not obtained in ", maxit,
			" iterations"))
	if(first) {
		smooth.frame$spar <- fit$spar
		first <- FALSE
	}
	names(fit$df) <- dimnames(s)[[2]]
	names(fit$beta) <- labels(x)[[2]]
        qrx <- structure(list(qr = fit$qr,qraux = fit$qraux,
                     rank = qrank, pivot = fit$qpivot,tol=1e-7),class="qr")
        effects<-fit$effects    #qr.qty(qrx,fit$y)
        r1 <- seq(len = qrx$rank)
        dn <- colnames(x)
        if (is.null(dn)) 
          dn <- paste("x", 1:p, sep = "")
        names(effects) <- c(dn[qrx$pivot[r1]], rep.int("", n - qrx$rank))

 	rl <- list(coefficients = fit$beta, residuals = fit$y - fit$eta, 
                   fitted.values = fit$eta,
                   effects=effects, weights=w, rank=qrank,
                   assign=attr(x,"assign"),
                   qr=qrx,
                   smooth = fit$s,
                   nl.df = fit$df - 1
                   )
	rl$df.residual <- n - qrank - sum(rl$nl.df) - sum(fit$w == 0.)
	if(se)
		rl <- c(rl, list(var = fit$var))
	c(list(smooth.frame = smooth.frame), rl)
}
"step.gam" <-
  function(object, scope, scale, direction = c("both", "backward", "forward"),
           trace = TRUE, keep = NULL, steps = 1000, ...)
{
  scope.char <- function(formula){
   tt<-terms(formula)
  tl<-attr(tt,"term.labels")
  if(attr(tt,"intercept"))c("1",tl)else tl
 }
  re.arrange <- function(keep)
    {
      namr <- names(k1 <- keep[[1]])
      namc <- names(keep)
      nc <- length(keep)
      nr <- length(k1)
      array(unlist(keep, recursive = FALSE), c(nr, nc), list(namr, namc))
    }
  untangle.scope <- function(terms, regimens)
    {
      a <- attributes(terms)
      response <- deparse(a$variables[[2]])
      term.labels <- a$term.labels
      if(!is.null(a$offset)) {
        off1 <- deparse(a$variables[[a$offset]])
      }
      nt <- length(regimens)
      select <- integer(nt)
      for(i in seq(nt)) {
        j <- match(regimens[[i]], term.labels, 0)
        if(any(j)) {
          if(sum(j > 0) > 1)
            stop(paste("The elements of a regimen",
                       i, 
                       "appear more than once in the initial model",
                       sep = " "))
          select[i] <- seq(j)[j > 0]
          term.labels <- term.labels[ - sum(j)]
        }
        else {
          if(!(j <- match("1", regimens[[i]], 0)))
            stop(paste("regimen", i, 
                       "does not appear in the initial model",
                       sep = " "))
          select[i] <- j
        }
      }
      if(length(term.labels))
        term.labels <- paste(term.labels, "+")
      if(!is.null(a$offset))
        term.labels <- paste(off1, term.labels, sep = " + ")
      return(list(response = paste(response, term.labels, sep = " ~ "
                    ), select = select))
    }
  make.step <- function(models, fit, scale, object)
    {
      chfrom <- sapply(models, "[[", "from")
      chfrom[chfrom == "1"] <- ""
      chto <- sapply(models, "[[", "to")
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
      aod <- data.frame(From = chfrom, To = chto, Df = ddf, Deviance
                        = ddev, "Resid. Df" = df, "Resid. Dev" = dev, AIC = 
			AIC, check.names = FALSE)
      fit$anova <- as.anova(aod, heading)
      fit
    }
  direction <- match.arg(direction)
  if(missing(scope))
    stop("you must supply a scope argument to step.gam(); the gam.scope() function might be useful"
         )
  if(!is.character(scope[[1]]))
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
  visited <- array(FALSE, term.lengths)
  visited[array(items, c(1, n.items))] <- TRUE
  if(!is.null(keep)) {
    keep.list <- vector("list", length(visited))
    nv <- 1
  }
  models <- vector("list", length(visited))
  nm <- 2
  form.vector <- character(n.items)
  for(i in seq(n.items))
    form.vector[i] <- scope[[i]][items[i]]
  form <- deparse(object$formula)
  if(trace)
    cat("Start: ", form)
  fit <- object
  n <- length(fit$fitted)
  if(missing(scale)) {
    famname <- family$family["name"]
    scale <- switch(famname,
                    Poisson = 1,
                    Binomial = 1,
                    deviance.lm(fit)/fit$df.resid)
  }
  else if(scale == 0)
    scale <- deviance.lm(fit)/fit$df.resid
###  bAIC <- deviance(fit) + 2 * (n - fit$df.resid) * scale
  bAIC<-fit$aic
  if(trace)
    cat("; AIC=", format(round(bAIC, 4)), "\n")
  models[[1]] <- list(deviance = deviance(fit), df.resid = fit$df.resid,
                      AIC = bAIC, from = "", to = "")
  if(!is.null(keep)) {
    keep.list[[nv]] <- keep(fit, bAIC)
    nv <- nv + 1
  }
  AIC <- bAIC + 1
  while(bAIC < AIC & steps > 0) {
    steps <- steps - 1
    AIC <- bAIC
    bitems <- items
    bfit <- fit
    for(i in seq(n.items)) {
      if(backward) {
        trial <- items
                                        # try go down a level
        trial[i] <- trial[i] - 1
        if(trial[i] > 0 && !visited[array(trial, c(
                                                   1, n.items))]) {
          visited[array(trial, c(1, n.items))] <-
            TRUE
          tform.vector <- form.vector
          tform.vector[i] <- scope[[i]][trial[
                                              i]]
          form <- paste(form.y, paste(
                                      tform.vector, collapse = " + ")
                        )
          if(trace)
            cat("Trial: ", form)
          tfit <- update(object, eval(parse(
                                            text = form)), trace = FALSE, ...)
          tAIC<-tfit$aic
###          tAIC <- deviance(tfit) + 2 * (n - tfit$
###                                        df.resid) * scale
          if(!is.null(keep)) {
            keep.list[[nv]] <- keep(tfit,
                                    tAIC)
            nv <- nv + 1
          }
          if(tAIC < bAIC) {
            bAIC <- tAIC
            bitems <- trial
            bfit <- tfit
            bform.vector <- tform.vector
            bfrom <- form.vector[i]
            bto <- tform.vector[i]
          }
          if(trace)
            cat("; AIC=", format(round(
                                       tAIC, 4)), "\n")
        }
      }
      if(forward) {
        trial <- items
        trial[i] <- trial[i] + 1
        if(trial[i] <= term.lengths[i] && !visited[
                  array(trial, c(1, n.items))]) {
          visited[array(trial, c(1, n.items))] <-
            TRUE
          tform.vector <- form.vector
          tform.vector[i] <- scope[[i]][trial[
                                              i]]
          form <- paste(form.y, paste(
                                      tform.vector, collapse = " + ")
                        )
          if(trace)
            cat("Trial: ", form)
          tfit <- update(object, eval(parse(
                                            text = form)), trace = FALSE, ...)
          tAIC<-tfit$aic
###          tAIC <- deviance(tfit) + 2 * (n - tfit$
###                                        df.resid) * scale
          if(!is.null(keep)) {
            keep.list[[nv]] <- keep(tfit,
                                    tAIC)
            nv <- nv + 1
          }
          if(tAIC < bAIC) {
            bAIC <- tAIC
            bitems <- trial
            bfit <- tfit
            bform.vector <- tform.vector
            bfrom <- form.vector[i]
            bto <- tform.vector[i]
          }
          if(trace)
            cat("; AIC=", format(round(
                                       tAIC, 4)), "\n")
        }
      }
    }
    if(bAIC >= AIC | steps == 0) {
      if(!is.null(keep))
        fit$keep <- re.arrange(keep.list[seq(nv - 1)])
      return(make.step(models[seq(nm - 1)], fit, scale, 
                       object))
    }
    else {
      if(trace)
        cat("Step : ", deparse(bfit$formula), "; AIC=",
            format(round(bAIC, 4)), "\n\n")
      items <- bitems
      models[[nm]] <- list(deviance = deviance(bfit), 
                           df.resid = bfit$df.resid, AIC = bAIC, from = 
                           bfrom, to = bto)
      nm <- nm + 1
      fit <- bfit
      form.vector <- bform.vector
    }
  }
}
"summary.gam" <-
  function(object, dispersion = NULL,...)
{
  save.na.action <- object$na.action
  object$na.action <- NULL
  fun <- function(assign, coeff)
    sum(!is.na(coeff[assign]))
  wt <- object$weights
  coef <- object$coef
  dresid <- residuals(object, "deviance")
  resid <- object$residuals
  n <- length(resid)
  s <- object$s
  nl.chisq <- object$nl.chisq
  assg <- object$assign
  if(is.null(assg))
    assg <- attributes(object$terms)$assign
  df<-rep(1,length(assg))
  df[is.na(object$coef)]<-0
  df<-tapply(df,assg,sum)
  dfnames<-attr(object$terms,"term.labels")
  if(attr(object$terms,"intercept"))dfnames<-c("(Intercept)",dfnames)
  names(df)<-dfnames
  df<-unlist(df)
  nldf <- object$nl.df
  n <- length(object$residuals)
  if(is.null(rdf <- object$df.resid)) {
    rdf <- n - sum(df)
    if(!is.null(nldf))
      rdf <- rdf - sum(nldf)
  }
  if(!is.null(wt)) {
    wt <- wt^0.5
    resid <- resid * wt
    excl <- wt == 0
    if(any(excl)) {
      warning(paste(sum(excl), 
                    "rows with zero weights not counted"))
      resid <- resid[!excl]
      dresid <- dresid[!excl]
      if(is.null(object$df.residual))
        rdf <- rdf - sum(excl)
    }
  }
  if(rdf > 0)
    phihat <- sum(resid^2)/rdf
  else {
    phihat <- Inf
    warning("Residual degrees of freedom are negative or zero.  This occurs when the sum of the parametric and nonparametric degrees of freedom exceeds the number of observations.  The model is probably too complex for the amount of data available."
            )
  }
  famname <- object$family[["family"]]
  if(is.null(famname))
    famname <- "gaussian"
  chiorf <- TRUE
  if(!is.null(dispersion) && dispersion == 0)
    dispersion <- phihat
  if(is.null(dispersion))
    dispersion <- switch(famname,
                         poisson = 1,
                         binomial = 1,
                         {
                           chiorf <- FALSE
                           phihat
                         }
                         )
  names(dispersion) <- famname
  if(length(df)) {
    aod <- as.matrix(round(df, 1))
    dimnames(aod) <- list(names(df), "Df")
    if(!is.null(nl.chisq)) {
      aod <- cbind(aod, NA, NA, NA)
      nl.chisq <- nl.chisq/dispersion
      snames <- names(nldf)
      aod[snames, 2] <- round(nldf, 1)
      aod[snames, 3] <- if(!chiorf) nl.chisq/nldf else 
      nl.chisq
      aod[snames, 4] <- if(chiorf) 1 - pchisq(nl.chisq, nldf)
      else if(rdf > 0)
        1 - pf(nl.chisq/nldf, nldf, rdf)
      else NA
      rnames <- c("Df", "Npar Df", "Npar Chisq", "P(Chi)")
      if(!chiorf)
        rnames[3:4] <- c("Npar F", "Pr(F)")
      dimnames(aod) <- list(names(df), rnames)
      heading <- if(chiorf) 
        "\nDF for Terms and Chi-squares for Nonparametric Effects\n"
      else "\nDF for Terms and F-values for Nonparametric Effects\n"
    }
    else heading <- "DF for Terms\n\n"
    aod <- as.anova(data.frame(aod, check.names = FALSE), heading)
  }
  else aod <- NULL
  structure(list(call = object$call, terms = object$terms, anova = aod,
                 dispersion = dispersion, df = c(sum(df) + sum(nldf), rdf),
                 deviance.resid = dresid, deviance = deviance(object), 
                 null.deviance = object$null.deviance, aic=object$aic,iter = object$iter, 
                 na.action = save.na.action), class = "summary.gam")
}
"ylim.scale" <-
function(ylim, scale = 0.)
{
	scale2 <- diff(ylim)
	if(scale2 < scale)
		rep(mean(ylim), 2.) + ((ylim - mean(ylim)) * scale)/scale2
	else ylim
}
