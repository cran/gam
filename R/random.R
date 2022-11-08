#' @rdname gam.random
#' @export random
"random" <-
    function (f, df = NULL, lambda = 0, intercept = TRUE)
{
    scall <- deparse(sys.call())
    if (!inherits(f, "factor"))
        stop("random() expects a factor or category as its first argument")
    newf=rep(0,length(f))
    attr(newf,"values")=f
    attr(newf, "call") <- substitute(gam.random(data[[scall]], z,
                                             w, df = df, lambda = lambda, intercept = intercept))
    oldClass(newf) <- "smooth"
    newf
}
