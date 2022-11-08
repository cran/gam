#' @method print Gamex
#' @export
#' @export print.Gamex
"print.Gamex" <-
  function(x,...)
  {
    print(x$coefficients)
    invisible()
  }

