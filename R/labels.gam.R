#' @method labels Gam
#' @export
#' @export labels.Gam
labels.Gam<-function(object,...){
      attr(object$terms, "term.labels")
    }
