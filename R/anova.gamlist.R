#' @method anova Gamlist
#' @export
#' @export anova.Gamlist
"anova.Gamlist" <-
function(object, ..., test = c("none", "Chisq", "F", "Cp")){
  test=match.arg(test)
  class(object)="glmlist"
  anova(object, test = test)
}
