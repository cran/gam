"anova.gamlist" <-
function(object, ..., test = c("none", "Chisq", "F", "Cp")){
  test=match.arg(test)
  anova.glmlist(object, test = test)
}
