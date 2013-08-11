"anova.gamlist" <-
function(object, ..., test = c("none", "Chisq", "F", "Cp")){
  test=match.arg(test)
  anova(structure(c(list(object), ...), class="glmlist"), test = test)
  ## anova.glmlist(object, test = test)
}
